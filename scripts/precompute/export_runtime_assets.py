from __future__ import annotations
import argparse
import json
import struct
from pathlib import Path

import h5py
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MAGIC = 0x4A525431  # JRT1
VERSION = 1


def write_header(f):
    f.write(struct.pack('<II', MAGIC, VERSION))


def write_string(f, s: str):
    b = s.encode('utf-8')
    f.write(struct.pack('<I', len(b)))
    f.write(b)


def downsample_indices(n: int, max_bins: int | None) -> np.ndarray:
    if max_bins is None or max_bins <= 0 or n <= max_bins:
        return np.arange(n, dtype=np.int64)
    return np.unique(np.round(np.linspace(0, n - 1, max_bins)).astype(np.int64))


def lambda_from_h5(h) -> np.ndarray:
    if 'wavelength_m' in h:
        return np.array(h['wavelength_m'], dtype=np.float64).reshape(-1)
    for key in ['wavenumber', 'wavenumber_cm1', 'wavenumbers_cm1']:
        if key in h:
            wn = np.array(h[key], dtype=np.float64).reshape(-1)
            wn = np.maximum(wn, 1.0e-300)
            return 1.0 / (100.0 * wn)
    raise KeyError('No wavelength_m or wavenumber axis found')


def read_opacity(path: Path):
    with h5py.File(path, 'r') as h:
        species = str(h.attrs.get('species', path.name.split('_')[0]))
        lam = lambda_from_h5(h)
        p = np.array(h.get('pressure', h.get('pressures_bar')), dtype=np.float64).reshape(-1)
        T = np.array(h.get('temperature', h.get('temperatures_K')), dtype=np.float64).reshape(-1)
        sigma = np.array(h['sigma'], dtype=np.float32)
    order = np.argsort(lam)
    lam = lam[order]
    sigma = sigma[:, :, order]
    return species, lam, p, T, sigma


def interpolate_to_master(lam_src, values_src, lam_master):
    # values_src shape (..., lambda)
    flat = values_src.reshape((-1, values_src.shape[-1]))
    out = np.empty((flat.shape[0], lam_master.size), dtype=np.float32)
    for i in range(flat.shape[0]):
        out[i, :] = np.interp(lam_master, lam_src, flat[i, :], left=0.0, right=0.0).astype(np.float32)
    return out.reshape(values_src.shape[:-1] + (lam_master.size,))


def read_emission(path: Path):
    with h5py.File(path, 'r') as h:
        lam = lambda_from_h5(h)
        T = np.array(h['temperatures_K'], dtype=np.float64).reshape(-1)
        key = 'F_lambda_W_m3_sr_per_particle_shape' if 'F_lambda_W_m3_sr_per_particle_shape' in h else 'F_lambda'
        F = np.array(h[key], dtype=np.float32)
    order = np.argsort(lam)
    return lam[order], T, F[:, order]


def read_cia(path: Path):
    with h5py.File(path, 'r') as h:
        pair = str(h.attrs.get('species_pair', path.name.split('_')[0]))
        lam = lambda_from_h5(h)
        T = np.array(h['temperatures_K'], dtype=np.float64).reshape(-1)
        key = 'cia_m5_per_molecule2' if 'cia_m5_per_molecule2' in h else 'alpha'
        coeff = np.array(h[key], dtype=np.float32)
    order = np.argsort(lam)
    return pair, lam[order], T, coeff[:, order]


def main() -> None:
    parser = argparse.ArgumentParser(description='Export HDF5 precomputed tables to compact binary runtime assets.')
    parser.add_argument('--precomputed', type=Path, default=PROJECT_ROOT / 'data/precomputed')
    parser.add_argument('--output', type=Path, default=PROJECT_ROOT / 'data/runtime')
    parser.add_argument('--max-bins', type=int, default=8192, help='Use 75000 for full grid; 4096/8192 for fast tests.')
    args = parser.parse_args()

    op_paths = sorted((args.precomputed / 'opacity').glob('*_opacity_grid.h5'))
    if not op_paths:
        raise SystemExit(f'No opacity HDF5 tables found in {args.precomputed / "opacity"}')

    opacity_rows = [read_opacity(p) for p in op_paths]
    master_lam = opacity_rows[0][1]
    idx = downsample_indices(master_lam.size, args.max_bins)
    master_lam = master_lam[idx]
    args.output.mkdir(parents=True, exist_ok=True)

    with open(args.output / 'spectral_grid.bin', 'wb') as f:
        write_header(f)
        f.write(struct.pack('<Q', master_lam.size))
        f.write(np.ascontiguousarray(master_lam, dtype='<f8').tobytes())

    manifest = {'wavelength_bins': int(master_lam.size), 'opacities': [], 'cia': [], 'emission': None}

    with open(args.output / 'opacities.bin', 'wb') as f:
        write_header(f)
        f.write(struct.pack('<I', len(opacity_rows)))
        for species, lam, p, T, sigma in opacity_rows:
            sig2 = interpolate_to_master(lam, sigma, master_lam)
            write_string(f, species)
            f.write(struct.pack('<QQQ', p.size, T.size, master_lam.size))
            f.write(np.ascontiguousarray(p, dtype='<f8').tobytes())
            f.write(np.ascontiguousarray(T, dtype='<f8').tobytes())
            f.write(np.ascontiguousarray(sig2, dtype='<f4').tobytes())
            manifest['opacities'].append({'species': species, 'pressure': int(p.size), 'temperature': int(T.size)})
            print(f'Exported opacity {species}: p={p.size}, T={T.size}, lambda={master_lam.size}')

    cia_paths = sorted((args.precomputed / 'cia').glob('*_cia_grid.h5'))
    if cia_paths:
        with open(args.output / 'cia.bin', 'wb') as f:
            write_header(f)
            f.write(struct.pack('<I', len(cia_paths)))
            for pth in cia_paths:
                pair, lam, T, coeff = read_cia(pth)
                coeff2 = interpolate_to_master(lam, coeff, master_lam)
                write_string(f, pair)
                f.write(struct.pack('<QQ', T.size, master_lam.size))
                f.write(np.ascontiguousarray(T, dtype='<f8').tobytes())
                f.write(np.ascontiguousarray(coeff2, dtype='<f4').tobytes())
                manifest['cia'].append({'pair': pair, 'temperature': int(T.size)})
                print(f'Exported CIA {pair}: T={T.size}, lambda={master_lam.size}')

    em_path = args.precomputed / 'emission' / 'H3plus_emission_grid.h5'
    if em_path.exists():
        lam, T, F = read_emission(em_path)
        F2 = interpolate_to_master(lam, F, master_lam)
        with open(args.output / 'emission.bin', 'wb') as f:
            write_header(f)
            write_string(f, 'H3plus')
            f.write(struct.pack('<QQ', T.size, master_lam.size))
            f.write(np.ascontiguousarray(T, dtype='<f8').tobytes())
            f.write(np.ascontiguousarray(F2, dtype='<f4').tobytes())
        manifest['emission'] = {'species': 'H3plus', 'temperature': int(T.size)}
        print(f'Exported H3plus emission: T={T.size}, lambda={master_lam.size}')

    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(f'Runtime assets written to {args.output}')


if __name__ == '__main__':
    main()
