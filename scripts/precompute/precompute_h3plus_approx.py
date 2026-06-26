from __future__ import annotations
import argparse
from pathlib import Path
import h5py
import numpy as np

from precompute_common import ensure_parent, make_master_grid, write_common_axes

HC_OVER_K = 1.438776877e-2
C1 = 1.191042972e-16


def planck_lambda(wavelength_m: np.ndarray, T: float) -> np.ndarray:
    x = HC_OVER_K / (wavelength_m * T)
    return C1 / (wavelength_m**5 * (np.exp(np.clip(x, 1e-12, 700.0)) - 1.0))


def gaussian_log_lambda(wavelength_m: np.ndarray, center_um: float, width_dex: float) -> np.ndarray:
    loglam = np.log10(wavelength_m / 1e-6)
    return np.exp(-0.5 * ((loglam - np.log10(center_um)) / width_dex) ** 2)


def h3plus_template(wavelength_m: np.ndarray) -> np.ndarray:
    template = (
        1.0 * gaussian_log_lambda(wavelength_m, 3.5, 0.055)
        + 0.55 * gaussian_log_lambda(wavelength_m, 3.9, 0.045)
        + 0.35 * gaussian_log_lambda(wavelength_m, 2.1, 0.060)
    )
    area = np.trapezoid(template, wavelength_m)
    return template / area if area > 0 else template


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a compact approximate H3+ emission template. The optional states file is recorded only for provenance.")
    parser.add_argument("states_file", nargs="?", default=None)
    parser.add_argument("--output", type=Path, default=Path("data/precomputed/emission/H3plus_emission_grid.h5"))
    parser.add_argument("--n-bins", type=int, default=75_000)
    parser.add_argument("--lambda-min-m", type=float, default=1e-11)
    parser.add_argument("--lambda-max-m", type=float, default=1e4)
    parser.add_argument("--temperatures", type=float, nargs="*", default=[300, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 4000, 6000])
    parser.add_argument("--scale", type=float, default=1.0)
    args = parser.parse_args()

    grid = make_master_grid(args.n_bins, args.lambda_min_m, args.lambda_max_m)
    T = np.array(args.temperatures, dtype=np.float64)
    template = h3plus_template(grid.wavelength_m)
    ensure_parent(args.output)
    with h5py.File(args.output, "w") as out:
        write_common_axes(out, grid)
        out.create_dataset("temperatures_K", data=T)
        ds = out.create_dataset("F_lambda_W_m3_sr_per_particle_shape", shape=(T.size, args.n_bins), dtype="f4", compression="gzip", compression_opts=4, chunks=(1, min(args.n_bins, 4096)))
        for i, temp in enumerate(T):
            ds[i, :] = (args.scale * template * planck_lambda(grid.wavelength_m, temp)).astype(np.float32)
        out.attrs["species"] = "H3plus"
        out.attrs["model"] = "smooth H3+ near-IR thermal band template; replace with line-summed compact table when ready"
        if args.states_file:
            out.attrs["states_file"] = str(args.states_file)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
