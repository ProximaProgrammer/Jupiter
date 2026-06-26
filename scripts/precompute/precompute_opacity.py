from __future__ import annotations
import argparse
from pathlib import Path
from datetime import datetime, timezone

import h5py
import numpy as np

from precompute_common import ensure_parent, find_dataset_by_candidates

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def read_units(dset, default="unknown"):
    return dset.attrs.get("units", default)


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert an ExoMol/TauREx HDF5 opacity file to the project precomputed format.")
    parser.add_argument("species", help="CH4, NH3, PH3, H2O, ...")
    parser.add_argument("input_file", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    inp = args.input_file if args.input_file.is_absolute() else PROJECT_ROOT / args.input_file
    out = args.output or PROJECT_ROOT / "data" / "precomputed" / "opacity" / f"{args.species}_opacity_grid.h5"
    ensure_parent(out)
    if not inp.exists():
        raise FileNotFoundError(inp)

    with h5py.File(inp, "r") as f:
        wn_key = find_dataset_by_candidates(f, ["bin_edges", "wavenumber", "wavenumber_cm1", "wn", "nu"])
        p_key = find_dataset_by_candidates(f, ["p", "pressure", "pressures", "pressure_bar"])
        t_key = find_dataset_by_candidates(f, ["t", "temperature", "temperatures", "temperature_K"])
        sig_key = find_dataset_by_candidates(f, ["xsecarr", "sigma", "xsec", "cross_section", "cross-sections"])
        if not all([wn_key, p_key, t_key, sig_key]):
            raise KeyError(f"Could not identify required datasets in {inp}. Found: {[r['name'] for r in []]}")
        wn = np.array(f[wn_key], dtype=np.float64).reshape(-1)
        p = np.array(f[p_key], dtype=np.float64).reshape(-1)
        t = np.array(f[t_key], dtype=np.float64).reshape(-1)
        sigma = np.array(f[sig_key], dtype=np.float64)
        sigma_units = str(read_units(f[sig_key]))

    # Normalize shape to (pressure, temperature, wavelength/wavenumber).
    if sigma.shape != (len(p), len(t), len(wn)):
        if sigma.shape == (len(t), len(p), len(wn)):
            sigma = np.transpose(sigma, (1, 0, 2))
        else:
            raise ValueError(f"Unrecognized sigma shape {sigma.shape}; expected {(len(p), len(t), len(wn))} or {(len(t), len(p), len(wn))}")

    # Common TauREx/ExoMol xsec is cm^2 molecule^-1. If metadata says m2, do not convert.
    if "m^2" in sigma_units or "m2" in sigma_units:
        sigma_m2 = sigma
    else:
        sigma_m2 = sigma * 1.0e-4

    with h5py.File(out, "w") as g:
        g.attrs["species"] = args.species
        g.attrs["source"] = str(inp)
        g.attrs["date_created"] = datetime.now(timezone.utc).isoformat()
        g.attrs["spectral_coordinate"] = "wavenumber"
        g.attrs["wavenumber_units"] = "cm^-1"
        g.attrs["pressure_units"] = "bar"
        g.attrs["temperature_units"] = "K"
        g.attrs["sigma_units"] = "m^2/molecule"
        g.create_dataset("wavenumber", data=wn, compression="gzip", compression_opts=4)
        g.create_dataset("pressure", data=p, compression="gzip", compression_opts=4)
        g.create_dataset("temperature", data=t, compression="gzip", compression_opts=4)
        g.create_dataset("sigma", data=sigma_m2.astype(np.float32), compression="gzip", compression_opts=4)
    print(f"Created {out}")
    print(f"  p={p.shape}, T={t.shape}, wn={wn.shape}, sigma={sigma_m2.shape}")


if __name__ == "__main__":
    main()
