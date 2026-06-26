from __future__ import annotations
import argparse
from pathlib import Path
from datetime import datetime, timezone

import h5py
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOSCHMIDT_M3 = 2.686780111e25  # molecules / m^3 at 1 amagat


def is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def parse_hitran_cia_blocks(inp: Path, species_pair: str):
    """
    HITRAN CIA files are usually block-structured:

        H2-H2  20.000 10000.000 9981  200.0 ...
        20.0000  2.668E-47
        21.0000  2.931E-47
        ...

    Header contains temperature. Data rows contain only wavenumber and CIA coefficient.
    """
    rows = []
    current_T = None

    with inp.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith("!"):
                continue

            parts = stripped.replace(",", " ").split()

            # Header line: contains species pair and enough metadata.
            # Example:
            # H2-H2 20.000 10000.000 9981 200.0 8.788E-45 -.999 6
            if parts[0].upper() == species_pair.upper():
                if len(parts) < 5:
                    raise RuntimeError(f"Malformed CIA header line: {stripped}")
                current_T = float(parts[4])
                continue

            # Some files may have headers where the pair token is not first because of spacing/labels.
            if species_pair.upper() in [p.upper() for p in parts]:
                numeric = [float(p) for p in parts if is_float(p)]
                if len(numeric) >= 4:
                    current_T = numeric[3]
                    continue

            # Data row: wavenumber alpha.
            if len(parts) >= 2 and is_float(parts[0]) and is_float(parts[1]):
                if current_T is None:
                    # Numeric line before any header: ignore instead of inventing T.
                    continue
                wn = float(parts[0])
                alpha = float(parts[1])
                rows.append((wn, current_T, alpha))

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert HITRAN CIA text to JupiterRT HDF5.")
    parser.add_argument("--species-pair", required=True, help="e.g. H2-H2 or H2-He")
    parser.add_argument("--input", type=Path, default=None, help="Default: data/raw/cia/<pair>/opacity_data.cia")
    parser.add_argument("--output", type=Path, default=None, help="Default: data/precomputed/cia/<pair>_cia_grid.h5")
    parser.add_argument("--raw-units", choices=["cm-1-amagat-2", "m5-molecule-2"], default="cm-1-amagat-2")
    args = parser.parse_args()

    inp = args.input or PROJECT_ROOT / "data" / "raw" / "cia" / args.species_pair / "opacity_data.cia"
    out = args.output or PROJECT_ROOT / "data" / "precomputed" / "cia" / f"{args.species_pair}_cia_grid.h5"

    inp = inp if inp.is_absolute() else PROJECT_ROOT / inp
    out = out if out.is_absolute() else PROJECT_ROOT / out

    if not inp.exists():
        raise FileNotFoundError(inp)

    out.parent.mkdir(parents=True, exist_ok=True)

    rows = parse_hitran_cia_blocks(inp, args.species_pair)

    if not rows:
        raise RuntimeError(f"No numeric CIA rows parsed from {inp}")

    data = np.array(rows, dtype=np.float64)

    wn_unique = np.unique(data[:, 0])
    t_unique = np.unique(data[:, 1])

    coeff = np.zeros((t_unique.size, wn_unique.size), dtype=np.float64)
    counts = np.zeros_like(coeff)

    wn_index = {v: i for i, v in enumerate(wn_unique)}
    t_index = {v: i for i, v in enumerate(t_unique)}

    for wn, T, alpha in data:
        i = t_index[T]
        j = wn_index[wn]
        coeff[i, j] += alpha
        counts[i, j] += 1.0

    mask = counts > 0
    coeff[mask] /= counts[mask]

    if args.raw_units == "cm-1-amagat-2":
        # HITRAN CIA coefficient:
        # cm^-1 amagat^-2
        #
        # Convert cm^-1 -> m^-1 by multiplying by 100.
        # Convert amagat^-2 -> (molecule m^-3)^-2 using Loschmidt number.
        coeff = coeff * 100.0 / (LOSCHMIDT_M3 ** 2)

    with h5py.File(out, "w") as h:
        h.attrs["species_pair"] = args.species_pair
        h.attrs["source"] = str(inp)
        h.attrs["date_created"] = datetime.now(timezone.utc).isoformat()
        h.attrs["coefficient_units"] = "m^5 molecule^-2"
        h.create_dataset("wavenumber_cm1", data=wn_unique, compression="gzip", compression_opts=4)
        h.create_dataset("temperatures_K", data=t_unique, compression="gzip", compression_opts=4)
        h.create_dataset("cia_m5_per_molecule2", data=coeff.astype(np.float32), compression="gzip", compression_opts=4)

    print(f"Created {out}: T={t_unique.size}, wn={wn_unique.size}, rows={len(rows)}")


if __name__ == "__main__":
    main()