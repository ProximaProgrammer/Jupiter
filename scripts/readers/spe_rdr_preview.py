#!/usr/bin/env python3
"""Preview a JIRAM-like SPE/RDR binary using its detached or inline PDS label.

This does not pretend to solve every PDS edge case; it is a guardrail around the
failure mode you described where the calibration/read returned empty-looking data.
It prints the interpreted shape, dtype, finite range, and the first valid spectrum.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import numpy as np


_NUM = r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?"


def read_label(path: Path) -> str:
    if path.suffix.lower() in {".lbl", ".xml"}:
        return path.read_text(errors="replace")
    for cand in [path.with_suffix(".LBL"), path.with_suffix(".lbl"), Path(str(path) + ".LBL"), Path(str(path) + ".lbl")]:
        if cand.exists():
            return cand.read_text(errors="replace")
    # Some SPE files carry an ASCII header. Read only the first 64 KiB.
    return path.read_bytes()[:65536].decode("latin-1", errors="ignore")


def scalar(label: str, names: list[str], default: Any = None) -> Any:
    for name in names:
        m = re.search(rf"\b{name}\b\s*=\s*(?:\"([^\"]+)\"|({ _NUM })|([^\s,]+))", label, flags=re.I | re.X)
        if m:
            v = next(g for g in m.groups() if g is not None)
            try:
                return int(float(v))
            except ValueError:
                return v.strip()
    return default


def dtype_from_label(label: str) -> np.dtype:
    sample_type = str(scalar(label, ["SAMPLE_TYPE", "DATA_TYPE"], "MSB_REAL")).upper()
    bits = int(scalar(label, ["SAMPLE_BITS", "BYTES_PER_PIXEL"], 32))
    if "REAL" in sample_type or "FLOAT" in sample_type:
        if bits in (4, 32):
            return np.dtype(">f4" if "MSB" in sample_type else "<f4")
        if bits in (8, 64):
            return np.dtype(">f8" if "MSB" in sample_type else "<f8")
    signed = "UNSIGNED" not in sample_type
    endian = ">" if "MSB" in sample_type else "<"
    if bits in (8,): return np.dtype("u1" if not signed else "i1")
    if bits in (16,): return np.dtype(endian + ("i2" if signed else "u2"))
    if bits in (32,): return np.dtype(endian + ("i4" if signed else "u4"))
    raise ValueError(f"Unsupported SAMPLE_TYPE/SAMPLE_BITS: {sample_type} / {bits}")


def main() -> None:
    p = argparse.ArgumentParser(description="Preview SPE/RDR data safely.")
    p.add_argument("spe", type=Path)
    p.add_argument("--offset-bytes", type=int, default=None)
    p.add_argument("--lines", type=int, default=None)
    p.add_argument("--samples", type=int, default=None)
    p.add_argument("--bands", type=int, default=None)
    args = p.parse_args()

    label = read_label(args.spe)
    dtype = dtype_from_label(label)
    lines = args.lines or int(scalar(label, ["LINES", "CORE_ITEMS"], 1))
    samples = args.samples or int(scalar(label, ["LINE_SAMPLES", "SAMPLES", "BAND_SAMPLES"], 336))
    bands = args.bands or int(scalar(label, ["BANDS"], 1))
    offset = args.offset_bytes
    if offset is None:
        ptr = scalar(label, ["^SPECTRAL_QUBE", "^QUBE", "^IMAGE", "^TABLE"], 1)
        rec = int(scalar(label, ["RECORD_BYTES"], 0) or 0)
        offset = max(0, (int(ptr) - 1) * rec) if isinstance(ptr, int) and rec > 0 else 0

    raw = np.fromfile(args.spe, dtype=dtype, offset=offset)
    expected = lines * samples * bands
    print(f"file={args.spe}")
    print(f"dtype={dtype}, offset={offset}, raw_values={raw.size}, expected={expected} ({lines}x{samples}x{bands})")
    if raw.size == 0:
        raise SystemExit("No values were read. Check offset, label pointer, and endian/sample type.")
    if raw.size >= expected and expected > 0:
        arr = raw[:expected].reshape((lines, bands, samples))
    else:
        arr = raw.reshape((1, 1, raw.size))
    finite = arr[np.isfinite(arr)]
    print(f"shape={arr.shape}, finite={finite.size}/{arr.size}")
    if finite.size:
        print(f"min={finite.min():.8g}, max={finite.max():.8g}, mean={finite.mean():.8g}")
    first = arr.reshape((-1, arr.shape[-1]))[0]
    print("first_spectrum_head=", np.array2string(first[: min(20, first.size)], precision=6, separator=", "))


if __name__ == "__main__":
    main()
