#!/usr/bin/env python3
"""Predict Galilean satellite vectors for later JIRAM/Juno ray geometry.

This is intentionally optional: it requires spiceypy and a user-supplied metakernel.
It uses the standard SPICE flow: furnsh(...) loads kernels, utc2et converts time,
spkpos returns target position relative to observer, and pxform can rotate frames.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict Galilean satellite positions from a SPICE metakernel.")
    parser.add_argument("--metakernel", required=True, help="Path to a SPICE metakernel (.tm) that loads LSK/SPK/PCK/FK/CK as available.")
    parser.add_argument("--utc", required=True, help="UTC time, e.g. 2016-08-27T12:00:00")
    parser.add_argument("--observer", default="JUNO", help="Observer body, default JUNO")
    parser.add_argument("--frame", default="J2000", help="Output frame, default J2000")
    parser.add_argument("--abcorr", default="LT+S", help="Aberration correction, default LT+S")
    parser.add_argument("--out", default="data/output/spice_satellite_vectors.csv")
    args = parser.parse_args()

    try:
        import spiceypy as spice
    except Exception as exc:  # pragma: no cover
        raise SystemExit("Install spiceypy first: python3 -m pip install spiceypy") from exc

    spice.furnsh(str(Path(args.metakernel)))
    try:
        et = spice.utc2et(args.utc)
        targets = ["IO", "EUROPA", "GANYMEDE", "CALLISTO"]
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["utc", "et", "target", "observer", "frame", "abcorr", "x_km", "y_km", "z_km", "light_time_s", "range_km"])
            for target in targets:
                pos, lt = spice.spkpos(target, et, args.frame, args.abcorr, args.observer)
                rng = float((pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2) ** 0.5)
                w.writerow([args.utc, et, target, args.observer, args.frame, args.abcorr, pos[0], pos[1], pos[2], lt, rng])
        print(f"Wrote {args.out}")
    finally:
        spice.kclear()


if __name__ == "__main__":
    main()
