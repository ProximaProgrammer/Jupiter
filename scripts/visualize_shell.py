#!/usr/bin/env python3
"""Create an interactive HTML view of the Patch 04 second-to-outermost shell.

Usage:
  python3 scripts/visualize_shell.py data/output/shells/outer_interior_shell_latest.csv

Requires plotly only for the HTML view:
  python3 -m pip install plotly
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


def read_rows(path: Path) -> dict[str, list[float]]:
    rows: dict[str, list[float]] = {}
    with path.open(newline="") as f:
        # skip leading comments while preserving the header line
        while True:
            pos = f.tell()
            line = f.readline()
            if not line:
                raise SystemExit(f"No CSV header found in {path}")
            if not line.startswith("#"):
                f.seek(pos)
                break
        reader = csv.DictReader(f)
        for name in reader.fieldnames or []:
            rows[name] = []
        for r in reader:
            for k, v in r.items():
                try:
                    rows[k].append(float(v))
                except (TypeError, ValueError):
                    rows[k].append(math.nan)
    return rows


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/output/shells/outer_interior_shell_latest.csv")
    if not path.exists():
        raise SystemExit(f"Missing shell CSV: {path}")
    rows = read_rows(path)
    try:
        import plotly.graph_objects as go
    except Exception as exc:  # pragma: no cover
        raise SystemExit("Plotly is needed for the interactive viewer: python3 -m pip install plotly") from exc

    x, y, z = rows["x_m"], rows["y_m"], rows["z_m"]
    T = rows["T_K"]
    vx, vy, vz = rows["vx_inertial_m_s"], rows["vy_inertial_m_s"], rows["vz_inertial_m_s"]

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z, mode="markers",
        marker={"size": 3, "color": T, "colorscale": "Turbo", "colorbar": {"title": "T [K]"}},
        text=[f"T={t:.2f} K<br>v=({a:.1f},{b:.1f},{c:.1f}) m/s" for t, a, b, c in zip(T, vx, vy, vz)],
        name="temperature shell",
    ))

    # Sparse velocity cone layer to keep the HTML light.
    stride = max(1, len(x) // 900)
    fig.add_trace(go.Cone(
        x=x[::stride], y=y[::stride], z=z[::stride],
        u=vx[::stride], v=vy[::stride], w=vz[::stride],
        sizemode="absolute", sizeref=1.8e6, anchor="tail", showscale=False,
        name="inertial velocity",
    ))
    fig.update_layout(
        title="JupiterRT Patch 04: second-to-outermost shell (outer layer is ghost)",
        scene={"aspectmode": "data", "xaxis_title": "x [m]", "yaxis_title": "y [m]", "zaxis_title": "z [m]"},
        margin={"l": 0, "r": 0, "t": 40, "b": 0},
    )
    out = path.with_suffix(".html")
    fig.write_html(out)
    latest = path.parent / "outer_interior_shell_latest.html"
    if out != latest:
        latest.write_text(out.read_text())
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
