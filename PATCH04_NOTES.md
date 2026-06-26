# JupiterRT Patch 04 Notes

Patch 04 converts the project from a static smoke-test style renderer into a time-evolving simulation loop with a larger dashboard visualizer.

## What changed

- Added conservative finite-volume face advection as the default transport step.
- Kept the older epsilon-shift transport as a fallback switch.
- Expanded source terms: gravity with a smooth enclosed-mass model, pressure gradients, rotating-frame centrifugal and Coriolis terms, convective buoyancy, velocity drag, effective viscosity, thermal diffusion, grey radiative relaxation, and a small H3+ chemistry hook.
- Added deterministic thermal and velocity perturbations to the initial conditions so the simulation visibly evolves instead of remaining a perfectly radial static gradient.
- Added boundary projection at the inner/outer radii and polar caps.
- Expanded diagnostics with energy, momentum, Mach number, and H3+ fraction.
- Reworked the visualizer to write a large multi-panel dashboard frame: `data/output/frames/dashboard_XXXXXX.ppm` and `dashboard_latest.ppm`.
- Added `scripts/render/make_movie.sh` to stitch dashboard frames into an MP4 when `ffmpeg` is installed.

## Quick run

```bash
cd ~/Documents/Jupiter_v2
cmake -S . -B build-mac-release
cmake --build build-mac-release
./build-mac-release/jupiter --nr 32 --ntheta 24 --nphi 48 --steps 120 --dt 0.1
open data/output/frames/dashboard_latest.ppm
```

To make an MP4:

```bash
scripts/render/make_movie.sh data/output 12
open data/output/jupiter_patch04_dashboard.mp4
```

## Important limitation

This patch adds the infrastructure and first-pass closures for the physics we discussed. It is not a validated Jupiter interior model yet. The EOS, turbulent transport, convective source, chemistry, and grey RT-relaxation terms are replaceable closures; the formal spectral ray calculation remains in `RadiativeTransfer.cpp` and still uses `I_out = I_in exp(-tau) + S(1-exp(-tau))`.
