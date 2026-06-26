# JupiterRT

Clean deployable Jupiter forward-model scaffold for the current project.

This is a fresh project because the old repository still builds the toy `playground` target. The intended workflow is:

1. create/unzip this as a new project folder;
2. symlink or copy your giant data folders into `data/raw/`;
3. build CPU on your Mac for correctness;
4. sync the project to the RTX 3060 machine and build with CUDA.

## Included physics/scaffold

- Spherical cell geometry with exact shell-sector volumes and face areas.
- Conserved state and primitive state separation.
- Layer-aware EOS closures:
  - upper atmosphere: H2/He ideal-gas-like closure;
  - molecular hydrogen layer;
  - supercritical hydrogen Mie-Gruneisen-like cold pressure;
  - metallic hydrogen cold pressure plus electron-degeneracy floor.
- Source updates:
  - gravity;
  - pressure gradients;
  - rotating-frame centrifugal and Coriolis sources.
- Conservative epsilon cell-to-cell transport through face neighbors only; corner terms are deliberately ignored as O(epsilon^2).
- Radiative transfer per cell:

  `I_out(lambda) = I_in(lambda) exp(-tau_lambda) + S_lambda (1 - exp(-tau_lambda))`

- Doppler transform into the cell comoving frame using the local matter velocity plus

  `v_rot = Omega_J x r`

  because the grid is rotating Eulerian.
- Runtime spectral assets from ExoMol/TauREx HDF5 and CIA text/HDF5 tables.
- CUDA/RTX 3060 build path and a CUDA RT layer-update kernel.

## What this is not yet

This is not a final publication-grade Jupiter interior model. The EOS is structurally correct and replaceable, but analytic. Later, the clean replacement is a tabulated H/He EOS behind the same `EOS::ConservedToPrimitive` interface.

## CPU build on Mac

```bash
cmake --preset mac-release
cmake --build --preset mac-release -j
./build-mac-release/jupiter --steps 1 --nr 32 --ntheta 24 --nphi 48 --assets data/runtime
```

If no runtime assets exist yet, the code uses analytic fallback opacities/emission.

## Runtime asset preprocessing

After linking/copying data into `data/raw`, install Python dependencies:

```bash
python3 -m pip install numpy h5py
```

Fast first pass:

```bash
python3 scripts/pipeline/build_all_assets.py --max-bins 8192
```

Full spectral pass, if memory/time is acceptable:

```bash
python3 scripts/pipeline/build_all_assets.py --max-bins 75000
```

Expected outputs:

```text
data/runtime/spectral_grid.bin
data/runtime/opacities.bin
data/runtime/cia.bin
data/runtime/emission.bin
data/runtime/manifest.json
```

Then run:

```bash
./build-mac-release/jupiter --steps 3 --nr 48 --ntheta 32 --nphi 64 --assets data/runtime
```

Output:

```text
data/output/synthetic_radial_ray.csv
```

## Remote RTX 3060 build

On the remote host, after syncing:

```bash
cmake --preset cuda-rtx3060-release
cmake --build --preset cuda-rtx3060-release -j
./build-cuda-release/jupiter --cuda-smoke-test --use-cuda-rt --steps 3 --nr 48 --ntheta 32 --nphi 64 --assets data/runtime
```

The CUDA kernel currently accelerates the per-layer spectral recurrence once opacities/source vectors have been computed. The next performance step is moving opacity interpolation and many rays onto the GPU, but this version is already build- and SSH-ready for an RTX 3060.
