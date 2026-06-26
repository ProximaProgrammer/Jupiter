# JupiterRT Patch File List

Patch 04 touches the simulation core and visualizer. Generated outputs, build directories, raw data, precomputed opacity tables, and runtime binary assets are intentionally not part of the patch.

## Modified

- `config/default.sim`
- `src/config/SimConfig.hpp`
- `src/config/SimConfig.cpp`
- `src/pipeline/InitialConditions.hpp`
- `src/pipeline/InitialConditions.cpp`
- `src/physics/Dynamics.hpp`
- `src/physics/Dynamics.cpp`
- `src/io/Visualizer.hpp`
- `src/io/Visualizer.cpp`
- `src/io/Diagnostics.cpp`
- `src/main.cpp`

## Added

- `PATCH04_NOTES.md`
- `scripts/render/make_movie.sh`
