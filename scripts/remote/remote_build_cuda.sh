#!/usr/bin/env bash
set -euo pipefail

# Run this on the RTX 3060 remote machine after syncing the project there.
# It assumes CUDA toolkit and a C++ compiler are installed.

nvidia-smi || true
nvcc --version
cmake --preset cuda-rtx3060-release
cmake --build --preset cuda-rtx3060-release -j"$(nproc)"
./build-cuda-release/jupiter --cuda-smoke-test --use-cuda-rt --steps 1 --nr 32 --ntheta 24 --nphi 48 --assets data/runtime
