#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-data/output}"
FPS="${2:-12}"
MOVIE="${3:-${OUT_DIR}/jupiter_patch04_dashboard.mp4}"
FRAMES="${OUT_DIR}/frames/dashboard_%06d.ppm"

if ! command -v ffmpeg >/dev/null 2>&1; then
  echo "ffmpeg is not installed. Install it or convert ${OUT_DIR}/frames/dashboard_*.ppm manually." >&2
  exit 1
fi

if ! ls "${OUT_DIR}"/frames/dashboard_*.ppm >/dev/null 2>&1; then
  echo "No dashboard frames found at ${OUT_DIR}/frames/dashboard_*.ppm" >&2
  exit 1
fi

ffmpeg -y \
  -framerate "${FPS}" \
  -i "${FRAMES}" \
  -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
  -pix_fmt yuv420p \
  "${MOVIE}"

echo "Wrote ${MOVIE}"
