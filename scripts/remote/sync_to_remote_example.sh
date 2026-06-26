#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   scripts/remote/sync_to_remote_example.sh user@host /remote/path/JupiterRT
# This intentionally excludes giant raw/precomputed/runtime data. Sync data separately or build assets on remote.

REMOTE="${1:?remote like user@host}"
DEST="${2:?remote destination path}"

rsync -avz --delete \
  --exclude 'build*/' \
  --exclude 'data/raw/' \
  --exclude 'data/precomputed/' \
  --exclude 'data/runtime/' \
  --exclude 'data/output/' \
  --exclude '.git/' \
  ./ "$REMOTE:$DEST/"

echo "Synced source to $REMOTE:$DEST"
echo "Next on remote: cd $DEST && scripts/remote/remote_build_cuda.sh"
