from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(cmd):
    print('\n$', ' '.join(str(c) for c in cmd))
    subprocess.check_call([str(c) for c in cmd], cwd=ROOT)


def main() -> None:
    parser = argparse.ArgumentParser(description='Run precompute conversion and export C++ runtime assets.')
    parser.add_argument('--max-bins', type=int, default=8192)
    args = parser.parse_args()
    run([sys.executable, 'scripts/precompute/run_precompute_all.py'])
    run([sys.executable, 'scripts/precompute/export_runtime_assets.py', '--max-bins', str(args.max_bins)])


if __name__ == '__main__':
    main()
