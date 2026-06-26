from __future__ import annotations
import argparse
from pathlib import Path
from precompute_common import print_h5_tree

parser = argparse.ArgumentParser(description='Print datasets inside an HDF5 file.')
parser.add_argument('file', type=Path)
args = parser.parse_args()
print_h5_tree(args.file)
