from pathlib import Path
import h5py

root = Path('data/precomputed')
for f in sorted(root.rglob('*.h5')):
    print('=' * 80)
    print(f)
    with h5py.File(f, 'r') as h:
        for k in h.keys():
            print(k, h[k].shape, h[k].dtype)
