from __future__ import annotations

from pathlib import Path
from typing import Iterable
import h5py


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def find_dataset_by_candidates(h: h5py.File | h5py.Group, candidates: Iterable[str]) -> str | None:
    wanted = {c.lower() for c in candidates}
    found: list[str] = []

    def visit(name: str, obj):
        if isinstance(obj, h5py.Dataset):
            base = name.split('/')[-1].lower()
            if base in wanted or name.lower() in wanted:
                found.append(name)

    h.visititems(visit)
    return found[0] if found else None


def print_h5_tree(path: Path) -> None:
    with h5py.File(path, 'r') as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f'{name}: shape={obj.shape}, dtype={obj.dtype}, attrs={dict(obj.attrs)}')
            else:
                print(f'{name}/')
        f.visititems(visit)

from dataclasses import dataclass
import numpy as np

@dataclass
class MasterGrid:
    wavelength_m: np.ndarray
    wavenumber_cm1: np.ndarray


def make_master_grid(n_bins: int, lambda_min_m: float, lambda_max_m: float) -> MasterGrid:
    wavelength_m = np.exp(np.linspace(np.log(lambda_min_m), np.log(lambda_max_m), n_bins)).astype(np.float64)
    wavenumber_cm1 = 1.0 / (100.0 * wavelength_m)
    return MasterGrid(wavelength_m=wavelength_m, wavenumber_cm1=wavenumber_cm1)


def write_common_axes(h, grid: MasterGrid) -> None:
    h.create_dataset('wavelength_m', data=grid.wavelength_m, compression='gzip', compression_opts=4)
    h.create_dataset('wavenumber_cm1', data=grid.wavenumber_cm1, compression='gzip', compression_opts=4)
