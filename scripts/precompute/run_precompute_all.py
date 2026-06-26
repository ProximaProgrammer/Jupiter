from __future__ import annotations
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(args):
    print('\n$', ' '.join(str(a) for a in args))
    subprocess.check_call([str(a) for a in args], cwd=ROOT)


def maybe_run(args, required_path: Path):
    if required_path.exists():
        run(args)
    else:
        print(f'Skip missing input: {required_path}')


def main() -> None:
    ex = ROOT / 'data/raw/exomol'
    species = ['CH4', 'NH3', 'PH3', 'H2O']
    for sp in species:
        inp = ex / sp / 'opacity_data.h5'
        maybe_run([sys.executable, 'scripts/precompute/precompute_opacity.py', sp, inp], inp)

    h3_states = ex / 'H3plus' / '1H3_p__MiZATeP.states'
    maybe_run([sys.executable, 'scripts/precompute/precompute_h3plus_approx.py', h3_states], h3_states)

    for pair in ['H2-H2', 'H2-He']:
        inp = ROOT / 'data/raw/cia' / pair / 'opacity_data.cia'
        maybe_run([sys.executable, 'scripts/precompute/precompute_cia_hitran.py', '--species-pair', pair], inp)


if __name__ == '__main__':
    main()
