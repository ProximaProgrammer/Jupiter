# Deployment instructions

## Recommended: create a new project

From the parent directory that contains your old `Jupiter` folder:

```bash
mkdir JupiterRT
cd JupiterRT
unzip /path/to/JupiterRT_deployable.zip
```

If the zip extracts into a folder named `JupiterRT_deployable`, either work inside it or rename it:

```bash
mv JupiterRT_deployable/* JupiterRT_deployable/.[!.]* . 2>/dev/null || true
rmdir JupiterRT_deployable 2>/dev/null || true
```

## Link your giant data instead of copying

Assuming old project is `../Jupiter`:

```bash
mkdir -p data/raw
ln -s ../../Jupiter/data/raw/exomol data/raw/exomol
ln -s ../../Jupiter/data/raw/cia data/raw/cia
```

If symlinks are confusing, copy the folders instead, but they are large:

```bash
mkdir -p data/raw
cp -R ../Jupiter/data/raw/exomol data/raw/exomol
cp -R ../Jupiter/data/raw/cia data/raw/cia
```

## Build local CPU

```bash
cmake --preset mac-release
cmake --build --preset mac-release -j
./build-mac-release/jupiter --steps 1 --nr 16 --ntheta 12 --nphi 24 --assets data/runtime
```

## Preprocess assets

```bash
python3 -m pip install numpy h5py
python3 scripts/pipeline/build_all_assets.py --max-bins 8192
```

## Run with assets

```bash
./build-mac-release/jupiter --steps 3 --nr 32 --ntheta 24 --nphi 48 --assets data/runtime
```

## Initialize Git

```bash
git init
git add .
git commit -m "Initial JupiterRT deployable scaffold"
```

Do not commit huge data. `.gitignore` already ignores `data/raw`, `data/precomputed`, `data/runtime`, and `data/output` except placeholder files.
