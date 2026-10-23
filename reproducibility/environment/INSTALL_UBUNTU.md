# Installing SZcerebellumModel on Ubuntu

**Release preparation draft.** These commands document the 2026 working configuration and
must be tested on a completely fresh Linux installation before claiming broad portability.
The replay runner has been verified in `--prepare-only` mode using representative source YAMLs;
the full NEST recompilation and 96-replay process is a separate acceptance test.

## A. Environment for reconstruction and NEST simulation

### 1. System packages and Python

On Ubuntu, ensure Python **3.11** and its `venv` module are available. Ubuntu releases differ in whether `python3.11` is distributed in the default apt repositories. Install Python 3.11.17 through your distro's supported packages or another managed Python distribution, and confirm it before proceeding.

```bash
python3.11 --version
sudo apt-get update
sudo apt-get install -y git build-essential cmake pkg-config \
  libopenmpi-dev openmpi-bin libgsl-dev libboost-dev \
  libreadline-dev libltdl-dev libhdf5-dev
```

Use a **Linux-native filesystem** for the venv, NEST build, and compiled libraries. An external NTFS drive is appropriate for input and output datasets, not for building or activating the virtual environment.

### 2. Clone and create an isolated venv

```bash
git clone https://github.com/albertoarturovergani/SZcerebellumModel.git
cd SZcerebellumModel
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip wheel setuptools
python --version
which python
```

If the repository instead adopts a pinned upstream Git submodule, the clone command must use `--recurse-submodules`. This document assumes the preferred **one-copy source snapshot at repo root**.

### 3. MPI Python binding and BSB/cerebellum dependencies

Match `mpi4py` to the MPI used to compile NEST. The previously working 2026 build used `mpi4py==3.1.6` and required a controlled setuptools build environment. On systems where installing the wheel succeeds, a locally compiled wheel is not necessary.

```bash
printf 'setuptools==80.9.0\n' > /tmp/sz-mpi-build-constraints.txt
MPICC=mpicc python -m pip install \
  --build-constraint /tmp/sz-mpi-build-constraints.txt \
  --no-cache-dir --no-binary=mpi4py mpi4py==3.1.6
```

The published release should pin the remaining package dependencies in `reproducibility/environment/requirements-observed.txt`, generated from the working environment and manually checked for machine-local editable paths. After that file exists and has been verified:

```bash
python -m pip install -e .
python -m pip check
```

The editable install is required to register the BSB plugin `my_components = cerebellum:classmap` from `pyproject.toml`. Do not install a second, unrelated `cerebellar-models` plugin into the same venv.

### 4. Build NEST 3.7 with Python and MPI

NEST 3.7 is not provided by the `cerebellum` pip package. Build NEST against **the same active Python 3.11 environment**, enabling its MPI and Python bindings. The official NEST 3.7 CMake instructions document `-Dwith-mpi=ON`, `-Dwith-python=ON`, and `nest_vars.sh`.

```bash
mkdir -p "$HOME/src" "$HOME/build"
git clone --branch v3.7 --depth 1 \
  https://github.com/nest/nest-simulator.git "$HOME/src/nest-simulator-3.7"
cmake -S "$HOME/src/nest-simulator-3.7" \
  -B "$HOME/build/nest-3.7-sz" \
  -DCMAKE_INSTALL_PREFIX="$VIRTUAL_ENV" \
  -Dwith-python=ON -Dwith-mpi=ON -Dwith-gsl=ON \
  -DPYTHON_EXECUTABLE="$(command -v python)"
cmake --build "$HOME/build/nest-3.7-sz" --parallel 8
cmake --install "$HOME/build/nest-3.7-sz"
source "$VIRTUAL_ENV/bin/nest_vars.sh"
python -c 'import nest; print("NEST:", nest.version())'
```

If CMake selects a different Python interpreter or MPI library, stop and fix the build configuration. Do not copy a compiled `cerebmodule.so` from someone else's machine.

### 5. Compile the cerebellar NESTML module

The source lives in `cerebellum/nest_models/`, **not** in the generated network HDF5. NESTML's `generate_target()` generates, compiles, and installs a NEST extension module.

```bash
mkdir -p build/nestml_cerebmodule
python - <<'PY'
from pathlib import Path
import os
from pynestml.frontend.pynestml_frontend import generate_target
v=Path(os.environ['VIRTUAL_ENV'])
generate_target(
    input_path=str(Path('cerebellum/nest_models').resolve()),
    target_platform='NEST',
    target_path=str(Path('build/nestml_cerebmodule').resolve()),
    install_path=str(v/'lib/nest'),
    module_name='cerebmodule',
)
PY
export LD_LIBRARY_PATH="$VIRTUAL_ENV/lib/nest:$VIRTUAL_ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export NEST_MODULE_PATH="$VIRTUAL_ENV/lib/nest"
python - <<'PY'
import cerebellum
import nest
print('cerebellum:', cerebellum.__version__)
print('NEST:', nest.version())
nest.Install('cerebmodule')
assert 'eglif_cond_alpha_multisyn' in nest.node_models
print('cerebmodule: READY')
PY
```

The NESTML build can take several minutes and requires working compiler/CMake tools. In some toolchains `install_path` and the NEST module search environment need adjustment. Verify `nest.Install('cerebmodule')` before attempting network replay.

### 6. Smoke test and provenance

```bash
python - <<'PY'
import importlib.metadata as md
import cerebellum, nest
for name in ('bsb','bsb-core','bsb-nest','bsb-hdf5','mpi4py','nestml'):
    try:
        print(name, md.version(name))
    except md.PackageNotFoundError:
        print(name, 'not installed')
print('cerebellum', cerebellum.__version__)
print('NEST', nest.version())
PY
python -m pip check
mpirun --version | head -2
```

Next, **from a fresh clone and using only tracked inputs**, run the default **three-level smoke test** (`k=1.000`, `0.760`, `0.600`) for all six experimental cases (18 simulations). The package copies archived input paths into a writable runtime tree; output files must remain outside the Git checkout. The smoke test checks pipeline execution, not manuscript statistics:

```bash
python reproducibility/scripts/replay.py --case all --output "$HOME/SZ_replay_outputs" --mpi 8
python reproducibility/scripts/analyze_quick_results.py --input "$HOME/SZ_replay_outputs" --strict
```

Use the separate `--prepare-only` test with a **different** output directory to avoid leaving unfinished runtime trees in the simulation output location:

```bash
python reproducibility/scripts/replay.py --case all --prepare-only --output "$HOME/SZ_prepared_inputs" --mpi 8
```

To reproduce all 16 available levels per case (96 simulations), add `--full`:

```bash
python reproducibility/scripts/replay.py --case all --full --output "$HOME/SZ_replay_outputs" --mpi 8
```

Successful completed smoke runs will be skipped; the full run adds levels not yet completed. Do not attempt a full run in a directory created solely by `--prepare-only`: those runs contain no completed status or NIO.

The smoke plots and table can be inspected interactively using `notebooks/quick_smoke_test.ipynb`.
Install supplemental packages with `python -m pip install -r reproducibility/environment/requirements-notebooks.txt`, then run `python -m jupyter lab` from the repository root.
For the original manuscript notebook workflows, `replay.py` produces HDF5 and NIO only:
use `reproducibility/scripts/rebuild_intermediates.py --input "$HOME/SZ_replay_outputs" --execute`
to reconstruct historical working folders and attempt their `.pkl` intermediates
with the project's original post-processing scripts. See `notebooks/README.md`
for execution order (BASAL before EVOKED), sigma requirements, and limitations.

## B. Lightweight environment to analyze the original files from Zenodo

Researchers who wish to reproduce figures from the **archived HDF5 and NIO** can download these files from [Zenodo record 23391996](https://zenodo.org/records/23391996). The three historical Jupyter notebooks are included in `notebooks/historical/`, but may require the archived precomputed `.pkl` files and original directory structure. The independent quick-analysis notebook operates directly on newly simulated HDF5/NIO instead. A complete analysis-only environment should ultimately have its own audited `requirements-analysis.txt`, with `numpy`, `scipy`, `pandas`, `matplotlib`, `neo`, `nixio`, and whatever other packages the actual analysis scripts import.

```bash
python3.11 -m venv .venv-analysis
source .venv-analysis/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy scipy pandas matplotlib neo nixio h5py pyyaml
```

**Release prerequisite:** a complete locked analysis environment must be validated against the archived data before declaring every figure reproduction supported. Some existing analysis scripts may import BSB/NEST unnecessarily and require minor refactoring.


## Reproducibility caveats

- Historical BSB reconstruction had stochastic placement/connectivity without an archived global reconstruction seed: rebuilding scaffolds is not expected to reproduce identical neuron IDs or synapses. Later BSB releases provide better explicit random seed controls, but that does not retroactively reseed the 2025 runs.
- NEST simulation seeds **were already recorded in the 2025 YAML**, using 1234 for Basal/Auto and case-specific values for Kase. Those are distinct from the structural seed and depend on the underlying scaffold, execution configuration, and runtime.
- Preserve the exact YAML burst schedules, scientific parameters, and original morphology files. Do not regenerate the stimulus for historical input replay.
- In 2026, 96/96 experiments completed, but 2025-vs-2026 numerical identity was not claimed. An archived 2025 NIO was unavailable for one Basal-4Hz condition: validation scripts must mark missing comparisons as unavailable, not synthesize them.

## References

- [cerebellar-models v0.5.0 source](https://github.com/dbbs-lab/cerebellar-models/tree/v0.5.0)
- [NEST 3.7 CMake and MPI options](https://nest-simulator.readthedocs.io/en/v3.7/installation/cmake_options.html)
- [NEST 3.7 release](https://github.com/nest/nest-simulator/releases/tag/v3.7)
- [NESTML target generation and installation](https://nestml.readthedocs.io/en/latest/running/index.html)
- [Zenodo data archive](https://zenodo.org/records/23391996)
