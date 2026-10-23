# SZcerebellumModel

A computational framework for investigating **schizophrenia-like structural and synaptic alterations** in a biologically grounded, spiking **cerebellar microcircuit model**.

This repository extends the [Brain Scaffold Builder (BSB)](https://github.com/dbbs-lab/bsb) and [cerebellum](https://github.com/dbbs-lab/cerebellar-models/tree/v0.4.0) computational frameworks to investigate how progressive neuronal degeneration and resilience mechanisms influence cerebellar network dynamics.

## Key Features

- **Progressive neuronal degeneration:** Morphological shrinking, dendritic pruning, and thinning.
- **Synaptic alterations:** Configurable modifications of synaptic connectivity and transmission.
- **Resilience mechanisms:** Computational implementations of mechanisms that may mitigate network dysfunction.
- **Automated workflows:** End-to-end pipelines for network construction, simulation, and analysis.
- **Multiple degeneration levels:** Systematic simulations across different degrees of schizophrenia-like alterations.
- **Multiscale analyses:** Firing rates, cross-correlation, transfer entropy (TE), phase-locking value (PLV), local field potentials (LFP), and structural network metrics.

The framework supports cerebellar neuronal populations and pathways involving **mossy fibers, granule cells, Golgi cells, basket cells, stellate cells, Purkinje cells, deep cerebellar nuclei (DCN), and inferior olive (IO)**, depending on the simulation configuration.

---

## Repository Structure

The repository includes the integrated 2026 source snapshot and historical
simulation and analysis workflows:

```text
SZcerebellumModel/
├── cerebellum/         # Cerebellar model source snapshot
├── configurations/     # Archived simulation configurations
├── morphologies/       # Neuronal morphologies
├── notebooks/          # Historical analysis notebooks
├── reproducibility/    # Environment, replay tools and manifests
├── validation/         # Smoke tests and numerical baselines
├── scripts_py/         # Python processing and analysis tools
├── scripts_sh/         # Simulation workflow scripts
├── pyproject.toml
└── README.md
```

See [reproducibility](reproducibility/README.md),
[Ubuntu installation](reproducibility/environment/INSTALL_UBUNTU.md), and
[validation](validation/README.md) for the current workflow and its limits.

### `scripts_sh/` — Shell Automation Scripts

Contains Bash scripts for orchestrating the computational workflow, including:

- Generation of morphologically altered neuronal networks
- Network construction using BSB
- Simulation configuration and execution
- Batch processing across degeneration levels
- Post-simulation analysis
- Automated end-to-end workflows

### `scripts_py/` — Python Processing Modules

Contains Python scripts and utilities for:

- Neuronal morphology transformations
- Synaptic and connectivity scaling
- YAML configuration generation
- Network construction and simulation setup
- Electrophysiological and structural analyses
- Extraction and aggregation of simulation results
- Visualization and comparison across degeneration levels

---

## 1. Prerequisites

The workflow builds upon the BSB and `cerebellum` frameworks and requires a compatible Python environment, MPI, and the NEST simulation backend.

### 1.1 Brain Scaffold Builder (BSB)

Official repository: [dbbs-lab/bsb](https://github.com/dbbs-lab/bsb)

The simulations described in the manuscript used **BSB 4.4** ([release `v4.4.0`](https://github.com/dbbs-lab/bsb/tree/v4.4.0)). In a dedicated Python environment, install the version-matched BSB core:

```bash
python -m pip install "bsb-core[parallel]==4.4.0"
```

Additional dependencies may be required depending on the simulation backend and analysis scripts.

### 1.2 Cerebellum Framework

Official repository: [dbbs-lab/cerebellar-models](https://github.com/dbbs-lab/cerebellar-models)

Clone the **cerebellum 0.4.0** release, which contains the upstream `cerebellum` package and directory structure used in the paper:

```bash
git clone --branch v0.4.0 https://github.com/dbbs-lab/cerebellar-models.git cerebellum
cd cerebellum
python -m pip install -e .
cd ..
```

For reproducing the analyses described in the associated manuscript, the reported software versions are **BSB 4.4** and **cerebellum 0.4.0**. Compatibility with newer releases should be verified.

---

## 2. Installation and Integration

**2026 integrated snapshot:** The instructions below describe the legacy integration of scripts into a separate upstream `cerebellum` checkout. For the integrated 2026 source snapshot, consult the [Ubuntu installation guide](reproducibility/environment/INSTALL_UBUNTU.md). Independent installation on a completely fresh Ubuntu system remains to be verified.

The scripts provided by `SZcerebellumModel` are designed to operate within the root directory of the upstream `cerebellum` repository.

### Step 1 — Clone SZcerebellumModel

From the parent directory containing your `cerebellum` installation, run:

```bash
git clone https://github.com/albertoarturovergani/SZcerebellumModel.git
```

### Step 2 — Copy the Scripts

Copy both directories into the upstream `cerebellum` repository:

```bash
cp -r SZcerebellumModel/scripts_sh cerebellum/
cp -r SZcerebellumModel/scripts_py cerebellum/
```

Alternatively, move the directories:

```bash
mv SZcerebellumModel/scripts_sh SZcerebellumModel/scripts_py cerebellum/
```

### Step 3 — Verify the Directory Structure

After integration, the working directory should have the following structure:

```text
cerebellum/
├── cerebellum/
├── scripts_sh/
├── scripts_py/
└── ...
```

All workflow commands should be executed from the **root of the upstream `cerebellum` repository**, unless otherwise specified.

---

## 3. Running the Workflow

The repository includes a Bash script for executing an automated computational workflow covering:

1. Morphological degeneration
2. Network construction
3. Simulation configuration and execution
4. Post-simulation analysis

### Example Execution

Navigate to the upstream `cerebellum` directory:

```bash
cd cerebellum
```

Execute the provided workflow:

```bash
bash scripts_sh/run_test.sh
```

The `run_test.sh` script invokes `scripts_sh/sim_setup_v5.sh` with predefined simulation parameters.

The current example configuration uses:

- **Stimulation protocol:** `kase_stim_vitro`
- **Degeneration level:** 30% (`start=30`, `stop=40`, and `n_steps=1` select only the start value)
- **Degeneration steps:** 1
- **Synaptic alteration:** `glutamate:all`
- **Compensation flag:** `0`
- **DCN–IO flag:** `0`
- **Resilience parameters:** `k0 = 0.5`, `sharpness = 2`

These parameters can be modified directly in `scripts_sh/run_test.sh`.

**Note:** This workflow launches computational simulations rather than a lightweight unit test. Execution time and computational requirements depend on the selected network configuration, degeneration range, and simulation backend.

---

## 4. Reproducibility

**Validation status (23 October 2026):** BASAL, EVOKED and AUTOR smoke checks passed. A repeat run matched 1,398 numerical entries against the versioned smoke-test baselines with maximum absolute difference zero on the tested three-level dataset. This demonstrates repeatability of these tests in the validated environment, not independent reproduction of all manuscript figures, complete 16-level curves or sigma sweeps.

### Reproducing the simulations reported in the paper

To reproduce the main simulation conditions described in the associated manuscript, execute the following Bash scripts from the **root directory of the upstream `cerebellum` framework**, after completing the installation and integration steps above.

The experiments combine two stimulation protocols, **basal activity** and **Kase-inspired saccadic stimulation**, with two resilience configurations: **HR (with resilience)** and **NR (no resilience)**.

| Stimulation protocol | Resilience configuration | Script |
| --- | --- | --- |
| Basal activity | HR (with resilience) | [`scripts_sh/batch_setups_HR_base.sh`](scripts_sh/batch_setups_HR_base.sh) |
| Basal activity | NR (no resilience) | [`scripts_sh/batch_setups_NR_base.sh`](scripts_sh/batch_setups_NR_base.sh) |
| Kase stimulation | HR (with resilience) | [`scripts_sh/batch_setups_HR_kase.sh`](scripts_sh/batch_setups_HR_kase.sh) |
| Kase stimulation | NR (no resilience) | [`scripts_sh/batch_setups_NR_kase.sh`](scripts_sh/batch_setups_NR_kase.sh) |

Run the simulations using:

```bash
bash scripts_sh/batch_setups_HR_base.sh
bash scripts_sh/batch_setups_NR_base.sh
bash scripts_sh/batch_setups_HR_kase.sh
bash scripts_sh/batch_setups_NR_kase.sh
```

Each launcher calls [`scripts_sh/sim_setup_v5.sh`](scripts_sh/sim_setup_v5.sh), which automatically performs morphological alterations, network construction, stimulus generation, simulation execution, and structural and functional analyses across the specified degeneration levels.

The HR and NR configurations use different resilience parameters (`reserve_sharpness = 15` and `2`, respectively, with `reserve_k0 = 0.5`). The basal simulations assess spontaneous network activity, whereas the Kase simulations assess responses to saccadic-like mossy-fiber stimulation. The basal HR and NR scripts cover 0–40% degeneration in 32 steps; Kase HR uses 32 steps and Kase NR uses 16 steps over the same range.

These launchers reproduce the main **simulation setups**. Exact numerical replication additionally depends on using the same random seeds, simulation environment, and analysis settings as in the original experiments. The current pipeline generates seeds at runtime, so repeated executions may differ numerically.

### Reproducibility Archive

The source code for network construction, simulation, and analysis is publicly available in this GitHub repository.

A dedicated **Zenodo reproducibility archive** is currently being prepared and will be made publicly available shortly.

The archive is intended to provide:

- Simulation configurations and parameter settings
- Numerical datasets underlying the main and supplementary figures
- Additional scripts and metadata required to reproduce the reported results
- Supporting materials for independent verification of the computational analyses

**Reproducibility Archive (Zenodo):** https://doi.org/10.5281/zenodo.23391996

The reproducibility archive will be maintained as a separate record from the manuscript preprint.

Until the archive is released, users can access the source code through this repository and consult the associated manuscript for a detailed description of the computational methods.

---

## 5. Preprint and Citation

The scientific methodology, computational framework, and main findings associated with this repository are described in the following preprint, publicly available on Zenodo:

**Vergani, A. A., Faris, P., De Grazia, M., Fusar-Poli, P., Casellato, C., & D'Angelo, E. U. (2026).** *Cerebellar Degeneration Induces Timeliness Disruption, Paradoxical Information Flow and Network Overdrive in a Digital Model of Schizophrenia* (Version 1.0). Zenodo.

**DOI:** https://doi.org/10.5281/zenodo.18813588

**Preprint:** [Access on Zenodo](https://zenodo.org/records/18813588)

### How to Cite

If you use this computational framework, its implementation, or its associated methodology in your research, please cite the preprint above.

**BibTeX:**

```bibtex
@misc{Vergani2026SZcerebellum,
  author = {Vergani, Alberto Arturo and Faris, Pawan and De Grazia, Marialaura and Fusar-Poli, Paolo and Casellato, Claudia and D'Angelo, Egidio Ugo},
  title = {Cerebellar Degeneration Induces Timeliness Disruption, Paradoxical Information Flow and Network Overdrive in a Digital Model of Schizophrenia},
  year = {2026},
  publisher = {Zenodo},
  version = {1.0},
  doi = {10.5281/zenodo.18813588},
  url = {https://doi.org/10.5281/zenodo.18813588},
  note = {Preprint}
}
```

---

## 6. Acknowledgments

This work was supported by the **Virtual Brain Twin project**, which received funding from the European Union's Horizon Europe research and innovation programme under Grant Agreement No. **101137289**.

Project information: [Virtual Brain Twin — CORDIS](https://cordis.europa.eu/project/id/101137289)

---

## 7. Contact

**Alberto Arturo Vergani, PhD**

Department of Brain and Behavioral Sciences  
University of Pavia, Italy

Email: [albertoarturo.vergani@unipv.it](mailto:albertoarturo.vergani@unipv.it)

GitHub: [albertoarturovergani](https://github.com/albertoarturovergani)
