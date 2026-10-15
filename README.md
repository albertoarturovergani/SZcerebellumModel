
# SZcerebellumModel

A computational framework for investigating **schizophrenia-like structural and synaptic alterations** in a biologically grounded, spiking **cerebellar microcircuit model**.

This repository extends the [Brain Scaffold Builder (BSB)](https://github.com/bsb-team/bsb) and [cerebellum](https://github.com/bsb-team/cerebellum) computational frameworks to investigate how progressive neuronal degeneration and resilience mechanisms influence cerebellar network dynamics.

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

The repository provides two main directories that integrate with the existing `cerebellum` framework:

```text
SZcerebellumModel/
├── scripts_sh/
│   └── Shell scripts for workflow automation
├── scripts_py/
│   └── Python modules for network processing and analysis
└── README.md
```

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

The workflow builds upon the BSB and `cerebellum` frameworks and requires a compatible Python environment.

### 1.1 Brain Scaffold Builder (BSB)

Official repository: [bsb-team/bsb](https://github.com/bsb-team/bsb)

Install the core BSB components and the cerebellar extension:

```bash
python -m pip install bsb-core bsb-cereb
```

Additional dependencies may be required depending on the simulation backend and analysis scripts.

### 1.2 Cerebellum Framework

Official repository: [bsb-team/cerebellum](https://github.com/bsb-team/cerebellum)

Clone the upstream repository:

```bash
git clone https://github.com/bsb-team/cerebellum.git
```

For reproducing the analyses described in the associated manuscript, the reported software versions are BSB 4.4 and cerebellum module 0.4.0. Compatibility with newer releases should be verified.

---

## 2. Installation and Integration

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
- **Degeneration range:** 30–40%
- **Degeneration steps:** 1
- **Synaptic alteration:** `glutamate:all`
- **Compensation flag:** `0`
- **DCN–IO flag:** `0`
- **Resilience parameters:** `k0 = 0.5`, `sharpness = 2`

These parameters can be modified directly in `scripts_sh/run_test.sh`.

**Note:** This workflow launches computational simulations rather than a lightweight unit test. Execution time and computational requirements depend on the selected network configuration, degeneration range, and simulation backend.

---

## 4. Reproducibility

The source code for network construction, simulation, and analysis is publicly available in this GitHub repository.

A dedicated **Zenodo reproducibility archive** is currently being prepared and will be made publicly available shortly.

The archive is intended to provide:

- Simulation configurations and parameter settings
- Numerical datasets underlying the main and supplementary figures
- Additional scripts and metadata required to reproduce the reported results
- Supporting materials for independent verification of the computational analyses

**Reproducibility Archive (Zenodo):** *[DOI and repository link to be added upon release]*

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
