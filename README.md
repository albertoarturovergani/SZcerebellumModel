# SZcerebellumModel

Code and reproducibility inputs for **Cerebellar Degeneration Induces Timeliness Disruption,
Paradoxical Information Flow and Network Overdrive in a Digital Model of Schizophrenia**.

This repository contains the historical source implementation of the cerebellum model,
configuration files, archived morphology inputs and Python/Bash analysis/generation scripts.
Large outputs (HDF5/NIO) are hosted separately: **[Zenodo record 23391996](https://zenodo.org/records/23391996)**.

## Archived scientific experiments

Six historical experiments, 16 degeneration levels each (96 network reconstructions plus NEST simulations):

| Case | Original 2025 experiment | Protocol | Resilience | Levels |
|---|---|---|---|---:|
| HR_auto | `20251120193505_basal_activity_15_0_glutamate_all` | auto | HR | 16 |
| NR_auto | `20251120172014_basal_activity_2_0_glutamate_all` | auto | NR | 16 |
| HR_basal4hz | `20250628122534_basal_activity_15_0_glutamate_all` | basal4hz | HR | 16 |
| NR_basal4hz | `20251004180523_basal_activity_2_0_glutamate_all` | basal4hz | NR | 16 |
| HR_kase | `20251113214858_kase_stim_vitro_15_0_glutamate_all` | kase | HR | 16 |
| NR_kase | `20251113161643_kase_stim_vitro_2_0_glutamate_all` | kase | NR | 16 |

`k=1.000` denotes no degeneration; `k=0.600` denotes 40% degeneration.

- **Auto-Rhythm / 0 Hz:** Poisson input to mossy fibers set to 0 Hz, weight 0; 5 s.
- **Basal 4 Hz:** Poisson mossy-fiber input at 4 Hz, weight 1; 5 s.
- **Kase:** 4 Hz mossy-fiber baseline and explicitly archived saccadic burst times/rates; 10 s.

The two Basal experiment families are different. Auto-Rhythm is not the 4 Hz Basal protocol.

## Software and randomness

The scientific package corresponds to **cerebellum (cerebellar-models) v0.5.0**.
The 2026 replay was performed with NEST 3.7 and a BSB 4.4-family distribution in Python 3.11.
The exact Python package inventory was captured during staging, if the correct environment was active.
Review `reproducibility/environment/python_package_inventory.json` before citing a complete software lock.

**Structural generation:** archived BSB configurations do not contain a global structural RNG seed,
so independent stochastic neuron placement and connectivity will not be bitwise identical.
BSB 8.0 introduces a reproducibility-oriented global RNG configuration, but this does not
retroactively determine the unrecorded structural random state of the BSB 4.4 experiments.

**NEST simulations:** archived Basal/Auto YAML use seed `1234`; Kase YAML use individual saved seeds.
A fixed NEST seed does not imply the same spikes when the structural network or environment differs.
Newer cerebellar-models releases, including 0.12.1, were **not** used for the historical replay.

The upstream `v0.5.0` structural YAML mistakenly declared `cerebellum~=0.4.0`. The operational
GitHub inputs update **only that dependency metadata** to `cerebellum~=0.5.0`, while preserving
scientific parameters. Changes and file hashes are documented under `reproducibility/manifests/`.

## Quick start: recreate networks and NEST outputs

1. Use Linux/Ubuntu, Python 3.11, MPI, NEST 3.7 with matching PyNEST, and a compiled `cerebmodule`.
   See **[Ubuntu environment instructions](reproducibility/environment/INSTALL_UBUNTU.md)**.
2. Create/activate a virtual environment and install the package from this checkout:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e .
python -m pip check
```

3. **Default: three-level smoke test.** This runs exactly `k=1.000`, `0.760`, and `0.600` for each of the six cases (18 simulations, not 96). Keep output files outside the Git checkout:

```bash
python reproducibility/scripts/replay.py --case all --output "$HOME/SZ_replay_outputs" --mpi 8
```

For just one case, replace `--case all` with `--case NR_basal4hz` (3 simulations).
For a different subset, use `--levels 1.000,0.760,0.600`.
For a **configuration-only check**, add `--prepare-only` (no HDF5/NIO simulation).

4. Generate diagnostic results from the new HDF5/NIO files:

```bash
python reproducibility/scripts/analyze_quick_results.py \
  --input "$HOME/SZ_replay_outputs" --strict
```

The script writes population firing rates and three diagnostic plots to `~/SZ_replay_outputs/analysis_smoke/`.
Open [`notebooks/quick_smoke_test.ipynb`](notebooks/quick_smoke_test.ipynb) to explore the same results interactively.

5. **Full published parameter sweep (explicit opt-in):** all 16 levels in all six cases, 96 simulations:

```bash
python reproducibility/scripts/replay.py --case all --full --output "$HOME/SZ_replay_outputs" --mpi 8
```

A default three-level test is sufficient for environment/pipeline verification but **not** for reproducing complete manuscript curves or their inferential statistics.

**The smoke test must be run on each target installation; package installation does not build NEST itself.**
The replay script copies configurations and a *working* set of morphologies and source code to a
separate output directory before calling `bsb compile`, `bsb reconfigure`, and `bsb simulate`.
It never edits the archived inputs or overwrites an existing unfinished run automatically.

## Reproduce analyses and explore notebooks

### New simulations (fast diagnostic)

[`notebooks/quick_smoke_test.ipynb`](notebooks/quick_smoke_test.ipynb) and
[`reproducibility/scripts/analyze_quick_results.py`](reproducibility/scripts/analyze_quick_results.py)
read newly generated HDF5/NIO directly, combine the **eight MPI NIO partitions** correctly,
and plot raw mean firing rates for Purkinje and DCN-P in each experimental protocol.
They require no previous `stats.yaml` or historical `.pkl` cache.
The 3-point outputs are useful for **installation verification**, not a substitute for
reproducing the full published analyses.

### Original 2025 manuscript analysis notebooks

Three real Jupyter notebooks used in the study are included under
[`notebooks/historical/`](notebooks/historical/) from the 2025
[Zenodo reproducibility archive](https://zenodo.org/records/23391996):

| Notebook | Experiment pair | Purpose |
| --- | --- | --- |
| `META_PAPER_MODELING_AUTOR_PREFIXED.ipynb` | `HR_auto` / `NR_auto` | Auto-Rhythm and supplements |
| `META_PAPER_MODELING_BASAL_PREFIXED.ipynb` | `HR_basal4hz` / `NR_basal4hz` | Basal responses and Fig. 2 |
| `META_PAPER_MODELING_EVOKED_v5_PREFIXED.ipynb` | `HR_kase` / `NR_kase` | Evoked Fig. 3–4 analyses |

The original notebook inputs include `B_collective_results/1.structures/df_summary.pkl`
and `B_collective_results/2.functions/df_summary.pkl`. The **original**
`scripts_py/do_network_analysis.py` and `scripts_py/do_completeAnalysis_v2.py`
generate these files; the Kase processor generates `sigma_*/<sigma>_df.pkl`.
The historical BASAL notebook itself creates `tempPlots/basal_df_all_str.pkl` and
`tempPlots/basal_df_all_weights.pkl`, which EVOKED subsequently reads.

**These intermediates can be rebuilt; they are not mandatory downloads.** The portable
`replay.py` stops after HDF5/NIO. The optional adapter below makes an analysis-only
workspace with the historical `1.structure`, `2.stimulus`, `3.function` folder layout,
links the generated HDF5/NIO without duplicating them, and calls the original processing
scripts to recreate intermediate `.pkl` files:

```bash
python reproducibility/scripts/rebuild_intermediates.py \
  --input "$HOME/SZ_replay_outputs" --require-smoke
python reproducibility/scripts/rebuild_intermediates.py \
  --input "$HOME/SZ_replay_outputs" --require-smoke --execute
```

The first command only stages the layout; the second executes legacy analysis. It
writes to `$HOME/SZ_replay_outputs/legacy_analysis/`, never modifying the repository's
archived configurations or the original Zenodo datasets. This legacy post-processing
path is **not yet validated end-to-end against a fresh 18-simulation NEST run**, so
inspect its logs and resulting DataFrame schemas before claiming exact reproduction
of the published plots.

The historical notebooks are copied into the legacy workspace if staged with
`python reproducibility/scripts/import_historical_notebooks.py`. Run BASAL there
before EVOKED so the `basal_*.pkl` dependencies are generated in the right order.
For the complete Fig. 2–4 analyses, run all 16 k levels plus the Kase sigma sweep
used in the original analysis. The 3-level smoke run is only for pipeline checks.
See [notebooks/README.md](notebooks/README.md) for exact steps and caveats.

Download original 2025 datasets and historical analysis resources from
[Zenodo record 23391996](https://zenodo.org/records/23391996).
Missing historical records must remain marked missing, never filled with new-simulation NIO.

## Repository layout

- `cerebellum/`: one local source snapshot of the v0.5.0 package (includes NESTML models).
- `configurations/SZ_repro_2025/experiments/`: 96 complete YAML configuration trees.
- `morphologies/SZ_repro_2025/`: 80 deduplicated SWCs plus morphology YAMLs.
- `reproducibility/`: portable runner (three-level default / full 16-level opt-in), raw-NIO diagnostic analysis, environment inventory and checksums.
- `notebooks/`: quick diagnostic notebook plus unchanged-code historical notebooks used for the paper.
- `data/`: instructions and ignore rules for optional data downloaded separately from Zenodo.
- `scripts_py/`, `scripts_sh/`: scientific generation and analysis pipelines from the SZ project.
- `validation/`: research validation materials; see availability notes.

The original scientific workflow scripts may generate new random Kase bursts;
**the archival replay above does not regenerate any burst or any YAML seed**.

## Source, data and citation

- **Cerebellar model code:** [dbbs-lab/cerebellar-models, version 0.5.0](https://github.com/dbbs-lab/cerebellar-models/tree/v0.5.0).
- **Brain Scaffold Builder (BSB):** [dbbs-lab/bsb](https://github.com/dbbs-lab/bsb).
- **Original simulation datasets (HDF5 and NIO):** [Zenodo reproducibility record 23391996](https://doi.org/10.5281/zenodo.23391996).
- **Scientific preprint:** [Zenodo preprint, version 1.0 (2026)](https://doi.org/10.5281/zenodo.18813588).

### Cite the scientific preprint

If you use this model, reproduce its computational experiments, or build upon its results, please cite:

**Vergani, A. A., Faris, P., De Grazia, M., Fusar-Poli, P., Casellato, C., & D'Angelo, E. U. (2026).** *Cerebellar Degeneration Induces Timeliness Disruption, Paradoxical Information Flow and Network Overdrive in a Digital Model of Schizophrenia* (Version 1.0) [Preprint]. Zenodo. https://doi.org/10.5281/zenodo.18813588

**BibTeX:**

```bibtex
@misc{Vergani2026SZcerebellum,
  author       = {Vergani, Alberto Arturo and Faris, Pawan and De Grazia, Marialaura and Fusar-Poli, Paolo and Casellato, Claudia and D'Angelo, Egidio Ugo},
  title        = {Cerebellar Degeneration Induces Timeliness Disruption, Paradoxical Information Flow and Network Overdrive in a Digital Model of Schizophrenia},
  year         = {2026},
  publisher    = {Zenodo},
  version      = {1.0},
  note         = {Preprint},
  doi          = {10.5281/zenodo.18813588},
  url          = {https://doi.org/10.5281/zenodo.18813588}
}
```

The [Zenodo dataset](https://doi.org/10.5281/zenodo.23391996) and the [scientific preprint](https://doi.org/10.5281/zenodo.18813588) are separate records. Cite the dataset when using archived original simulation outputs, in addition to citing the preprint when referring to the scientific study.

The source is distributed under this repository's GPL-3.0 license, with upstream authorship and licensing information in `reproducibility/UPSTREAM_NOTICE.md`.

## Funding and acknowledgments

The following funding acknowledgments concern the work of **Alberto Arturo Vergani (AAV)** and follow the attribution used in the manuscript's financial-disclosure draft:

- **EBRAINS 2.0:** This research received funding from the European Union's Horizon Europe Programme under Specific Grant Agreement No. **101147319** (*EBRAINS 2.0: A Research Infrastructure to Advance Neuroscience and Brain Health*; attributed in the manuscript to **Egidio D'Angelo (ED), Claudia Casellato (CC), and Alberto Arturo Vergani (AAV)**).
- **Virtual Brain Twin:** Work by **Alberto Arturo Vergani (AAV)** was supported by the European Union's Horizon Europe Research and Innovation Programme under Grant Agreement No. **101137289** (*Virtual Brain Twin for personalised treatment of Psychiatric Disorders*; project support attributed in the manuscript to **Egidio D'Angelo (ED)**).

European Commission project records: [EBRAINS 2.0 (101147319)](https://cordis.europa.eu/project/id/101147319) and [Virtual Brain Twin (101137289)](https://cordis.europa.eu/project/id/101137289).

Funding specific to other manuscript coauthors is detailed in the manuscript's Financial Disclosure. The final author/grant attribution should be consistent with the approved submission statement.

## Contact

Alberto Arturo Vergani — University of Pavia, Italy.

For questions about scientific reproducibility: albertoarturo.vergani@unipv.it.
