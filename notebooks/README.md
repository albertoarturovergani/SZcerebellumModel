# Jupyter notebooks for the SZ cerebellum model

## Fast, independent check of newly generated simulations

Use [`quick_smoke_test.ipynb`](quick_smoke_test.ipynb). It reads **fresh HDF5 + 8 NIO MPI partitions** produced by `reproducibility/scripts/replay.py` directly, without requiring historical result caches.

Recommended workflow (run from repository root, with the NEST/BSB Python venv activated):

```bash
python reproducibility/scripts/replay.py --case all --mpi 8 --output "$HOME/SZ_replay_outputs"
python reproducibility/scripts/analyze_quick_results.py --input "$HOME/SZ_replay_outputs" --strict
python -m pip install -r reproducibility/environment/requirements-notebooks.txt
python -m jupyter lab notebooks/quick_smoke_test.ipynb
```

The default replay selects `k=1.000`, `0.760`, `0.600` per case (18 simulations total). `--full` reproduces all 16 levels per case (96 simulations total). `--levels 1.000,0.760,0.600` explicitly selects any archived subset.

Output data are in `~/SZ_replay_outputs/analysis_smoke/` by default. The smoke analysis produces raw population firing rates, with HR/NR Purkinje and DCN-P plots for the three protocols. It intentionally does **not** claim to regenerate the exact Figure 2–4 analyses or transfer entropy: those depend on the original statistical processing and full levels.

## Historical manuscript analysis notebooks

The author's three notebooks from the [Zenodo reproducibility record](https://zenodo.org/records/23391996) are staged in [`historical/`](historical/) by:

```bash
python reproducibility/scripts/import_historical_notebooks.py
```

They are located in the local 2025 archive `~/cerebellum/cerebellum/ref_sims/paperModelingSZ1_reproducibility/` and retain **unchanged code cell sources**, with execution outputs cleared. Original SHA256 hashes and a dependency scan are stored in `historical/NOTEBOOK_PROVENANCE.json`.

| Original notebook | Historical experiment pairs | Manuscript context |
| --- | --- | --- |
| `META_PAPER_MODELING_AUTOR_PREFIXED.ipynb` | HR/NR 0 Hz runs from November 2025 | Auto-Rhythm and supplementary figures |
| `META_PAPER_MODELING_BASAL_PREFIXED.ipynb` | HR/NR 4 Hz runs from June/October 2025 | Figure 2 and structural summaries |
| `META_PAPER_MODELING_EVOKED_v5_PREFIXED.ipynb` | HR/NR Kase runs from November 2025 | Figures 3 and 4 |

## Rebuild the original analysis intermediates from new HDF5/NIO

The old workflow **does generate** `B_collective_results/1.structures/df_summary.pkl`
(using `do_network_analysis.py`) and `B_collective_results/2.functions/df_summary.pkl`
(using `do_completeAnalysis_v2.py`). Its Kase analysis generates
`B_collective_results/2.functions/sigma_*/<sigma>_df.pkl` using
`do_completeAnalysis_vbtkase_v2.py`. The BASAL *notebook*, not BSB/NEST alone,
creates `tempPlots/basal_df_all_str.pkl` and `tempPlots/basal_df_all_weights.pkl`;
EVOKED reads those BASAL notebook outputs.

`replay.py` produces HDF5 and NIO only. The adapter `rebuild_intermediates.py`
stages a **separate analysis workspace** matching the historical layout, with
`1.structure/`, `2.stimulus/`, and `3.function/` directories under the original run
names, plus per-experiment `B_collective_results/` summaries. It links, rather
than copies, NIO/HDF5 and copies configuration YAML so legacy processors can edit
only their working copies.

After the 18-simulation smoke test and after copying the original notebooks:

```bash
python reproducibility/scripts/import_historical_notebooks.py
python reproducibility/scripts/rebuild_intermediates.py \
  --input "$HOME/SZ_replay_outputs" --require-smoke
python reproducibility/scripts/rebuild_intermediates.py \
  --input "$HOME/SZ_replay_outputs" --require-smoke --execute \
  --sigma-baseline 100 --sigma-kase 10
```

The first adapter invocation checks and stages input files; the `--execute` command
runs the **original** post-processing for the available levels, generating its
`.pkl` summaries under `~/SZ_replay_outputs/legacy_analysis/<original-run-name>/`.
Logs are created per run. This execution may expose additional incompatibilities
with modern packages: it has been statically and synthetically checked, but has
not yet been verified against a fresh independent 18-simulation dataset.

For the full scientific analyses, use **all 16 levels** per case and specify the
Kase sigma grid used in the paper, for example:

```bash
python reproducibility/scripts/replay.py --case all --full \
  --output "$HOME/SZ_replay_outputs" --mpi 8
python reproducibility/scripts/rebuild_intermediates.py \
  --input "$HOME/SZ_replay_outputs" --execute \
  --sigma-kase 10,25,50,75,100
```

Enter the legacy workspace and start Jupyter *there*:

```bash
cd "$HOME/SZ_replay_outputs/legacy_analysis"
python -m jupyter lab
```

Open `META_PAPER_MODELING_BASAL_PREFIXED.ipynb` first (it creates the
`tempPlots/basal_*.pkl` files), then AUTOR and finally EVOKED. The original
notebooks retain their historical code and may need controlled adaptations to
run under a new environment or with only three levels. Their historical schemas
and plotting settings must be checked before claiming Figure 2–4 equivalence.

You can alternatively download and analyze the original 2025 record from
[Zenodo 23391996](https://zenodo.org/records/23391996). Never invent, interpolate,
or replace unavailable historical NIOs with new-replay recordings. In particular,
the original HR Basal 4 Hz run at `k=0.733` had no matching archived NIO, while
the 2026 replay successfully generated its own recording.
