#!/usr/bin/env python3
"""Validate SZ BASAL 3-level bridge PKLs and run only the ORIGINAL delta-FR formula.

No NEST, reconstruction, historical notebook side effects, or TE calculations.
Does not attempt the paper's 16-point regression or figures needing morphometric DCI.

Example:
  python validate_SZ_basal_bridge_v1.py \
    --bridge /mnt/bsb_repro/SZ_NOTEBOOK_BRIDGE_SMOKE_V2 \
    --notebook /mnt/bsb_repro/SZcerebellumModel-TEST/notebooks/historical/META_PAPER_MODELING_BASAL_PREFIXED.ipynb \
    --output /mnt/bsb_repro/SZ_NOTEBOOK_BASAL_VALIDATION_V1
"""
import argparse
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd

REQUIRED = {
    "structure": ("case", "k", "connection", "synapses", "weight", "presyn_pop", "firing_rate_mean", "dfType", "SZ Factor (%)"),
    "weights": ("case", "k", "connection", "weight"),
    "functional": ("case", "k", "pop", "firing_rate_mean", "firing_rate_std", "firing_rate_sem", "dfType", "SZ Factor (%)"),
}


def read_frame(path, keys, required):
    if not path.is_file():
        raise FileNotFoundError(path)
    df = pd.read_pickle(path)
    if df.empty:
        raise ValueError(f"Empty DataFrame: {path}")
    absent = sorted(set(required) - set(df.columns))
    if absent:
        raise ValueError(f"Missing required columns in {path}: {absent}")
    if df.duplicated(keys).any():
        raise ValueError(f"Duplicate key rows in {path}: {keys}")
    df = df.copy()
    df["k"] = pd.to_numeric(df["k"], errors="raise").round(3)
    return df


def load_historical_delta_function(path):
    """Extract only the original BASAL function definition, not notebook cells."""
    nb = json.loads(path.read_text(encoding="utf-8"))
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        if "def compute_delta_fr_with_percent_and_zscore" not in source:
            continue
        tree = ast.parse(source)
        definitions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "compute_delta_fr_with_percent_and_zscore"]
        if len(definitions) != 1:
            continue
        only_definition = ast.fix_missing_locations(ast.Module(body=definitions, type_ignores=[]))
        namespace = {"np": np, "pd": pd}
        exec(compile(only_definition, str(path), "exec"), namespace)
        return namespace["compute_delta_fr_with_percent_and_zscore"]
    raise ValueError(f"Original BASAL delta-FR function not found in {path}")


def audit(bridge, notebook):
    structure = read_frame(bridge / "basal_df_all_str.pkl", ("case", "k", "connection"), REQUIRED["structure"])
    weights = read_frame(bridge / "basal_df_all_weights.pkl", ("case", "k", "connection"), REQUIRED["weights"])
    functional = read_frame(bridge / "basal_df_all_fr.pkl", ("case", "k", "pop"), REQUIRED["functional"])

    cases = set(structure["case"].unique())
    if cases != {"HR_basal4hz", "NR_basal4hz"}:
        raise ValueError(f"Unexpected cases: {sorted(cases)}")
    keys_s = set(zip(structure.case, structure.k))
    keys_f = set(zip(functional.case, functional.k))
    if keys_s != keys_f:
        raise ValueError(f"Inconsistent case/k coverage (structure={keys_s}, functional={keys_f})")
    if not set(zip(weights.case, weights.k)).issubset(keys_s):
        raise ValueError("Weights contain unknown case/k levels")

    # Verify duplicated weight values in the structural bridge match the YAML weight table.
    audit_weights = structure[["case", "k", "connection", "weight"]].merge(
        weights[["case", "k", "connection", "weight"]],
        how="left", on=["case", "k", "connection"],
        validate="one_to_one", suffixes=("_bridge", "_yaml"), indicator=True)
    has_yaml = audit_weights["_merge"].eq("both")
    equal_values = np.isclose(
        audit_weights.loc[has_yaml, "weight_bridge"].astype(float),
        audit_weights.loc[has_yaml, "weight_yaml"].astype(float),
        equal_nan=True)
    if not bool(equal_values.all()):
        raise ValueError("The structural bridge does not match YAML synaptic weights")
    missing_weight_connections = sorted(structure.loc[structure.weight.isna(), "connection"].dropna().unique().tolist())
    missing_presyn = structure.loc[structure.firing_rate_mean.isna(), ["connection", "presyn_pop"]].drop_duplicates()

    # Apply the ORIGINAL notebook's delta-FR function to each case independently,
    # as the notebook originally processes two experiments separately.
    original = load_historical_delta_function(notebook)
    parts = []
    for case, g in functional.groupby("case", sort=True):
        if 1.000 not in g.k.values:
            raise ValueError(f"Missing healthy baseline k=1.000 for {case}")
        parts.append(original(g.copy()))
    deltas = pd.concat(parts, ignore_index=True)
    expected_delta_columns = {"delta_fr_mean", "delta_fr_sem", "delta_fr_percent", "zscore_fr"}
    if not expected_delta_columns.issubset(deltas.columns):
        raise ValueError("Historical delta function did not produce the required columns")
    baseline = deltas.loc[deltas.k == 1.000, "delta_fr_mean"].dropna()
    if not np.allclose(baseline, 0., atol=1e-9):
        raise ValueError("Delta-FR values at k=1.000 are not zero")

    report = {
        "status": "PASS",
        "scope": "BASAL bridge integrity and original delta-FR formula only, not original notebook/figures",
        "cases": sorted(cases),
        "k_levels": sorted(float(x) for x in structure.k.unique()),
        "rows": {"structure": len(structure), "weights": len(weights), "functional": len(functional)},
        "matched_yaml_weight_rows": int(has_yaml.sum()),
        "without_chemical_yaml_weight_connections": missing_weight_connections,
        "without_recorded_presynaptic_fr_rows": int(structure.firing_rate_mean.isna().sum()),
        "without_recorded_presynaptic_fr_connections": [
            {"connection": str(x.connection), "presyn_pop": (None if pd.isna(x.presyn_pop) else str(x.presyn_pop))}
            for x in missing_presyn.itertuples(index=False)
        ],
        "missing_morphometric_fields": [x for x in ("DCI", "cell_name") if x not in structure.columns],
        "historical_delta_formula_used": "compute_delta_fr_with_percent_and_zscore",
        "method_notes": [
            "Electrical gap_goc has no chemical YAML weight; leave missing.",
            "IO presynaptic firing rate was not recorded; leave missing.",
            "DCI and cell_name are not inferred from synaptic connectivity.",
            "No fit of degree five on a three-level smoke test.",
            "Historical transfer entropy is not changed or re-estimated.",
        ],
    }
    return deltas, report


def plots(deltas, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for pop in ("purkinje", "dcn_p"):
        sub = deltas[deltas["pop"].eq(pop)]
        if sub.empty:
            print("SKIP (no population):", pop)
            continue
        fig, ax = plt.subplots(figsize=(7, 4))
        for case, group in sub.groupby("case"):
            ordered = group.sort_values("SZ Factor (%)")
            ax.plot(ordered["SZ Factor (%)"], ordered["delta_fr_mean"], marker="o", label=case)
        ax.axhline(0, linewidth=0.8, linestyle="--")
        ax.set_xlabel("SZ factor (%)")
        ax.set_ylabel("Δ firing rate relative to k=1.000 (Hz)")
        ax.set_title(f"BASAL {pop}: three-level DIAGNOSTIC, not a paper figure")
        ax.legend()
        fig.tight_layout()
        p = output / f"smoke_basal_delta_fr_{pop}.png"
        fig.savefig(p, dpi=180)
        plt.close(fig)
        print("PLOT:", p)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bridge", required=True, type=Path)
    p.add_argument("--notebook", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    bridge = args.bridge.expanduser().resolve()
    notebook = args.notebook.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if not bridge.is_dir():
        p.error(f"Bridge directory not found: {bridge}")
    if not notebook.is_file():
        p.error(f"Historical notebook not found: {notebook}")
    if output == bridge or output in bridge.parents or bridge in output.parents:
        p.error("Output must be separate from bridge input")
    deltas, report = audit(bridge, notebook)
    print("BASAL BRIDGE CHECK: PASS")
    print("Cases:", ', '.join(report['cases']))
    print("K levels:", report["k_levels"])
    print("Rows:", report["rows"])
    print("YAML weights matched:", report["matched_yaml_weight_rows"])
    print("Chemical weights unavailable:", report["without_chemical_yaml_weight_connections"])
    print("Presynaptic FR unavailable:", report["without_recorded_presynaptic_fr_rows"])
    print("DCI/cell_name unavailable:", report["missing_morphometric_fields"])
    if args.dry_run:
        print("DRY RUN: no files written")
        return
    if output.exists() and any(output.iterdir()):
        p.error(f"Refusing to write into a nonempty output directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    deltas.to_csv(output / "basal_delta_firing_smoke.csv", index=False)
    (output / "basal_bridge_check.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    plots(deltas, output)
    print("BASAL DIAGNOSTIC: PASS (original notebook not executed)")

if __name__ == "__main__":
    main()
