#!/usr/bin/env python3
"""Non-destructive Auto-Rhythm three-k diagnostic using the historical notebook FR method.

Reads PKLs generated from the replay; extracts only the historical function
compute_delta_fr_with_percent_and_zscore from the supplied original notebook,
then uses it separately for HR and NR. It does not execute other notebook cells.
No TE, BSB, NEST or HDF5/NIO data are changed.
"""
import argparse
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd

CASES = {
    "HR_auto": "20251120193505_basal_activity_15_0_glutamate_all",
    "NR_auto": "20251120172014_basal_activity_2_0_glutamate_all",
}
EXPECTED_K = {1.000, 0.760, 0.600}
REQUIRED = {"k", "pop", "firing_rate_mean", "firing_rate_std", "firing_rate_sem"}
ORIGINAL_FUNCTION = "compute_delta_fr_with_percent_and_zscore"


def load_original_function(notebook: Path):
    doc = json.loads(notebook.read_text(encoding="utf-8"))
    candidates = []
    for cell in doc.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        source = cell.get("source", "")
        source = "".join(source) if isinstance(source, list) else source
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == ORIGINAL_FUNCTION:
                candidates.append(node)
    if len(candidates) != 1:
        raise ValueError(f"Expected one {ORIGINAL_FUNCTION} definition, found {len(candidates)} in {notebook}")
    isolated = ast.fix_missing_locations(ast.Module(body=[candidates[0]], type_ignores=[]))
    namespace = {"pd": pd, "np": np}
    exec(compile(isolated, str(notebook), "exec"), namespace)
    return namespace[ORIGINAL_FUNCTION]


def load_case(input_root: Path, label: str, folder: str, delta_function):
    summary = input_root / folder / "B_collective_results"
    structural = summary / "1.structures/df_summary.pkl"
    conn = summary / "1.structures/summary_connections.pkl"
    func = summary / "2.functions/df_summary.pkl"
    for f in (structural, conn, func):
        if not f.is_file():
            raise FileNotFoundError(f"{label}: missing {f}")
    fs = pd.read_pickle(func)
    ds = pd.read_pickle(structural)
    cs = pd.read_pickle(conn)
    for name, df in (("functional", fs), ("structure", ds), ("connections", cs)):
        if df.empty:
            raise ValueError(f"{label}: empty {name} frame")
        if "k" not in df.columns:
            raise ValueError(f"{label}: missing k in {name}")
        ks = {round(float(x), 3) for x in df["k"].dropna().unique()}
        if ks != EXPECTED_K:
            raise ValueError(f"{label}: {name} has k={ks}, expected {EXPECTED_K}")
    missing = sorted(REQUIRED - set(fs.columns))
    if missing:
        raise ValueError(f"{label}: missing functional columns {missing}")
    if fs.duplicated(["k", "pop"]).any():
        raise ValueError(f"{label}: duplicate (k,pop) functional rows")
    if fs[["k", "firing_rate_mean"]].isna().any().any():
        raise ValueError(f"{label}: NaN in k or firing_rate_mean")

    # The historical notebook calls this function on each experiment BEFORE
    # merging the HR and NR frames. Keep that order to avoid mixing baselines.
    result = delta_function(fs)
    result["case"] = label
    result["resilience"] = "With Resilience" if label.startswith("HR") else "Without Resilience"
    result["SZ Factor (%)"] = 100.0 * (1.0 - result["k"].astype(float))
    for target in ("purkinje", "dcn_p"):
        sub = result[result["pop"] == target]
        assert {round(float(v), 3) for v in sub["k"]} == EXPECTED_K, f"{label}: {target} incomplete"
        v = sub.loc[np.isclose(sub["k"], 1.0), "delta_fr_mean"].to_numpy(dtype=float)
        if len(v) != 1 or not np.isclose(v[0], 0):
            raise ValueError(f"{label}: baseline for {target} is not delta=0: {v}")
    no_baseline = sorted(set(result["pop"]) - set(result.loc[np.isclose(result["k"], 1.), "pop"]))
    zero_baseline = sorted(result.loc[np.isclose(result["k"], 1.) & (result["firing_rate_mean"] == 0), "pop"].tolist())
    print(f"PASS {label}: structure={ds.shape}, connections={cs.shape}, functional={fs.shape}")
    print("   k:", sorted(EXPECTED_K), "populations:", result["pop"].nunique())
    if no_baseline:
        print("   Populations lacking k=1 baseline:", no_baseline)
    if zero_baseline:
        print("   Zero-baseline populations (percent delta is undefined):", zero_baseline)
    return result


def make_plots(df, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for pop in ("purkinje", "dcn_p"):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for case, group in df.loc[df["pop"] == pop].groupby("case", sort=True):
            group = group.sort_values("SZ Factor (%)")
            ax.plot(group["SZ Factor (%)"], group["delta_fr_mean"], "o-", label=case)
        ax.axhline(0, color="gray", linestyle="--", linewidth=1)
        ax.set(xlabel="Simulated degeneration (%)", ylabel="Delta mean firing rate (Hz)",
               title=f"Auto-Rhythm: {pop} firing-rate change (three-point smoke test)")
        ax.legend()
        fig.tight_layout()
        target = out / f"autor_smoke_delta_fr_{pop}.png"
        fig.savefig(target, dpi=180)
        plt.close(fig)
        print("PLOT:", target)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=Path, help=".../legacy_analysis from clean postprocessing")
    ap.add_argument("--notebook", required=True, type=Path, help="Original AUTOR notebook")
    ap.add_argument("--output", required=True, type=Path, help="Diagnostic output directory")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    try:
        historic_delta = load_original_function(args.notebook)
        print("HISTORICAL FR FUNCTION FOUND:", ORIGINAL_FUNCTION)
        rows = [load_case(args.input, label, folder, historic_delta) for label, folder in CASES.items()]
        df = pd.concat(rows, ignore_index=True)
        print("TOTAL:", df.shape, "rows across HR and NR")
        print("=== PURKINJE / DCN-P PREVIEW ===")
        cols = ["case", "pop", "k", "firing_rate_mean", "delta_fr_mean", "delta_fr_percent"]
        print(df.loc[df["pop"].isin(["purkinje", "dcn_p"]), cols].sort_values(["case", "pop", "k"], ascending=[True, True, False]).to_string(index=False))
        print("AUTOR DIAGNOSTIC CHECK: PASS")
        if args.dry_run:
            print("DRY RUN: no files written.")
            return
        args.output.mkdir(parents=True, exist_ok=True)
        csv_path = args.output / "autor_smoke_delta_firing.csv"
        df.to_csv(csv_path, index=False)
        report = {
            "status": "PASS", "kind": "three-level diagnostic, NOT manuscript reproduction",
            "source_notebook": str(args.notebook),
            "source_function": ORIGINAL_FUNCTION,
            "cases": sorted(CASES), "k": sorted(EXPECTED_K),
            "rows": len(df), "metrics": ["firing_rate_mean", "delta_fr_mean", "delta_fr_percent", "zscore_fr"],
            "notes": ["Missing and zero-baseline populations are not imputed.",
                      "No fifth-degree polynomial fit is applied.",
                      "Original notebook is not executed in full or modified."]
        }
        (args.output / "autor_smoke_check.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        make_plots(df, args.output)
        print("CSV:", csv_path)
        print("AUTOR DIAGNOSTIC: PASS (original notebook not executed)")
    except Exception as e:
        ap.exit(1, f"AUTOR DIAGNOSTIC: FAIL: {type(e).__name__}: {e}\n")

if __name__ == "__main__":
    main()
