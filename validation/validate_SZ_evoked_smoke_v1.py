#!/usr/bin/env python3
"""Diagnose EVOKED (Kase) HR/NR smoke-test results from regenerated PKLs.

This is NOT a full reproduction of the historical EVOKED notebook. It reads
existing sigma-specific Kase intermediate data and optionally checks the
BASAL bridge files needed in a later historical notebook stage. All scientific
values, particularly the historical transfer-entropy estimator, are preserved.

WARNING: Only read trusted .pkl files; Python pickle can execute code.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

CASES = {
    "HR_kase": "20251113214858_kase_stim_vitro_15_0_glutamate_all",
    "NR_kase": "20251113161643_kase_stim_vitro_2_0_glutamate_all",
}
COLS = ["k", "optimal_lag_ms", "te_mossy_fibers_to_purk", "te_purk_to_mossy_fibers"]
BASAL_NAMES = ["basal_df_all_str.pkl", "basal_df_all_weights.pkl", "basal_df_all_fr.pkl"]


def get_k_values(df: pd.DataFrame) -> set[float]:
    if "k" not in df.columns:
        raise ValueError("Colonna 'k' mancante")
    ks = pd.to_numeric(df["k"], errors="raise").to_numpy(dtype=float)
    if not np.isfinite(ks).all():
        raise ValueError("Valori 'k' non finiti")
    return {round(float(k), 3) for k in ks}


def read_inputs(root: Path, sigma: int, expected: set[float], bridge: Path | None):
    dfs = []
    audit = {"sigma": sigma, "cases": {}, "basal_bridge": None}
    for case, folder in CASES.items():
        p = root / folder / "B_collective_results" / "2.functions" / f"sigma_{sigma}" / f"{sigma}_df.pkl"
        if not p.is_file():
            raise FileNotFoundError(f"File EVOKED mancante: {p}")
        df = pd.read_pickle(p)
        if df.empty:
            raise ValueError(f"PKL EVOKED vuoto: {p}")
        missing = sorted(set(COLS) - set(df.columns))
        if missing:
            raise ValueError(f"Colonne mancanti in {p}: {missing}")
        observed = get_k_values(df)
        if observed != expected:
            raise ValueError(f"Livelli k inattesi per {case}: {observed}, attesi {expected}")
        if df["k"].duplicated().any():
            raise ValueError(f"Livelli k duplicati per {case} in {p}")
        value_block = df[COLS[1:]].apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
        if not np.isfinite(value_block).all():
            raise ValueError(f"Lag o TE non finite per {case}")
        if 1.000 not in observed:
            raise ValueError(f"Manca il baseline k=1.000 per {case}")

        out = pd.DataFrame({
            "case": case,
            "k": df["k"].to_numpy(dtype=float),
            "sz_percent": (1.0 - df["k"].to_numpy(dtype=float)) * 100.0,
            "optimal_lag_ms": df["optimal_lag_ms"].to_numpy(dtype=float),
            "te_mf_pc_legacy": df["te_mossy_fibers_to_purk"].to_numpy(dtype=float),
            "te_pc_mf_legacy": df["te_purk_to_mossy_fibers"].to_numpy(dtype=float),
        })
        baseline = float(out.loc[np.isclose(out["k"], 1), "te_mf_pc_legacy"].iloc[0])
        out["delta_te_mf_pc_legacy_vs_k1"] = out["te_mf_pc_legacy"] - baseline
        out["te_directionality_legacy"] = out["te_mf_pc_legacy"] - out["te_pc_mf_legacy"]
        if "burst_onsets_ms" in df.columns:
            out["n_bursts"] = df["burst_onsets_ms"].apply(
                lambda x: len(x) if isinstance(x, (list, tuple, np.ndarray)) else np.nan
            ).to_numpy()
        dfs.append(out)
        audit["cases"][case] = {
            "file": str(p), "rows": len(df), "k_values": sorted(observed, reverse=True),
            "te_baseline_legacy": baseline,
            "columns_verified": COLS,
        }
        print(f"PASS {case}: {len(df)} livelli k, lag e TE finite; baseline TE={baseline:.6f}")

    if bridge is not None:
        files = {}
        for filename in BASAL_NAMES:
            p = bridge / filename
            if not p.is_file():
                raise FileNotFoundError(f"Bridge BASAL mancante: {p}")
            df = pd.read_pickle(p)
            if df.empty:
                raise ValueError(f"Bridge BASAL vuoto: {p}")
            found = get_k_values(df)
            if found != expected:
                raise ValueError(f"Livelli k nel BASAL bridge {filename}: {found}, attesi {expected}")
            files[filename] = {"rows": len(df), "columns": len(df.columns), "k_values": sorted(found)}
            print(f"PASS BASAL bridge {filename}: {df.shape}")
        audit["basal_bridge"] = files

    result = pd.concat(dfs, ignore_index=True).sort_values(["case", "k"], ascending=[True, False])
    return result, audit


def plot_series(data: pd.DataFrame, column: str, ylabel: str, title: str, target: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    for label, group in data.groupby("case", sort=True):
        g = group.sort_values("sz_percent")
        ax.plot(g["sz_percent"], g[column], marker="o", label=label)
    ax.set_xlabel("Degenerazione SZ (%)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(target, dpi=220)
    plt.close(fig)
    print("PLOT:", target)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True,
                    help="Directory legacy_analysis con i PKL Kase rigenerati")
    ap.add_argument("--output", type=Path, required=True, help="Cartella OUTPUT separata")
    ap.add_argument("--basal-bridge", type=Path, help="Cartella con i tre PKL bridge BASAL")
    ap.add_argument("--sigma", type=int, default=10, help="Sigma diagnostico già calcolato (default 10)")
    ap.add_argument("--levels", default="1.000,0.760,0.600",
                    help="Tre livelli attesi separati da virgole")
    ap.add_argument("--dry-run", action="store_true", help="Controlla, senza scrivere file")
    a = ap.parse_args()
    expected = {round(float(s.strip()), 3) for s in a.levels.split(",")}
    if len(expected) != 3:
        ap.error("Lo smoke test EVOKED richiede esattamente 3 livelli")
    if a.sigma < 0:
        ap.error("--sigma deve essere >=0")

    data, audit = read_inputs(a.input, a.sigma, expected, a.basal_bridge)
    audit["status"] = "PASS"
    audit["note"] = (
        "Diagnostic only. Uses historical TE numbers as stored. No TE recalculation. "
        "One sigma and three k levels do not replicate paper figures or sigma selection. "
        "BASAL bridge presence does not imply that the original EVOKED notebook runs unchanged."
    )
    print("\n=== EVOKED DIAGNOSTIC PREVIEW ===")
    print(data[["case", "k", "optimal_lag_ms", "te_mf_pc_legacy", "te_pc_mf_legacy"]].to_string(index=False))
    if a.dry_run:
        print("\nEVOKED BRIDGE CHECK: PASS (dry run; nessun file scritto)")
        return

    a.output.mkdir(parents=True, exist_ok=True)
    data.to_csv(a.output / "evoked_smoke_metrics.csv", index=False)
    with (a.output / "evoked_bridge_check.json").open("w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=2, ensure_ascii=False)
    specs = [
        ("optimal_lag_ms", "Lag MF→PC (ms)", "Kase: lag MF→PC", "evoked_smoke_lag.png"),
        ("te_mf_pc_legacy", "TE MF→PC (scala storica)", "Kase: TE legacy MF→PC", "evoked_smoke_te_legacy.png"),
        ("delta_te_mf_pc_legacy_vs_k1", "ΔTE legacy vs k=1.000", "Kase: variazione TE vs baseline", "evoked_smoke_delta_te_legacy.png"),
        ("te_directionality_legacy", "TE(MF→PC) − TE(PC→MF), scala storica", "Kase: direzionalità TE legacy", "evoked_smoke_te_directionality.png"),
    ]
    for col, ylabel, title, filename in specs:
        plot_series(data, col, ylabel, title, a.output / filename)
    print("\nEVOKED DIAGNOSTIC: PASS. Dati originali e TE storica invariati.")


if __name__ == "__main__":
    main()
