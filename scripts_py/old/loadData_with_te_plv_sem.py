
import numpy as np
import pandas as pd
from collections import defaultdict
from pathlib import Path
import neo
import gc
import yaml

from your_module import (
    compute_firing_rate,
    compute_cross_correlation,
    compute_optimal_lag,
    compute_transfer_entropy_binned,
    compute_plv,
    mean_pairwise_corr,
    analyze_fooof,
)

def sem_from_percentiles(values):
    percentiles = [0, 25, 50, 75, 100]
    extracted = np.percentile(values, percentiles)
    return float(np.std(extracted, ddof=1) / np.sqrt(len(extracted)))

def loadData(data_path: Path, condition_name: str, use_top_quartile: bool = False,
             selection_mode: str = "percentile", selection_percentile: float = 0.99) -> pd.DataFrame:
    dt = 0.1
    sigma = 100
    trim_start = 0.01
    trim_end = 0.15
    summary_rows = []

    for path in data_path.rglob("results_*/basal_activity"):
        match = [part for part in path.parts if part.startswith("results_")]
        if not match:
            print(f"⚠️ No k found in: {path}")
            continue
        k = float(match[0].replace("results_", ""))
        atrophy = 100 * (1 - k)
        print(f"📁 {condition_name} – Analisi per k={k:.3f}")

        burst_onsets_ms, burst_durations_ms = [], []
        yaml_file = path.parent / "reconstruction" / f"basal_vitro_{k:.3f}.yaml"
        if yaml_file.exists():
            with open(yaml_file, "r") as f:
                yaml_data = yaml.safe_load(f)
            devices = yaml_data.get("simulations", {}).get("basal_activity", {}).get("devices", {})
            for name, props in devices.items():
                if name.startswith("burst_") and props.get("device") == "poisson_generator":
                    burst_onsets_ms.append(float(props["start"]))
                    burst_durations_ms.append(float(props["stop"]) - float(props["start"]))

        spikes = []
        for f in sorted(path.glob("*.nio")):
            reader = neo.io.NixIO(filename=str(f), mode="ro")
            block = reader.read_block()
            for seg_idx, segment in enumerate(block.segments):
                for st_idx, st in enumerate(segment.spiketrains):
                    times = st.times.rescale("ms").magnitude
                    senders = st.annotations.get("senders", [f"unit{st_idx}"])
                    pop = st.annotations.get("device", "unknown").split("_")[0]
                    for time, sid in zip(times, senders):
                        spikes.append({
                            "time_ms": time,
                            "sender_id": sid,
                            "pop": pop,
                            "segment": seg_idx,
                            "file": f.name
                        })
            del reader, block
            gc.collect()

        df = pd.DataFrame(spikes)
        df = df[df['pop'].isin(['mossy', 'purkinje'])]
        if df.empty:
            print("⚠️ No data after filtering")
            continue

        duration = df["time_ms"].max()
        grouped = df.groupby(["pop", "sender_id"])

        data = defaultdict(list)
        for (pop, sid), g in grouped:
            t, r = compute_firing_rate(g["time_ms"].values, duration, dt=dt, sigma=sigma)
            data["pop"].append(pop)
            data["sender_id"].append(sid)
            data["rate"].append(r)
            data["time"].append(t)

        mossy = np.vstack([r for i, r in enumerate(data["rate"]) if data["pop"][i] == "mossy"])
        purk = np.vstack([r for i, r in enumerate(data["rate"]) if data["pop"][i] == "purkinje"])

        if mossy.size == 0 or purk.size == 0:
            print("⚠️ Dati insufficienti")
            continue

        n = len(data["time"][0])
        start = int(n * trim_start)
        end = int(n * (1 - trim_end))

        m_mean = mossy[:, start:end].mean(0)
        p_mean = purk[:, start:end].mean(0)

        te_m2p = np.array([compute_transfer_entropy_binned(m, p_mean) for m in mossy[:, start:end]])
        te_p2m = np.array([compute_transfer_entropy_binned(p, m_mean) for p in purk[:, start:end]])
        plv_m2p = np.array([compute_plv(m, p_mean) for m in mossy[:, start:end]])
        plv_p2m = np.array([compute_plv(m_mean, p) for p in purk[:, start:end]])
        plv_all = np.concatenate([plv_m2p, plv_p2m])

        row = {
            "condition": condition_name,
            "k": k,
            "atrophy_percent": atrophy,
            "mossy_mean_rate": float(m_mean.mean()),
            "mossy_quantile_rate": float(np.quantile(m_mean, selection_percentile)),
            "mossy_sem_rate": float(mossy[:, start:end].std(0).mean() / np.sqrt(mossy.shape[0])),
            "purkinje_mean_rate": float(p_mean.mean()),
            "purkinje_quantile_rate": float(np.quantile(p_mean, selection_percentile)),
            "purkinje_sem_rate": float(purk[:, start:end].std(0).mean() / np.sqrt(purk.shape[0])),
            "mossy_timeseries": m_mean.tolist(),
            "purkinje_timeseries": p_mean.tolist(),
            "cosine_similarity": float(np.dot(m_mean, p_mean) / (np.linalg.norm(m_mean) * np.linalg.norm(p_mean))),
            "cross_corr": compute_cross_correlation(m_mean, p_mean),
            "optimal_lag_ms": compute_optimal_lag(m_mean, p_mean, dt)[0],
            "te_mossy_to_purk": float(te_m2p.mean()),
            "te_mossy_to_purk_sem": sem_from_percentiles(te_m2p),
            "te_purk_to_mossy": float(te_p2m.mean()),
            "te_purk_to_mossy_sem": sem_from_percentiles(te_p2m),
            "plv": float(plv_all.mean()),
            "plv_sem": sem_from_percentiles(plv_all),
            "intra_corr_mossy": mean_pairwise_corr(mossy[:, start:end]),
            "intra_corr_purk": mean_pairwise_corr(purk[:, start:end]),
            "fooof_mossy": analyze_fooof(m_mean, fs=1000/dt),
            "fooof_purk": analyze_fooof(p_mean, fs=1000/dt),
            "burst_onsets_ms": burst_onsets_ms,
            "burst_durations_ms": burst_durations_ms
        }

        summary_rows.append(row)

    df_summary = pd.DataFrame(summary_rows)
    if not df_summary.empty:
        df_summary = df_summary.sort_values("k")
        postproc_path = Path("postProcessing")
        postproc_path.mkdir(exist_ok=True)
        df_summary.to_pickle(postproc_path / f"df_summary_{condition_name}.pkl")

    return df_summary
