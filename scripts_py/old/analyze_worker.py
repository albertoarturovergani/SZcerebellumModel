
import numpy as np
import pandas as pd
from scipy.signal import convolve, welch, hilbert
from scipy.signal.windows import triang
from fooof import FOOOF
from math import log2
from pathlib import Path
from collections import defaultdict, Counter
from scipy.stats import pearsonr
import neo
import gc
import yaml

def analyze_single_path(args):
    (path, condition_name, dt, sigma, freq_range, trim_start, trim_end, selection_percentile) = args
    try:
        match = [part for part in path.parts if part.startswith("results_")]
        if not match:
            print(f"⚠️ No k found in: {path}")
            return None
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
                    start = float(props["start"])
                    stop = float(props["stop"])
                    burst_onsets_ms.append(start)
                    burst_durations_ms.append(stop - start)

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

        if not df.empty:
            # Seleziona solo il 10% randomico dei granuli
            granule_ids = df[df["pop"] == "granule"]["sender_id"].unique()
            if len(granule_ids) > 0:
                np.random.seed(42)  # per riproducibilità
                keep_ids = set(np.random.choice(granule_ids, size=max(1, int(0.1 * len(granule_ids))), replace=False))
                df = df[(df["pop"] != "granule") | (df["sender_id"].isin(keep_ids))]

        if df.empty:
            print("⚠️ No data after filtering")
            return None

        duration = df["time_ms"].max()
        grouped = df.groupby(["pop", "sender_id"])

        data = defaultdict(list)
        for (pop, sid), g in grouped:
            t = np.arange(0, duration, dt)
            spike_train = np.zeros_like(t)
            indices = (np.array(g["time_ms"].values) / dt).astype(int)
            indices = indices[indices < len(spike_train)]
            spike_train[indices] = 1
            kernel_size = max(int(2 * sigma / dt), 3)
            kernel = triang(kernel_size)
            kernel /= kernel.sum()
            rate = convolve(spike_train, kernel, mode='same') * (1000.0 / dt)
            data["pop"].append(pop)
            data["sender_id"].append(sid)
            data["rate"].append(rate)
            data["time"].append(t)

        n = len(data["time"][0])
        start = int(n * trim_start)
        end = int(n * (1 - trim_end))

        row = {
            "condition": condition_name,
            "k": k,
            "atrophy_percent": atrophy,
            "burst_onsets_ms": burst_onsets_ms,
            "burst_durations_ms": burst_durations_ms,
        }

        pop_rates = {}
        for pop in sorted(set(data["pop"])):
            rates = np.vstack([r for i, r in enumerate(data["rate"]) if data["pop"][i] == pop])
            if rates.size == 0:
                continue
            mean_rate = rates[:, start:end].mean(0)
            sem_rate = rates[:, start:end].std(0) / np.sqrt(rates.shape[0])
            quantile_value = float(np.quantile(mean_rate, selection_percentile))

            row[f"{pop}_mean_rate"] = float(mean_rate.mean())
            row[f"{pop}_quantile_rate"] = quantile_value
            row[f"{pop}_sem_rate"] = float(sem_rate.mean())
            row[f"{pop}_timeseries"] = mean_rate.tolist()

            pop_rates[pop] = mean_rate

        pop_names = list(pop_rates.keys())
        for i in range(len(pop_names)):
            for j in range(i + 1, len(pop_names)):
                a, b = pop_rates[pop_names[i]], pop_rates[pop_names[j]]
                row[f"cross_corr_{i}_{j}"] = float(np.max(np.abs(np.correlate((a - np.mean(a))/np.std(a), (b - np.mean(b))/np.std(b), mode='full')) / len(a)))

        if pop_rates:
            all_activity = np.vstack(list(pop_rates.values())).mean(axis=0)
            row["network_mean_rate"] = float(all_activity.mean())
            row["network_quantile_rate"] = float(np.quantile(all_activity, selection_percentile))
            row["network_sem_rate"] = float(np.std(all_activity) / np.sqrt(len(all_activity)))
            row["network_timeseries"] = all_activity.tolist()

        gc.collect()
        return row
    except Exception as e:
        print(f"❌ Error while processing {path}: {e}")
        return None
