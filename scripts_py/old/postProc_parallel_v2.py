import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from collections import defaultdict, Counter
from pathlib import Path
from fooof import FOOOF
import matplotlib.pyplot as plt
from math import log2
import pickle
import gc
import neo
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict
from scipy.signal import convolve
from scipy.signal.windows import triang
from pylatex import Document, Section, Figure, NoEscape
import argparse
import re
import yaml
import numpy as np
from collections import Counter
from math import log2
from fooof import FOOOF
from scipy.signal import welch
from scipy.signal import hilbert
from scipy.stats import pearsonr
import seaborn as sns
import argparse

def plot_time_series(df_all, pop="mossy", dt=0.1):
    plt.figure(figsize=(10, 5))
    for cond in df_all["condition"].unique():
        sub = df_all[df_all["condition"] == cond]
        for i, row in sub.iterrows():
            k = row["k"]
            t = np.arange(len(row[f"{pop}_timeseries"])) * dt
            y = np.array(row[f"{pop}_timeseries"])
            label = f"{cond} (k={k:.2f})"
            plt.plot(t, y, label=label, alpha=0.6)

    plt.title(f"Firing rate istantaneo medio – {pop.capitalize()}")
    plt.xlabel("Time (ms)")
    plt.ylabel("Firing rate (Hz)")
    plt.legend(loc="upper right", bbox_to_anchor=(1.25, 1.0))
    plt.grid(True)
    plt.tight_layout()
    plt.show()



def compute_firing_rate(spike_times, duration, dt=1.0, sigma=100):
    time = np.arange(0, duration, dt)
    spike_train = np.zeros_like(time)
    indices = (np.array(spike_times) / dt).astype(int)
    indices = indices[indices < len(spike_train)]
    spike_train[indices] = 1
    kernel_size = max(int(2 * sigma / dt), 3)
    kernel = triang(kernel_size)
    kernel /= kernel.sum()
    rate = convolve(spike_train, kernel, mode='same') * (1000.0 / dt)
    return time, rate

def compute_population_psd(array, fs):
    if array.size == 0: return [], []
    psds = [welch(unit, fs=fs, nperseg=1024)[1] for unit in array]
    return welch(array[0], fs=fs, nperseg=1024)[0].tolist(), np.mean(psds, axis=0).tolist()

def analyze_fooof(signal, fs, freq_range=(1, 100)):
    freqs, psd = welch(signal, fs=fs, nperseg=1024)
    fm = FOOOF(peak_width_limits=[1, 12], max_n_peaks=5, verbose=False)
    fm.fit(freqs, psd, freq_range)
    return {
        "aperiodic_params": fm.aperiodic_params_.tolist(),
        "peak_params": [p.tolist() for p in fm.peak_params_],
        "r_squared": float(fm.r_squared_),
        "error": float(fm.error_),
        "n_peaks": int(fm.n_peaks_)
    }

def compute_cross_correlation(x, y):
    x = (x - np.mean(x)) / np.std(x)
    y = (y - np.mean(y)) / np.std(y)
    return float(np.max(np.abs(np.correlate(x, y, mode='full')) / len(x)))

def compute_optimal_lag(x, y, dt=1.0, max_lag_ms=200):
    max_lag_steps = int(max_lag_ms / dt)
    x, y = (x - np.mean(x)) / np.std(x), (y - np.mean(y)) / np.std(y)
    corr = np.correlate(y, x, mode='full') / len(x)
    lags = np.arange(-len(x)+1, len(x)) * dt
    center = len(corr) // 2
    window = slice(center - max_lag_steps, center + max_lag_steps + 1)
    best_idx = np.argmax(np.abs(corr[window]))
    return float(lags[window][best_idx]), float(corr[window][best_idx])

def compute_plv(sig1, sig2):
    phase1, phase2 = np.angle(hilbert(sig1)), np.angle(hilbert(sig2))
    return float(np.abs(np.mean(np.exp(1j * (phase1 - phase2)))))

def mean_pairwise_corr(array):
    n = array.shape[0]
    return float(np.nanmean([pearsonr(array[i], array[j])[0]
                             for i in range(n) for j in range(i+1, n)])) if n > 1 else np.nan

def entropy(probs):
    return -sum(p * log2(p) for p in probs if p > 0)

def conditional_entropy(p_joint, p_conditioned):
    total_ce = 0.0
    for joint_key, p_joint_val in p_joint.items():
        if not isinstance(joint_key, tuple): continue
        cond_key = joint_key[1] if len(joint_key) == 2 else (joint_key[1], joint_key[2]) if len(joint_key) == 3 else None
        if cond_key in p_conditioned and p_joint_val > 0:
            total_ce += p_joint_val * log2(p_conditioned[cond_key] / p_joint_val)
    return total_ce

def compute_transfer_entropy_binned(source, target, bins=16):
    s_d = np.digitize(source, np.histogram(source, bins=bins)[1]) - 1
    t_d = np.digitize(target, np.histogram(target, bins=bins)[1]) - 1
    x_t, y_t, y_t1 = s_d[:-1], t_d[:-1], t_d[1:]
    p_y1y, p_y = Counter(zip(y_t1, y_t)), Counter(y_t)
    p_y1yx, p_yx = Counter(zip(y_t1, y_t, x_t)), Counter(zip(y_t, x_t))
    total = len(y_t1)
    for d in [p_y1y, p_y, p_y1yx, p_yx]: d.update({k: v/total for k, v in d.items()})
    return conditional_entropy(p_y1y, p_y) - conditional_entropy(p_y1yx, p_yx)



from analyze_worker import analyze_single_path

def loadData(data_path: Path, condition_name: str, use_top_quartile: bool = False,
             selection_mode: str = "percentile", selection_percentile: float = 0.99) -> pd.DataFrame:
    dt = 0.1
    sigma = 100
    freq_range = (1, 100)
    trim_start = 0.01
    trim_end = 0.15

    summary_rows = []
    from concurrent.futures import ProcessPoolExecutor
    paths = list(data_path.rglob("results_*/basal_activity"))
    tasks = [
        (path, condition_name, 0.1, 100, (1,100), 0.01, 0.15, selection_percentile)
        for path in paths if path.parts[-1] == "basal_activity"
    ]
    with ProcessPoolExecutor(max_workers=4) as executor:
        for result in executor.map(analyze_single_path, tasks):
            if result:
                summary_rows.append(result)
    # fine parallelizzazione
        match = [part for part in path.parts if part.startswith("results_")]
        if not match:
            print(f"⚠️ No k found in: {path}")
            continue
        k = float(match[0].replace("results_", ""))
        atrophy = 100 * (1 - k)
        print(f"📁 {condition_name} – Analisi per k={k:.3f}")

        # Carica eventi burst mossy dal file YAML associato (es. per k=1.000)
        burst_onsets_ms, burst_durations_ms = [], []
        yaml_file = path.parent / "reconstruction" / f"basal_vitro_{k:.3f}.yaml"
        if yaml_file.exists():
            import yaml
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

        # Per ogni popolazione, calcola metriche e aggiungi al row
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

        # Connettività tra tutte le coppie di popolazioni
        pop_names = list(pop_rates.keys())
        for i in range(len(pop_names)):
            for j in range(i + 1, len(pop_names)):
                a, b = pop_names[i], pop_names[j]
                a_mean, b_mean = pop_rates[a], pop_rates[b]
                row[f"cosine_similarity_{a}_{b}"] = float(np.dot(a_mean, b_mean) / (np.linalg.norm(a_mean) * np.linalg.norm(b_mean)))
                row[f"cross_corr_{a}_{b}"] = compute_cross_correlation(a_mean, b_mean)
                row[f"optimal_lag_ms_{a}_{b}"] = compute_optimal_lag(a_mean, b_mean, dt)[0]
                row[f"te_{a}_to_{b}"] = compute_transfer_entropy_binned(a_mean, b_mean)
                row[f"te_{b}_to_{a}"] = compute_transfer_entropy_binned(b_mean, a_mean)
                row[f"plv_{a}_{b}"] = compute_plv(a_mean, b_mean)
                
        # Metriche complessive di attività della rete
        if pop_rates:
            all_activity = np.vstack(list(pop_rates.values())).mean(axis=0)
            row["network_mean_rate"] = float(all_activity.mean())
            row["network_quantile_rate"] = float(np.quantile(all_activity, selection_percentile))
            row["network_sem_rate"] = float(np.std(all_activity) / np.sqrt(len(all_activity)))
            row["network_timeseries"] = all_activity.tolist()

        summary_rows.append(row)


    df_summary = pd.DataFrame(summary_rows)
    if not df_summary.empty:
        df_summary = df_summary.sort_values("k")

        # Salvataggio automatico in postProcessing
        postproc_path = Path(f"{Path}/postProcessing")
        postproc_path.mkdir(exist_ok=True)
        df_summary.to_pickle(postproc_path / f"df_summary_{condition_name}.pkl")

    return df_summary

def main():
    parser = argparse.ArgumentParser(description="Post-processing cerebellum analysis")
    parser.add_argument("folder", type=str, help="Path to main data folder")
    parser.add_argument("--condition", type=str, default="cond", help="Condition name")
    args = parser.parse_args()

    df = loadData(Path(args.folder), condition_name=args.condition)

    if args.plot:
        plot_time_series(df, pop="mossy")
        plot_time_series(df, pop="purkinje")

if __name__ == "__main__":
    
