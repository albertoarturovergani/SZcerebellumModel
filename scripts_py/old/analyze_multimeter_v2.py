#!/usr/bin/env python3

import gc
import neo
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import yaml
from collections import defaultdict
import argparse
import re
import pandas as pd
from scipy.stats import sem


def extract_atrophy_level(path):
    match = re.search(r"_deg_(\d+\.\d+)", str(path))
    if match:
        ratio = float(match.group(1))
        return round((1 - ratio) * 100, 1)  # percentuale di atrofia
    return None


def find_peak_features(signal, time, win):
    mask = (time >= win[0]) & (time <= win[1])
    if np.any(mask):
        segment = signal[mask]
        t_segment = time[mask]
        if len(segment) > 0:
            idx_min = np.argmin(segment)
            return segment[idx_min], t_segment[idx_min]
    return np.nan, np.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Cartella principale della simulazione")
    args = parser.parse_args()

    t_pre = 0
    t_post = 6
    Z_MOhm = 1.0
    root_dir = Path(args.root)

    summary = []
    all_lfp_by_atrophy = defaultdict(list)
    time_epoch = None

    for result_dir in sorted(root_dir.glob('*_deg_*/3.function')):
        print(f"\n📁 Entrata in {result_dir}")
        atrophy = extract_atrophy_level(result_dir)

        stim_dir = result_dir.parent / "2.stimulus"
        yaml_files = list(stim_dir.glob("*.yaml"))
        assert len(yaml_files) == 1
        with open(yaml_files[0], 'r') as f:
            stim_dict = yaml.safe_load(f)

        devices = stim_dict.get("simulations", {}).get("pawan_stim_vitro", {}).get("devices", {})
        burst_events = [dev["start"] for name, dev in devices.items()
                        if name.startswith("burst_") and isinstance(dev, dict) and "start" in dev]
        burst_events = np.array(sorted(burst_events))
        print(f"✅ {len(burst_events)} eventi trovati")

        epoched_lfp = []

        for file in sorted(result_dir.glob("*.nio")):
            reader = neo.io.NixIO(filename=str(file), mode="ro")
            block = reader.read_block()
            for segment in block.segments:
                analogs = {a.name: a for a in segment.analogsignals if a.name.lower().startswith("i_syn")}
                if not analogs:
                    continue
                fs = 1.0 / list(analogs.values())[0].sampling_period.rescale("s").magnitude
                t = list(analogs.values())[0].times.rescale("ms").magnitude
                n_pre = int(t_pre / 1000 * fs)
                n_post = int(t_post / 1000 * fs)
                for a in analogs.values():
                    s = a.rescale("pA").magnitude.squeeze()
                    lfp = -s * Z_MOhm / 1000.0
                    for ev in burst_events:
                        idx = np.searchsorted(t, ev)
                        if idx - n_pre >= 0 and idx + n_post < len(s):
                            epoched_lfp.append(lfp[idx - n_pre: idx + n_post])
            reader.close()
            del reader, block
            gc.collect()

        if not epoched_lfp:
            continue

        stack_lfp = np.stack(epoched_lfp)
        if time_epoch is None:
            time_epoch = np.linspace(-t_pre, t_post, stack_lfp.shape[1])

        all_lfp_by_atrophy[atrophy].append(stack_lfp)

        for trial in stack_lfp:
            amp_N2a, lat_N2a = find_peak_features(trial, time_epoch, (2.0, 3.0))
            amp_N2b, lat_N2b = find_peak_features(trial, time_epoch, (3.0, 4.0))
            summary.append({
                "atrophy": atrophy,
                "amp_N2a": amp_N2a,
                "lat_N2a": lat_N2a,
                "amp_N2b": amp_N2b,
                "lat_N2b": lat_N2b
            })

    df = pd.DataFrame(summary)

    fig, axs = plt.subplots(2, 2, figsize=(12, 10))

    from matplotlib.cm import get_cmap
    from matplotlib.colors import Normalize
    
    # Colormap viridis e normalizzatore
    atrophies_sorted = sorted(all_lfp_by_atrophy.keys())
    norm = Normalize(vmin=min(atrophies_sorted), vmax=max(atrophies_sorted))
    cmap = get_cmap("viridis_r")
    
    for a in atrophies_sorted:
        lfp_concat = np.concatenate(all_lfp_by_atrophy[a], axis=0)
        mean_lfp = np.nanmean(lfp_concat, axis=0)
        color = cmap(norm(a))
        axs[0, 0].plot(time_epoch, mean_lfp, label=f"{a:.1f}%", color=color)


    axs[0, 0].set_title("LFP Media per Livello di Atrofia")
    axs[0, 0].set_xlabel("Tempo (ms)")
    axs[0, 0].set_ylabel("LFP (µV)")
    axs[0, 0].legend(title="Atrofia")

    df.groupby("atrophy")[["amp_N2a", "amp_N2b"]].mean().plot(kind='bar', ax=axs[0, 1])
    axs[0, 1].set_title("Ampiezza Media N2a/N2b")
    axs[0, 1].set_ylabel("Ampiezza (µV)")

    df.groupby("atrophy")[["lat_N2a", "lat_N2b"]].mean().plot(kind='bar', ax=axs[1, 0])
    axs[1, 0].set_title("Latenza Media N2a/N2b")
    axs[1, 0].set_ylabel("Latenza (ms)")

    df_count = df.copy()
    df_count['has_N2a'] = df_count['amp_N2a'].notna()
    df_count['has_N2b'] = df_count['amp_N2b'].notna()
    df_count = df_count.groupby("atrophy")[["has_N2a", "has_N2b"]].sum()
    df_count.plot(kind='bar', ax=axs[1, 1])
    axs[1, 1].set_title("Conteggio Picchi N2a/N2b")

    plt.tight_layout()
    plt.savefig(root_dir / "summary_N2_peaks.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()