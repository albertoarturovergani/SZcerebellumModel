#!/usr/bin/env python3
import os
import sys
import yaml
import matplotlib.pyplot as plt
import numpy as np
import re
from collections import defaultdict

import os
import glob
import yaml
import matplotlib.pyplot as plt

def generate_event_aligned_histogram_from_yaml(stimuli, bin_width=5, window=(-200, 200)):
    bursts = defaultdict(list)
    for stim in stimuli:
        match = re.match(r"(pre|post)_burst_(\d+)_\d+", stim["label"])
        if match:
            bursts[int(match.group(2))].append(stim)

    t_min, t_max = window
    time_bins = np.arange(t_min, t_max, bin_width)
    firing_rate = np.zeros(len(time_bins))

    for stim_list in bursts.values():
        t0_list = [s["start"] for s in stim_list if "post_burst" in s["label"]]
        if not t0_list:
            continue
        t0 = min(t0_list)

        for s in stim_list:
            rel_start = s["start"] - t0
            rel_stop = s["stop"] - t0
            bin_indices = (time_bins >= rel_start) & (time_bins < rel_stop)
            firing_rate[bin_indices] += s["rate"]

    firing_rate /= max(len(bursts), 1)
    return time_bins + bin_width / 2, firing_rate

def generate_stimulus_plots(sim_dir: str):
    if not os.path.isdir(sim_dir):
        print(f"❌ Directory non trovata: {sim_dir}")
        return

    yaml_files = [f for f in os.listdir(sim_dir) if f.endswith(".yaml")]
    if not yaml_files:
        print(f"⚠️ Nessun file YAML trovato in: {sim_dir}")
        return

    stimuli = []
    for yfile in yaml_files:
        with open(os.path.join(sim_dir, yfile), 'r') as f:
            ydata = yaml.safe_load(f)
            try:
                devices = ydata['simulations']['basal_activity']['devices']
                for label, stim in devices.items():
                    if stim.get("device") == "poisson_generator":
                        stimuli.append({
                            "label": label,
                            "start": stim["start"],
                            "stop": stim["stop"],
                            "rate": stim["rate"],
                            "type": "burst" if "burst" in label else "stim"
                        })
            except Exception:
                continue

    if not stimuli:
        print(f"⚠️ Nessuno stimolo valido trovato in {sim_dir}")
        return

    stimuli.sort(key=lambda x: x["start"])
    stim_dir = os.path.join(sim_dir, "stimulus")
    os.makedirs(stim_dir, exist_ok=True)

    # === Plot 2: Rate vs Time ===
    midpoints = [(s["start"] + s["stop"]) / 2 for s in stimuli]
    rates = [s["rate"] for s in stimuli]
    plt.figure(figsize=(8, 5))
    plt.scatter(midpoints, rates, marker="o", color="darkgreen")
    plt.xlabel("Time (ms)")
    plt.ylabel("Rate (Hz)")
    plt.title("Rate vs Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(stim_dir, "stimulus_rate_vs_time.png"))
    plt.close()
    print(f"✅ Salvato: stimulus_rate_vs_time.png")

    # === Plot 4: Istogramma peri-saccadico reale ===
    time_bins, rate_profile = generate_event_aligned_histogram_from_yaml(stimuli)
    plt.figure(figsize=(10, 5))
    plt.bar(time_bins, rate_profile, width=5, align='center', color='gray', edgecolor='black')
    plt.axvline(0, color='red', linestyle='--', label='Saccade onset')
    plt.title("Istogramma peri-saccadico (da YAML)")
    plt.xlabel("Tempo rispetto alla saccade (ms)")
    plt.ylabel("Firing rate medio (Hz)")
    plt.grid(True, linestyle=':')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(stim_dir, "peri_saccadic_histogram.png"))
    plt.close()
    print(f"✅ Salvato: peri_saccadic_histogram.png")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_stim_timeline_v2.py <sim_dir>")
        sys.exit(1)

    sim_dir = sys.argv[1]
    generate_stimulus_plots(sim_dir)
