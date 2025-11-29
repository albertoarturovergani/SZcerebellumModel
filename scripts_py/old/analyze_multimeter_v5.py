#!/usr/bin/env python3

import gc
import neo
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import yaml
from collections import defaultdict
import argparse
from scipy.stats import sem
import re
from scipy.signal import butter, filtfilt, find_peaks
import matplotlib.colors as mcolors
import pandas as pd

lowcut = 0.1
highcut = 300
order = 4

Z_Ohm = 10.0 / 2  # Impedenza in MOhm

def bandpass_filter(data, fs, lowcut, highcut, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return filtfilt(b, a, data)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Cartella principale della simulazione")
    args = parser.parse_args()

    t_pre = 0.1
    t_post = 15
    root_dir = Path(args.root)

    for result_dir in sorted(root_dir.glob('*_deg_*/3.function')):
        print(f"\n\U0001F4C1 Entrata in {result_dir}")
        output_dir = result_dir / "plots_multimeter"
        output_dir.mkdir(exist_ok=True)

        stim_dir = result_dir.parent / "2.stimulus"
        yaml_files = list(stim_dir.glob("*.yaml"))
        assert len(yaml_files) == 1, f"❌ Attesi 1 solo file YAML, trovati {len(yaml_files)} in {stim_dir}"
        stim_yaml_path = yaml_files[0]

        with open(stim_yaml_path, 'r') as f:
            stim_dict = yaml.safe_load(f)

        devices = stim_dict.get("simulations", {}).get("pawan_stim_vitro", {}).get("devices", {})
        burst_events = sorted([dev["start"] for name, dev in devices.items() if name.startswith("burst_") and isinstance(dev, dict) and "start" in dev])
        burst_events = np.array(burst_events)
        print(f"✅ Trovati {len(burst_events)} eventi burst da {stim_yaml_path.name}")

        epoched_signals = defaultdict(list)
        fs = None

        for file in sorted(result_dir.glob("*.nio")):
            print(f"\n\U0001F4C2 Processing {file.name}")
            reader = neo.io.NixIO(filename=str(file), mode="ro")
            block = reader.read_block()

            for segment in block.segments:
                print(f"  ↪️ Segment: {segment.name}")
                analogs = {a.name: a for a in segment.analogsignals if a.name.lower().startswith("i_syn")}
                if not analogs:
                    print(f"⚠️ Nessun segnale i_syn nel segmento {segment.name} di {file.name}")
                    continue

                if fs is None:
                    fs = 1.0 / list(analogs.values())[0].sampling_period.rescale("s").magnitude
                    t = list(analogs.values())[0].times.rescale("ms").magnitude
                    n_pre = int(t_pre / 1000 * fs)
                    n_post = int(t_post / 1000 * fs)

                for k, a in analogs.items():
                    s = a.rescale("pA").magnitude.squeeze()
                    for ev in burst_events:
                        if t[0] <= ev <= t[-1]:
                            idx = np.searchsorted(t, ev)
                            if idx - n_pre >= 0 and idx + n_post < len(s):
                                epoch = s[idx - n_pre: idx + n_post]
                                baseline = np.mean(epoch[:n_pre])
                                epoch = epoch - baseline
                                epoched_signals[k].append(epoch)

            reader.close()
            del reader, block
            gc.collect()

        if fs is None:
            print("❌ Nessun segnale i_syn trovato in nessun file .nio. Salto questa cartella.")
            continue

        time_epoch = np.linspace(-t_pre, t_post, n_pre + n_post, endpoint=False)
        np.save(output_dir / "epoched_signals.npy", epoched_signals)

        all_trials_sum = []
        for k in sorted(epoched_signals):
            stack = np.stack(epoched_signals[k])
            if stack.size == 0:
                continue
            mean = np.mean(stack, axis=0)
            plt.figure(figsize=(6, 4))
            for trial in stack:
                lfp = -trial / (4 * np.pi * Z_Ohm) / 1e6
                plt.plot(time_epoch, lfp, color='gray', alpha=0.3)
                all_trials_sum.append(lfp)
            mean_lfp = -mean / (4 * np.pi * Z_Ohm) / 1e6
            plt.plot(time_epoch, mean_lfp, label=f"{k} (mean)", color='black', lw=2)
            plt.title(f"LFP butterfly (cellulare) - {k}")
            plt.xlabel("Tempo (ms)")
            plt.ylabel("LFP (V)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig(output_dir / f"{k}_lfp_butterfly_cells.png", dpi=150)
            plt.close()

        # Plot della somma totale (butterfly globale)
        if all_trials_sum:
            all_trials_sum = np.stack(all_trials_sum)
            total_mean = np.mean(all_trials_sum, axis=0)
            plt.figure(figsize=(6, 4))
            for trial in all_trials_sum:
                plt.plot(time_epoch, trial, color='gray', alpha=0.3)
            plt.plot(time_epoch, total_mean, color='black', lw=2, label='mean')
            plt.title("LFP butterfly - somma totale")
            plt.xlabel("Tempo (ms)")
            plt.ylabel("LFP (V)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig(output_dir / "lfp_butterfly_sum_all_cells.png", dpi=150)
            plt.close()

if __name__ == "__main__":
    main()
