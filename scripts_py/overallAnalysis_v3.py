#!/usr/bin/env python3

import os
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import argparse
import seaborn as sns
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
from matplotlib import colormaps
cmap = colormaps.get_cmap("viridis")

def flatten_dict(d, parent_key='', sep='__'):
    """Flatten a nested dictionary."""
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def plot_overlapped_psd(df_all, out_dir):
    """Plot overlapped PSD curves for mossy and Purkinje populations."""
    for pop in ['mossy', 'purkinje']:
        plt.figure(figsize=(10, 6))
        valid_k = df_all['k'].dropna().values
        norm = Normalize(vmin=np.min(valid_k), vmax=np.max(valid_k))

        for _, row in df_all.iterrows():
            k = row.get('k', None)
            freqs = row.get(f'psd__{pop}__frequencies_hz', None)
            power = row.get(f'psd__{pop}__power', None)

            if isinstance(freqs, list) and isinstance(power, list) and k is not None:
                color = cmap(norm(k))
                plt.plot(freqs, power, color=color, alpha=0.9)

        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power")
        plt.xlim(0, 100)
        plt.title(f"Overlapped PSD – {pop.capitalize()} population")
        plt.grid(True)

        # Sort legend by increasing k
        handles, labels = plt.gca().get_legend_handles_labels()
        if handles:
            sorted_labels_handles = sorted(zip(labels, handles), key=lambda x: float(x[0][2:]))
            labels, handles = zip(*sorted_labels_handles)
            plt.legend(handles, labels, title='k', bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.tight_layout()
        save_path = os.path.join(out_dir, f"overlapped_psd_{pop}.png")
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"📈 Saved overlapped PSD plot: {save_path}")

def analyze_yaml_plots(base_dir):
    pattern = os.path.join(base_dir, "*_deg_*/3.function/plots_firing_rate/summary_*.yaml")
    yaml_files = glob(pattern)

    print(f"\n📁 Found YAML files: {len(yaml_files)}")
    if not yaml_files:
        print("⚠️  No files found. Check path or folder structure.")
        return

    all_data = []
    for yfile in yaml_files:
        try:
            with open(yfile, 'r') as f:
                ydata = yaml.safe_load(f)
            flat = flatten_dict(ydata)

            basename = os.path.basename(yfile).replace('.yaml', '')
            parts = basename.split('_')

            #if len(parts) >= 3:
            #    flat['time_window'] = 'Overall Time'
            #    flat['k'] = float(parts[2])
            #else:
            #    print(f"⚠️ Unrecognized file: {basename}")

            all_data.append(flat)

        except Exception as e:
            print(f"⚠️ Error reading {yfile}: {e}")

    df_all = pd.DataFrame(all_data)
    df_all = df_all.apply(pd.to_numeric, errors='ignore')

    if df_all.empty or 'k' not in df_all.columns:
        print("⚠️  No valid data found in YAML files.")
        return

    out_dir = os.path.join(base_dir, "B_collective_results/2.functions")
    os.makedirs(out_dir, exist_ok=True)

    for window in sorted(df_all['time_window'].unique()):
        df_window = df_all[df_all['time_window'] == window]
        df_window = df_window.dropna(subset=['k'])

        print(f"\n📊 Analysis for {window} – {len(df_window)} rows")

        out_dir_window = out_dir  # Non creare sottocartella!

        for var in df_window.columns:
            if var in ['k', 'atrophy_percent', 'time_window']:
                continue

            valid_values = df_window[var].dropna()
            if valid_values.apply(lambda x: isinstance(x, (list, dict))).any():
                continue
            if valid_values.empty:
                continue

            plt.figure(figsize=(8, 5))
            sns.regplot(x='k', y=var, data=df_window, order=2)
            plt.xlabel('k')
            plt.ylabel(var)
            plt.title(f'{var} vs Atrophy – {window}')
            plt.grid(True)
            plt.tight_layout()
            filepath = os.path.join(out_dir_window, f"{var}_vs_atrophy.png")
            plt.savefig(filepath)
            plt.close()
            print(f"✅ Saved: {filepath}")

        # Save CSV summary
        df_window.to_csv(os.path.join(out_dir_window, f"summary_data.csv"), index=False)

        # Plot PSD only for Overall Time
        #if window == 'Overall Time':
        #    plot_overlapped_psd(df_window, out_dir_window)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze and plot YAML summary files (summary_k_*.yaml)")
    parser.add_argument("base_dir", type=str, help="Base folder containing the simulation results")
    args = parser.parse_args()
    analyze_yaml_plots(args.base_dir)
