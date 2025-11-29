import os
import re
import argparse
import yaml
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def mirrored_sigmoid(k_array, midpoint, sharpness=20):
    raw = 1 / (1 + np.exp(-sharpness * (k_array - midpoint)))
    return 0.5 + (raw - raw.min()) / (raw.max() - raw.min()) * 0.5

def collect_structure_summaries(base_path):
    output_dir = os.path.join(base_path, "B_collective_results", "1.structures")
    os.makedirs(output_dir, exist_ok=True)
    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    info_file = os.path.join(base_path, "info.yaml")
    if not os.path.exists(info_file):
        print("❌ File info.yaml non trovato.")
        return

    with open(info_file, "r") as f:
        info_data = yaml.safe_load(f)

    reserve_k0 = info_data.get("modeling_parameters", {}).get("reserve_k0", 0.75)
    reserve_sharpness = info_data.get("modeling_parameters", {}).get("reserve_sharpness", 20)

    pattern = re.compile(r"mouse_cerebellar_cortex_(\d+\.\d+)_structure_summary\.yaml")

    summary_rows = []
    detailed_placement = []
    detailed_connectivity = []

    for root, dirs, files in os.walk(base_path):
        for file in files:
            match = pattern.match(file)
            if match:
                k_str = match.group(1)
                k = float(k_str)

                file_path = os.path.join(root, file)
                with open(file_path, "r") as f:
                    data = yaml.safe_load(f)

                summary_rows.append({
                    "k": k,
                    "atrophy_pct": (1 - k) * 100,
                    "total_cells": data["placement"]["total_cells"],
                    "mean_density": data["placement"]["mean_density"],
                    "max_density": data["placement"]["max_density"],
                    "total_synapses": data["connectivity"]["total_synapses"],
                    "mean_syn_per_pair": data["connectivity"]["mean_syn_per_pair"],
                    "mean_convergence": data["connectivity"]["mean_convergence"],
                    "mean_divergence": data["connectivity"]["mean_divergence"],
                })

                for cell_type, vals in data["placement"]["cell_types"].items():
                    detailed_placement.append({
                        "k": k,
                        "atrophy_pct": (1 - k) * 100,
                        "cell_type": cell_type,
                        "count": vals["count"],
                        "density": vals["density"]
                    })

                for conn, vals in data["connectivity"]["connections"].items():
                    detailed_connectivity.append({
                        "k": k,
                        "atrophy_pct": (1 - k) * 100,
                        "connection": conn,
                        "synapses": vals["synapses"],
                        "syn_per_pair": vals["syn_per_pair"],
                        "convergence": vals["convergence"],
                        "divergence": vals["divergence"]
                    })

    if not summary_rows:
        print("❌ Nessun file YAML trovato.")
        return

    df_summary = pd.DataFrame(summary_rows).sort_values("k").reset_index(drop=True)
    k_array = df_summary["k"].values
    df_summary["cell_strength"] = mirrored_sigmoid(k_array, midpoint=reserve_k0, sharpness=reserve_sharpness)

    df_cells = pd.DataFrame(detailed_placement)
    df_conns = pd.DataFrame(detailed_connectivity)
    df_cells = df_cells.merge(df_summary[["k", "cell_strength"]], on="k")
    df_conns = df_conns.merge(df_summary[["k", "cell_strength"]], on="k")

    df_summary.to_pickle(os.path.join(output_dir, "summary_structures.pkl"))
    df_cells.to_pickle(os.path.join(output_dir, "summary_cells.pkl"))
    df_conns.to_pickle(os.path.join(output_dir, "summary_connections.pkl"))

    # === PLOT GLOBALI MULTIPANNELLO ===
    def plot_metrics(xvar, label_suffix):
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        metrics = [
            ("total_cells", "Total cells"),
            ("mean_density", "Mean density [μm⁻³]"),
            ("total_synapses", "Total synapses"),
            ("mean_syn_per_pair", "Synapses per pair"),
            ("mean_convergence", "Mean convergence"),
            ("mean_divergence", "Mean divergence")
        ]
        for ax, (col, title) in zip(axes.flat, metrics):
            sns.lineplot(data=df_summary, x=xvar, y=col, ax=ax, marker="o")
            ax.set_title(title)
            ax.set_xlabel(label_suffix)
            ax.grid(True)
        plt.tight_layout()
        fig.savefig(os.path.join(plot_dir, f"summary_structures_all_panels_by_{xvar}.png"), dpi=300)
        plt.close()

        for metric in ["count", "density"]:
            g = sns.FacetGrid(df_cells, col="cell_type", col_wrap=4, sharey=False, sharex=False, height=4)
            g.map_dataframe(sns.lineplot, x=xvar, y=metric, marker="o")
            g.set_titles("{col_name}")
            g.set_axis_labels(label_suffix, metric.capitalize())
            for ax in g.axes.flat:
                ax.grid(True)
            g.fig.suptitle(f"Cell types – {metric} by {label_suffix}", fontsize=16)
            g.fig.tight_layout()
            g.fig.subplots_adjust(top=0.9)
            g.savefig(os.path.join(plot_dir, f"cells_{metric}_by_type_by_{xvar}.png"))
            plt.close()

        for metric in ["synapses", "syn_per_pair", "convergence", "divergence"]:
            g = sns.FacetGrid(df_conns, col="connection", col_wrap=4, sharey=False, sharex=False, height=4)
            g.map_dataframe(sns.lineplot, x=xvar, y=metric, marker="o")
            g.set_titles("{col_name}")
            g.set_axis_labels(label_suffix, metric.replace("_", " ").capitalize())
            for ax in g.axes.flat:
                ax.grid(True)
            g.fig.suptitle(f"Connectivity – {metric} by {label_suffix}", fontsize=16)
            g.fig.tight_layout()
            g.fig.subplots_adjust(top=0.9)
            g.savefig(os.path.join(plot_dir, f"connections_{metric}_by_type_by_{xvar}.png"))
            plt.close()

    plot_metrics("atrophy_pct", "Atrophy level (%)")
    plot_metrics("cell_strength", "Cell strength")

    # === PLOT RELAZIONE cell_strength vs atrophy_pct + identità ===
    plt.figure(figsize=(6, 5))
    sns.lineplot(data=df_summary, x="atrophy_pct", y="cell_strength", marker="o", label="Cell strength")
    # Linea identità normalizzata
    x = df_summary["atrophy_pct"]
    x_norm = 1 - (x - x.min()) / (x.max() - x.min()) * 0.5
    plt.plot(x, x_norm, '--', color='gray', label="Identity")
    plt.title(f"Cell Strength vs Atrophy Level \n sigmoid(k; midpoint={reserve_k0}, sharpness={reserve_sharpness})")
    plt.xlabel("Atrophy level (%)")
    plt.ylabel("Cell strength (normalized)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "cell_strength_vs_atrophy.png"), dpi=300)
    plt.close()
    print("📈 Salvato: cell_strength_vs_atrophy.png")


    
    return df_summary, df_cells, df_conns




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Colleziona YAML strutturali e salva tabelle + grafici")
    parser.add_argument("--path", required=True, help="Percorso alla cartella esperimento")
    args = parser.parse_args()

    collect_structure_summaries(args.path)
