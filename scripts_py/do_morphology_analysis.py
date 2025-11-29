#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import spearmanr

def collect_morphology_csvs(base_path: Path):
    output_path = base_path / "B_collective_results" / "1.structures"
    output_path.mkdir(parents=True, exist_ok=True)

    all_rows = []

    for csv_path in base_path.rglob("*/1.structure/morphologies/*.csv"):
        parts = csv_path.parts
        k_match = [p for p in parts if "_deg_" in p]
        if not k_match:
            print(f"⚠️ Skip (no k found): {csv_path}")
            continue
        try:
            k_str = k_match[0].split("_")[2]
            k = float(k_str)
        except Exception as e:
            print(f"⚠️ Errore parsing k da '{k_match[0]}': {e}")
            continue

        try:
            df = pd.read_csv(csv_path)
            df["k"] = k
            df["Atrophy Level (%)"] = 100 * (1 - k)
            df["source_file"] = str(csv_path.relative_to(base_path))
            all_rows.append(df)
            print(f"✅ {csv_path.name} → k={k}")
        except Exception as e:
            print(f"❌ Errore nel file {csv_path}: {e}")

    if not all_rows:
        print("⚠️ Nessun file valido trovato.")
        return

    df_all = pd.concat(all_rows, ignore_index=True)
    output_file = output_path / "all_morphologies.pkl"
    df_all.to_pickle(output_file)
    print(f"\n💾 File unificato salvato in formato PKL: {output_file}")

    plot_morpho_vs_atrophy(output_file)

def plot_morpho_vs_atrophy(pkl_path: Path):
    df = pd.read_pickle(pkl_path)
    output_dir = pkl_path.parent
    base_name = pkl_path.stem

    exclude_cols = ['k', 'Atrophy Level (%)', 'source_file', 'cell_name', 'num_branch',
'n_axon', 'n_soma', 'n_dend', 'k_shrinking',
'k_radius', 'k_pruning', 'k_gen', 'Atrophy Factor (%)',
'num_sections', 'num_tip_orders', 'num_branch_tips',
                   ]
    morpho_cols = [col for col in df.columns if col not in exclude_cols]

    if "cell_name" not in df.columns:
        print("⚠️ La colonna 'cell_name' (tipo cellulare) è mancante.")
        return

    cell_types = sorted(df["cell_name"].unique())

    for morpho_col in morpho_cols:
        fig, axes = plt.subplots(nrows=1, ncols=len(cell_types), figsize=(5 * len(cell_types), 4), sharey=False)
        if len(cell_types) == 1:
            axes = [axes]

        for i, cell in enumerate(cell_types):
            ax = axes[i]
            sub = df[df["cell_name"] == cell]
            x = sub["Atrophy Level (%)"]
            y = sub[morpho_col]
            sns.scatterplot(x=x, y=y, ax=ax, color='black', alpha=0.6)
            #sns.regplot(x=x, y=y, ax=ax, scatter=False, color='red', ci=None)

            if len(x) >= 2 and not np.all(np.isnan(y)):
                rho, pval = spearmanr(x, y, nan_policy='omit')
                ax.set_title(f"{cell}\nρ = {rho:.2f}, p = {pval:.1e}")
            else:
                ax.set_title(cell)

            ax.set_xlabel("Atrophy Level (%)")
            ax.set_ylabel(morpho_col)

        fig.suptitle(f"{morpho_col} vs Atrophy Level", fontsize=16, fontweight="bold", y=1.1)
        fig.tight_layout()
        plt.subplots_adjust(top=0.85)

        save_path = output_dir / f"{base_name}_atrophy_vs_{morpho_col}.png"
        fig.savefig(save_path, dpi=300)
        print(f"💾 Figura salvata: {save_path}")
        plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatena tutti i CSV di morfologie da sottocartelle /1.structure/morphologies/ e salva in .pkl")
    parser.add_argument("--path", required=True, help="Percorso alla directory base")

    args = parser.parse_args()
    base_path = Path(args.path).expanduser().resolve()

    if not base_path.exists():
        print(f"❌ Il percorso non esiste: {base_path}")
    else:
        collect_morphology_csvs(base_path)
