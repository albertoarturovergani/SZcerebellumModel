import os
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import argparse
import seaborn as sns

import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize

def plot_overlapped_psd(df_all, out_dir):
    for pop in ['mossy', 'purkinje']:
        plt.figure(figsize=(10, 6))

        # Estrai tutti i valori k validi per la colormap
        valid_k = df_all['k'].dropna().values
        norm = Normalize(vmin=np.min(valid_k), vmax=np.max(valid_k))
        cmap = get_cmap("viridis")

        for _, row in df_all.iterrows():
            k = row.get('k', None)
            freqs = row.get(f'psd__{pop}__frequencies_hz', None)
            power = row.get(f'psd__{pop}__power', None)

            if isinstance(freqs, list) and isinstance(power, list) and k is not None:
                color = cmap(norm(k))
                plt.plot(freqs, power, label=f'k={k:.2f}', color=color, alpha=0.9)

        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power")
        plt.xlim(0, 100)
        plt.title(f"Overlapped PSD – {pop.capitalize()} population")
        plt.grid(True)

        # Legenda ordinata per k crescente
        handles, labels = plt.gca().get_legend_handles_labels()
        if handles:
            sorted_labels_handles = sorted(zip(labels, handles), key=lambda x: float(x[0][2:]))  # "k=0.80" → 0.80
            labels, handles = zip(*sorted_labels_handles)
            plt.legend(handles, labels, title='k', bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.tight_layout()
        save_path = os.path.join(out_dir, f"overlapped_psd_{pop}.png")
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"📈 Salvato plot PSD sovrapposte: {save_path}")



def flatten_dict(d, parent_key='', sep='__'):
    """Flatten dizionari annidati"""
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def analyze_yaml_plots(base_dir):
    pattern = os.path.join(base_dir, "results_*/basal_activity/plots_firing_rate_*/summary_k_*.yaml")
    yaml_files = glob(pattern)

    print(f"\n📁 YAML trovati: {len(yaml_files)}")
    if len(yaml_files) == 0:
        print("⚠️  Nessun file trovato. Controlla il path o la struttura delle cartelle.")
        return

    all_data = []

    for yfile in yaml_files:
        try:
            with open(yfile, 'r') as f:
                ydata = yaml.safe_load(f)
            flat = flatten_dict(ydata)
            # Ricava k e finestra temporale dal nome del file
            basename = os.path.basename(yfile).replace('.yaml', '')
            parts = basename.split('_')
            if len(parts) >= 4:
                flat['time_window'] = parts[-1]  # es. 'T1'
                flat['k'] = float(parts[2])
            else:
                flat['time_window'] = 'Overall Time'
                flat['k'] = float(parts[2])
            all_data.append(flat)
        except Exception as e:
            print(f"Errore nel file {yfile}: {e}")

    df_all = pd.DataFrame(all_data)
    df_all = df_all.apply(pd.to_numeric, errors='ignore')

    if df_all.empty or 'time_window' not in df_all.columns:
        print("⚠️  Nessun dato valido trovato nei file YAML.")
        return

    out_dir = os.path.join(base_dir, "overallAnalysis")
    os.makedirs(out_dir, exist_ok=True)

    for window in sorted(df_all['time_window'].unique()):
        df_window = df_all[df_all['time_window'] == window]
        df_window = df_window.dropna(subset=['atrophy_percent'])

        print(f"\n📊 Analisi per {window} – {len(df_window)} righe")

        out_dir_window = os.path.join(out_dir, window)
        os.makedirs(out_dir_window, exist_ok=True)

        for var in df_window.columns:
            if var in ['k', 'atrophy_percent', 'time_window']:
                continue
        
            # Filtra righe valide (non NaN o non lista)
            valid_values = df_window[var].dropna()
            
            # Salta se contiene liste o dizionari
            if valid_values.apply(lambda x: isinstance(x, (list, dict))).any():
                continue
        
            if valid_values.empty:
                continue
        
            plt.figure(figsize=(8, 5))
            sns.regplot(x='atrophy_percent', y=var, data=df_window, order=2)
            plt.xlabel('Atrophy (%)')
            plt.ylabel(var)
            plt.title(f'{var} vs atrophy – {window}')
            plt.grid(False)
            plt.tight_layout()
            filepath = os.path.join(out_dir_window, f"{var}_vs_atrophy.png")
            plt.savefig(filepath)
            plt.close()
            print(f"✅ Salvato: {filepath}")


        # Salva CSV per la finestra
        df_window.to_csv(os.path.join(out_dir_window, f"summary_data_{window}.csv"), index=False)
        
        if window == 'Overall Time':
            plot_overlapped_psd(df_window, out_dir_window)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analizza e plottizza i file summary_k_*.yaml per T1-T5")
    parser.add_argument("base_dir", type=str, help="Cartella base che contiene results_*/")
    args = parser.parse_args()
    analyze_yaml_plots(args.base_dir)
