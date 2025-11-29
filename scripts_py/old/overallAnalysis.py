import os
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import argparse

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

    data = []
    for yfile in yaml_files:
        try:
            with open(yfile, 'r') as f:
                ydata = yaml.safe_load(f)
            basename = os.path.basename(yfile)
            k_val = ydata.get('k')
            if k_val is None:
                k_val = float(basename.replace('summary_k_', '').replace('.yaml', ''))

            flat = flatten_dict(ydata)
            flat['k'] = k_val
            data.append(flat)
        except Exception as e:
            print(f"Errore nel file {yfile}: {e}")

    df = pd.DataFrame(data)
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna(subset=['k']).sort_values(by='k')

    if df.empty:
        print("⚠️  Nessun dato valido trovato nei file YAML.")
        return

    out_dir = os.path.join(base_dir, "overallAnalysis")
    os.makedirs(out_dir, exist_ok=True)

    for var in df.columns:
        if var == 'atrophy_percent':
            continue
        if df[var].dropna().empty:
            continue

        plt.figure(figsize=(8, 5))
        plt.plot(df['atrophy_percent'], df[var], marker='o')
        plt.xlabel('atrophy_percent')
        plt.ylabel(var)
        plt.title(f'{var} vs k')
        plt.grid(True)
        plt.tight_layout()
        filepath = os.path.join(out_dir, f"{var}_vs_k.png")
        plt.savefig(filepath)
        plt.close()
        print(f"✅ Salvato: {filepath}")

    csv_path = os.path.join(out_dir, "summary_data.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n📄 CSV salvato in: {csv_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analizza e plottizza i file summary_k_*.yaml")
    parser.add_argument("base_dir", type=str, help="Percorso base della cartella results_*/")
    args = parser.parse_args()
    analyze_yaml_plots(args.base_dir)
