import argparse
import yaml
import matplotlib.pyplot as plt
from pathlib import Path
import glob

def generate_stimulus_timeline(yaml_path, stimName='basal_activity', outPath=None):
    yaml_path = Path(yaml_path)
    if outPath is None:
        outPath = yaml_path.parent / f"{yaml_path.stem}_timeline.png"
    else:
        outPath = Path(outPath) / f"{yaml_path.stem}_timeline.png"

    try:
        with open(yaml_path, 'r') as f:
            ydata = yaml.safe_load(f)
        devices = ydata['simulations'][stimName]['devices']
    except Exception as e:
        print(f"⚠️ Errore parsing {yaml_path}: {e}")
        return False

    stimuli = []
    for label, stim in devices.items():
        if stim.get("device") == "poisson_generator":
            stimuli.append({
                "label": label,
                "start": float(stim["start"]),
                "stop": float(stim["stop"]),
                "rate": float(stim["rate"]),
                "type": "burst" if "burst" in label else "stim"
            })

    if not stimuli:
        print(f"⚠️ Nessuno stimolo trovato in {yaml_path}")
        return False

    stimuli.sort(key=lambda x: x["start"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [2, 1]})
    for i, stim in enumerate(stimuli):
        color = "firebrick" if stim["type"] == "burst" else "steelblue"
        ax1.barh(i, stim["stop"] - stim["start"], left=stim["start"], height=0.8, color=color)
        ax1.text(stim["start"], i, f'{stim["rate"]} Hz', va='center', ha='right', fontsize=8, color='white')

    ax1.set_yticks(range(len(stimuli)))
    ax1.set_yticklabels([s["label"] for s in stimuli])
    ax1.set_xlabel("Time (ms)")
    ax1.set_title("Stimulus Timeline")

    midpoints = [(s["start"] + s["stop"]) / 2 for s in stimuli]
    rates = [s["rate"] for s in stimuli]
    ax2.scatter(midpoints, rates, marker="o", color="darkgreen")
    ax2.set_title("Rate vs Time")
    ax2.set_xlabel("Time (ms)")
    ax2.set_ylabel("Rate (Hz)")
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig(outPath, dpi=300)
    plt.close()
    print(f"✅ Salvato: {outPath}")
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genera timeline stimoli da YAML.")
    parser.add_argument('--stimDir', help="Glob path dei file YAML, es. './folder/*.yaml'")
    parser.add_argument('--stimName', default='basal_activity', help="Nome stimolo nel file YAML")
    parser.add_argument('--outPath', default=None, help="Cartella output (default: stessa del file YAML)")

    args = parser.parse_args()

    # supporto per stimDir glob
    yaml_files = glob.glob(args.stimDir)
    if not yaml_files:
        print(f"❌ Nessun file trovato per {args.stimDir}")
    for yaml_path in yaml_files:
        generate_stimulus_timeline(yaml_path, stimName=args.stimName, outPath=args.outPath)
