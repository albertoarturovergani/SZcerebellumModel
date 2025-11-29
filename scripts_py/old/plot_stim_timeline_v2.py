import os
import sys
import yaml
import matplotlib.pyplot as plt

def generate_stimulus_timeline(sim_dir: str, output_filename: str = "stimulus_timeline.png"):
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

    # 📊 Plot con due pannelli affiancati
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [2, 1]})

    # Pannello sinistro: timeline orizzontale
    for i, stim in enumerate(stimuli):
        color = "firebrick" if stim["type"] == "burst" else "steelblue"
        ax1.barh(i, stim["stop"] - stim["start"], left=stim["start"], height=0.8, color=color)
        ax1.text(stim["start"], i, f'{stim["rate"]} Hz', va='center', ha='right', fontsize=8, color='white')
    ax1.set_yticks(range(len(stimuli)))
    ax1.set_yticklabels([s["label"] for s in stimuli])
    ax1.set_xlabel("Time (ms)")
    ax1.set_title("Stimulus Timeline")

    # Pannello destro: rate vs time
    midpoints = [(s["start"] + s["stop"]) / 2 for s in stimuli]
    rates = [s["rate"] for s in stimuli]
    ax2.scatter(midpoints, rates, marker="o", linestyle="-", color="darkgreen")
    ax2.set_title("Rate vs Time")
    ax2.set_xlabel("Time (ms)")
    ax2.set_ylabel("Rate (Hz)")
    ax2.grid(True)

    plt.tight_layout()
    out_path = os.path.join(sim_dir, output_filename)
    plt.savefig(out_path)
    print(f"✅ Stimulus timeline salvata in: {out_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_stim_timeline.py <sim_dir>")
        sys.exit(1)

    sim_dir = sys.argv[1]
    generate_stimulus_timeline(sim_dir)
