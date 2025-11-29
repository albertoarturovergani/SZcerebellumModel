import os
import glob
import yaml
import matplotlib.pyplot as plt

def generate_stimulus_timeline(yaml_path):
    with open(yaml_path, 'r') as f:
        ydata = yaml.safe_load(f)

    stimuli = []
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
        return False

    if not stimuli:
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
    ax2.scatter(midpoints, rates, marker="o", linestyle="-", color="darkgreen")
    ax2.set_title("Rate vs Time")
    ax2.set_xlabel("Time (ms)")
    ax2.set_ylabel("Rate (Hz)")
    ax2.grid(True)

    plt.tight_layout()
    out_path = os.path.join(os.path.dirname(yaml_path), "stimulus_timeline.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"✅ Salvato: {out_path}")
    return True


def process_reconstruction_yamls(main_dir):
    pattern = os.path.join(main_dir, "results_*/reconstruction/*.yaml")
    yaml_paths = glob.glob(pattern)
    if not yaml_paths:
        print(f"⚠️ Nessun file YAML trovato con pattern: {pattern}")
        return

    count = 0
    for yaml_path in yaml_paths:
        if generate_stimulus_timeline(yaml_path):
            count += 1

    print(f"\n🎯 Generati {count} timeline plot da {len(yaml_paths)} YAML trovati.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python plot_stim_timelines_recursive.py <main_dir>")
        sys.exit(1)

    main_dir = sys.argv[1]
    process_reconstruction_yamls(main_dir)
