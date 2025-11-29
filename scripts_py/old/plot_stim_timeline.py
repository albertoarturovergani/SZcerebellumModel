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

    plt.figure(figsize=(10, 5))
    for i, stim in enumerate(stimuli):
        plt.barh(i, stim["stop"] - stim["start"], left=stim["start"],
                 height=0.8, color="firebrick" if stim["type"] == "burst" else "steelblue")
        plt.text(stim["start"], i, f'{stim["rate"]} Hz', va='center', ha='right', fontsize=8, color='white')

    plt.yticks(range(len(stimuli)), [s["label"] for s in stimuli])
    plt.xlabel("Time (ms)")
    plt.title("Stimulus Timeline")
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
