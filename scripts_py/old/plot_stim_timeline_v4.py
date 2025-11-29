#!/usr/bin/env python3
import os
import sys
import yaml
import matplotlib.pyplot as plt

def generate_stimulus_plots(sim_dir: str):
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
    stim_dir = os.path.join(sim_dir, "stimulus")
    os.makedirs(stim_dir, exist_ok=True)

    # === Plot 1: Stimulus Timeline ===
    plt.figure(figsize=(10, 6))
    for i, stim in enumerate(stimuli):
        color = "firebrick" if stim["type"] == "burst" else "steelblue"
        plt.barh(i, stim["stop"] - stim["start"], left=stim["start"], height=0.8, color=color)
        plt.text(stim["start"], i, f'{stim["rate"]} Hz', va='center', ha='right', fontsize=8, color='white')
    plt.yticks(range(len(stimuli)), [s["label"] for s in stimuli])
    plt.xlabel("Time (ms)")
    plt.title("Stimulus Timeline")
    plt.tight_layout()
    plt.savefig(os.path.join(stim_dir, "stimulus_timeline.png"))
    plt.close()
    print(f"✅ Salvato: stimulus_timeline.png")

    # === Plot 2: Rate vs Time ===
    midpoints = [(s["start"] + s["stop"]) / 2 for s in stimuli]
    rates = [s["rate"] for s in stimuli]
    plt.figure(figsize=(8, 5))
    plt.scatter(midpoints, rates, marker="o", color="darkgreen")
    plt.xlabel("Time (ms)")
    
    plt.ylabel("Rate (Hz)")
    plt.title("Rate vs Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(stim_dir, "stimulus_rate_vs_time.png"))
    plt.close()
    print(f"✅ Salvato: stimulus_rate_vs_time.png")

    # === Plot 3: Durata stimolo vs Rate ===
    durations = [s["stop"] - s["start"] for s in stimuli]
    plt.figure(figsize=(8, 5))
    plt.scatter(durations, rates, marker="o")
    plt.xlabel("Durata stimolo (ms)")
    plt.ylabel("Rate (Hz)")
    plt.title("Durata stimolo vs Rate")
    plt.xlim(-10, 200)

    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(stim_dir, "stimulus_duration_vs_rate.png"))
    plt.close()
    print(f"✅ Salvato: stimulus_duration_vs_rate.png")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_stim_timeline_v2.py <sim_dir>")
        sys.exit(1)

    sim_dir = sys.argv[1]
    generate_stimulus_plots(sim_dir)
