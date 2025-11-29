#!/usr/bin/env python3
import os
import sys
import yaml
import matplotlib.pyplot as plt

def generate_stimulus_plots(yaml_file_path: str):
    if not os.path.isfile(yaml_file_path):
        print(f"❌ YAML file non trovato: {yaml_file_path}")
        return

    with open(yaml_file_path, 'r') as f:
        ydata = yaml.safe_load(f)

    stimuli = []
    try:
        devices = ydata['simulations']['basal_activity']['devices']
        for label, stim in devices.items():
            if stim.get("device") == "poisson_generator":
                stimuli.append({
                    #"label": label,
                    #"start": stim["start"],
                    #"stop": stim["stop"],
                    "rate": stim["rate"],
                    #"type": "burst" if "burst" in label else "stim"
                })
    except Exception:
        print(f"⚠️ Nessun stimolo valido trovato in {yaml_file_path}")
        return

    if not stimuli:
        print(f"⚠️ Nessun stimolo valido trovato in {yaml_file_path}")
        return

    stimuli.sort(key=lambda x: x["start"])

    output_dir = os.path.join(os.path.dirname(yaml_file_path), "stimulus")
    os.makedirs(output_dir, exist_ok=True)

    # Estrai prefisso completo del file YAML senza estensione
    base_name = os.path.basename(yaml_file_path).replace(".yaml", "")

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
    save_path = os.path.join(output_dir, f"{base_name}_stimulus_timeline.png")
    plt.savefig(save_path)
    plt.close()
    print(f"✅ Salvato: {save_path}")

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
    save_path = os.path.join(output_dir, f"{base_name}_stimulus_rate_vs_time.png")
    plt.savefig(save_path)
    plt.close()
    print(f"✅ Salvato: {save_path}")

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
    save_path = os.path.join(output_dir, f"{base_name}_stimulus_duration_vs_rate.png")
    plt.savefig(save_path)
    plt.close()
    print(f"✅ Salvato: {save_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_stim_timeline_v2.py <stimulus_yaml_file>")
        sys.exit(1)

    yaml_file_path = sys.argv[1]
    generate_stimulus_plots(yaml_file_path)
