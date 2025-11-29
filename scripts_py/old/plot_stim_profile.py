import yaml
import matplotlib.pyplot as plt
import argparse
import os

def plot_stimulus_yaml_colored(yaml_file, save_path=None):
    with open(yaml_file, 'r') as file:
        data = yaml.safe_load(file)

    sim_key = list(data['simulations'].keys())[0]
    devices = data['simulations'][sim_key]['devices']

    times = []
    rates = []
    labels = []
    colors = []

    for name, device in devices.items():
        if device['device'] == 'poisson_generator':
            # Ricostruzione robusta di start e stop
            start = device.get('start', 0.0)  # Se manca start, usa 0.0
            if 'stop' in device:
                stop = device['stop']
            elif 'duration' in device:
                stop = start + device['duration']
            else:
                print(f"⚠️  Device {name} has no stop or duration, skipping...")
                continue

            rate = device['rate']
            times.append((start, stop))
            rates.append(rate)
            labels.append(name)

            if 'baseline' in name:
                colors.append('blue')
            elif 'burst' in name:
                colors.append('red')
            else:
                colors.append('gray')

    fig, ax = plt.subplots(figsize=(12, 6))

    for (start, stop), rate, label, color in zip(times, rates, labels, colors):
        ax.hlines(rate, start, stop, colors=color, linewidth=2)
        if color == 'red':
            ax.text((start + stop) / 2, rate + 10, label, ha='center', fontsize=8, color=color)
    
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Firing rate (Hz)')
    ax.set_title('Stimulus profile: Baseline (blue) and Bursts (red)')
    ax.grid(True)
    plt.tight_layout()

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        base_name = os.path.splitext(os.path.basename(yaml_file))[0]
        save_file = os.path.join(save_path, f"{base_name}.png")
        fig.savefig(save_file)
        print(f"✅ Plot salvato in: {save_file}")
    else:
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot stimulus YAML profile.")
    parser.add_argument('yaml_file', type=str, help='Path to the stimulus YAML file')
    parser.add_argument('--savePath', type=str, default=None, help='Path to save the plot (optional)')
    args = parser.parse_args()

    plot_stimulus_yaml_colored(args.yaml_file, args.savePath)
