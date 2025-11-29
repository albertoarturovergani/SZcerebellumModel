#!/usr/bin/env python3

def main():
    import gc
    import neo
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from pathlib import Path
    from collections import defaultdict
    from scipy.signal import convolve
    from scipy.signal.windows import triang
    from pylatex import Document, Section, Figure, NoEscape
    import argparse
    import re

    def compute_firing_rate(spike_times, duration, dt=1.0, sigma=100):
        time = np.arange(0, duration, dt)
        spike_train = np.zeros_like(time)
        indices = (np.array(spike_times) / dt).astype(int)
        indices = indices[indices < len(spike_train)]
        spike_train[indices] = 1

        kernel_size = int(2 * sigma / dt)
        kernel_size = max(kernel_size, 3)
        kernel = triang(kernel_size)
        kernel /= kernel.sum()

        rate = convolve(spike_train, kernel, mode='same') * (1000.0 / dt)
        return time, rate

    filter_mode = True
    output_dir = Path("plots_firing_rate")
    output_dir.mkdir(exist_ok=True)
    image_paths = []
    trim_ratio = 0.005
    dt = 0.1
    sigma = 100

    parser = argparse.ArgumentParser(description="Process neuronal firing rates.")
    parser.add_argument("root", nargs="?", default=".", help="Directory to search for results_*/basal_activity")
    args = parser.parse_args()
    root_path = Path(args.root)

    for result_dir in root_path.rglob('results_*/basal_activity'):
        # Trova il pezzo di path che corrisponde a "results_FLOAT"
        match = [part for part in result_dir.parts if re.match(r"results_\d+\.\d+", part)]
        if not match:
            print(f"⚠️ Impossibile determinare k da: {result_dir}")
            continue
        kvalue = float(match[0].split('_')[1])

        if kvalue > 0:
            output_dir = result_dir / f"plots_firing_rate_{kvalue:.3f}"
            output_dir.mkdir(parents=True, exist_ok=True)
            atrophy = 100 * (1 - kvalue)
            spike_rows = []
            for file in result_dir.glob("*.nio"):
                print(f"Reading: {file.name} from {result_dir}")
                reader = neo.io.NixIO(filename=str(file), mode="ro")
                block = reader.read_block()

                for seg_idx, segment in enumerate(block.segments):
                    for st_idx, spiketrain in enumerate(segment.spiketrains):
                        spike_times = spiketrain.times.rescale("ms").magnitude
                        senders = spiketrain.annotations.get("senders", None)
                        device = spiketrain.annotations.get("device", "unknown")
                        if senders is not None:
                            senders = np.atleast_1d(senders)
                            spike_times = np.atleast_1d(spike_times)
                            if len(senders) == len(spike_times):
                                for time, sender in zip(spike_times, senders):
                                    spike_rows.append({
                                        "time_ms": time,
                                        "sender_id": sender,
                                        "pop": device.split('_')[0],
                                        "segment": seg_idx,
                                        "file": file.name
                                    })
                            else:
                                print(f"⚠️ Mismatch senders/spike_times in {file.name}, segment {seg_idx}")
                        else:
                            for time in spike_times:
                                spike_rows.append({
                                    "time_ms": time,
                                    "sender_id": f"unit{st_idx}",
                                    "pop": device.split('_')[0],
                                    "segment": seg_idx,
                                    "file": file.name
                                })


                        """
                        
                        if senders and len(senders) == len(spike_times):
                            for time, sender in zip(spike_times, senders):
                                spike_rows.append({
                                    "time_ms": time,
                                    "sender_id": sender,
                                    "pop": device.split('_')[0],
                                    "segment": seg_idx,
                                    "file": file.name
                                })
                        else:
                            for time in spike_times:
                                spike_rows.append({
                                    "time_ms": time,
                                    "sender_id": f"unit{st_idx}",
                                    "pop": device.split('_')[0],
                                    "segment": seg_idx,
                                    "file": file.name
                                })
                        """
                del reader, block
                gc.collect()

            df_spikes = pd.DataFrame(spike_rows)
            if df_spikes.empty:
                continue

            if filter_mode:
                df_spikes = df_spikes[df_spikes['pop'].isin(['mossy', 'purkinje'])]
            else:
                df_spikes = df_spikes[~df_spikes['pop'].isin(['stim', 'baseline'])]
                granule_name = 'granule'
                granule_ids = df_spikes[df_spikes['pop'] == granule_name]['sender_id'].unique()
                if len(granule_ids) > 0:
                    subset_ids = np.random.choice(granule_ids, size=int(0.001 * len(granule_ids)), replace=False)
                    df_spikes = df_spikes[(df_spikes['pop'] != granule_name) | (df_spikes['sender_id'].isin(subset_ids))]

            duration_ms = df_spikes["time_ms"].max()
            grouped = df_spikes.groupby(["pop", "sender_id"])

            firing_data = defaultdict(list)
            for (pop, sender_id), group in grouped:
                spike_times = group["time_ms"].values
                time, rate = compute_firing_rate(spike_times, duration=duration_ms, dt=dt, sigma=sigma)
                firing_data["pop"].append(pop)
                firing_data["sender_id"].append(sender_id)
                firing_data["time"].append(time)
                firing_data["rate"].append(rate)

            mossy_rates = []
            purkinje_rates = []
            time_vector = firing_data["time"][0] if firing_data["time"] else None

            for i in range(len(firing_data["pop"])):
                pop = firing_data["pop"][i]
                rate = firing_data["rate"][i]
                if pop == "mossy":
                    mossy_rates.append(rate)
                elif pop == "purkinje":
                    purkinje_rates.append(rate)

            if time_vector is None or not mossy_rates or not purkinje_rates:
                print(f"⚠️ Skipping k={kvalue}: insufficient data.")
                continue

            n = len(time_vector)
            start_idx = int(n * trim_ratio)
            end_idx = int(n * (1 - trim_ratio))

            t = time_vector[start_idx:end_idx]
            mossy_array = np.vstack(mossy_rates)[:, start_idx:end_idx]
            purkinje_array = np.vstack(purkinje_rates)[:, start_idx:end_idx]

            mossy_mean = np.mean(mossy_array, axis=0)
            mossy_std = np.std(mossy_array, axis=0)
            purkinje_mean = np.mean(purkinje_array, axis=0)
            purkinje_std = np.std(purkinje_array, axis=0)

            similarity = np.dot(mossy_mean, purkinje_mean) / (np.linalg.norm(mossy_mean) * np.linalg.norm(purkinje_mean))

            fig = plt.figure(figsize=(10, 4.5))
            ax_main = fig.add_subplot(111)
            ax_main.plot(t, mossy_mean, label="Mossy fibers (input)", linewidth=1.5)
            ax_main.fill_between(t, mossy_mean - mossy_std, mossy_mean + mossy_std, alpha=0.3)
            ax_main.plot(t, purkinje_mean, label="Purkinje cells (output)", linewidth=1.5)
            ax_main.fill_between(t, purkinje_mean - purkinje_std, purkinje_mean + purkinje_std, alpha=0.3)
            ax_main.set_ylim(0, 400)
            ax_main.set_xlabel("Time (ms)")
            ax_main.set_ylabel("Firing rate (Hz)")
            ax_main.set_title(f"Mean firing rate ± SD – {atrophy:.0f}% atrophy\nCosine similarity = {similarity:.3f}")
            ax_main.grid(True)
            ax_main.legend(loc='upper left')

            filename_base = output_dir / f"fr_plot_k_{kvalue:.3f}"
            fig.savefig(f"{filename_base}.png", dpi=300)
            fig.savefig(f"{filename_base}.svg")
            if filename_base.with_suffix(".png").exists():
                image_paths.append(filename_base.with_suffix(".png"))
            plt.close(fig)

            print(f"✅ Saved: {filename_base}.png / .svg")

if __name__ == "__main__":
    main()
