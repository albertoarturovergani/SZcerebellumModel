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
    import yaml
    import numpy as np
    from collections import Counter
    from math import log2
    from fooof import FOOOF
    from scipy.signal import welch
    from scipy.signal import hilbert
    from scipy.stats import pearsonr

    def compute_population_psd(array, fs):
        if array.size == 0:
            return [], []
        psds = []
        for unit in array:
            f, p = welch(unit, fs=fs, nperseg=1024)
            psds.append(p)
        psd_mean = np.mean(psds, axis=0)
        return f.tolist(), psd_mean.tolist()


    
    # Phase Locking Value (PLV) – via Hilbert transform
    def compute_plv(sig1, sig2):
        phase1 = np.angle(hilbert(sig1))
        phase2 = np.angle(hilbert(sig2))
        phase_diff = phase1 - phase2
        plv = np.abs(np.mean(np.exp(1j * phase_diff)))
        return float(plv)
    
    # Correlazione media tra tutte le unità di un gruppo
    def mean_pairwise_corr(array):
        n_units = array.shape[0]
        if n_units < 2:
            return np.nan
        corrs = []
        for i in range(n_units):
            for j in range(i+1, n_units):
                r, _ = pearsonr(array[i], array[j])
                corrs.append(r)
        return float(np.nanmean(corrs)) if corrs else np.nan

    def analyze_and_plot_fooof(signal, fs, title, out_path, freq_range=(1, 100)):
        from fooof.plts.spectra import plot_spectrum
        freqs, psd = welch(signal, fs=fs, nperseg=1024)
        fm = FOOOF(peak_width_limits=[1, 12], max_n_peaks=5, verbose=False)
        fm.fit(freqs, psd, freq_range)
        fm.plot(
            plot_peaks='shade-line',
            plot_aperiodic=True,
            freqs=freqs,
            power_spectrum=psd,
            freq_range=freq_range,
            plt_log=False,
            add_legend=True,
            save_fig=True,
            file_name=out_path.stem,
            file_path=str(out_path.parent),
            file_format='png'
        )
        return {
            "aperiodic_params": fm.aperiodic_params_.tolist(),
            "peak_params": [p.tolist() for p in fm.peak_params_],
            "r_squared": float(fm.r_squared_),
            "error": float(fm.error_),
            "n_peaks": int(fm.n_peaks_)
        }



    def compute_optimal_lag(x, y, dt=1.0, max_lag_ms=200):
        max_lag_steps = int(max_lag_ms / dt)
        x_norm = (x - np.mean(x)) / np.std(x)
        y_norm = (y - np.mean(y)) / np.std(y)
        corr = np.correlate(y_norm, x_norm, mode='full') / len(x)
        lags = np.arange(-len(x) + 1, len(x)) * dt
        center = len(corr) // 2
        lag_window = slice(center - max_lag_steps, center + max_lag_steps + 1)
        corr_window = corr[lag_window]
        lags_window = lags[lag_window]
        best_idx = np.argmax(np.abs(corr_window))
        best_lag = lags_window[best_idx]
        best_corr = corr_window[best_idx]
        return float(best_lag), float(best_corr)

    
    def entropy(probs):
        return -sum(p * log2(p) for p in probs if p > 0)
    
    def conditional_entropy(p_joint, p_conditioned):
        total_ce = 0.0
        for joint_key, p_joint_val in p_joint.items():
            if not isinstance(joint_key, tuple):
                continue
            # Se la chiave è 2D: (y1, y0)
            if len(joint_key) == 2:
                y1, y0 = joint_key
                cond_key = y0
            # Se la chiave è 3D: (y1, y0, x0)
            elif len(joint_key) == 3:
                y1, y0, x0 = joint_key
                cond_key = (y0, x0)
            else:
                continue  # ignora se ha forma inattesa
    
            if cond_key in p_conditioned and p_joint_val > 0:
                p_cond = p_conditioned[cond_key]
                total_ce += p_joint_val * log2(p_cond / p_joint_val)
        return total_ce

    
    def compute_transfer_entropy_binned(source, target, bins=16):
        # Discretizza i segnali
        source_d = np.digitize(source, np.histogram(source, bins=bins)[1]) - 1
        target_d = np.digitize(target, np.histogram(target, bins=bins)[1]) - 1
    
        # Costruisci sequenze temporali
        x_t = source_d[:-1]
        y_t = target_d[:-1]
        y_t1 = target_d[1:]
    
        # Calcola distribuzioni con Counter
        p_y1y = Counter(zip(y_t1, y_t))
        p_y = Counter(y_t)
        p_y1yx = Counter(zip(y_t1, y_t, x_t))
        p_yx = Counter(zip(y_t, x_t))
    
        # Normalizza
        total = len(y_t1)
        for d in [p_y1y, p_y, p_y1yx, p_yx]:
            for k in d:
                d[k] /= total
    
        # H(Y_{t+1} | Y_t)
        h_y1_given_y = conditional_entropy(p_y1y, p_y)
    
        # H(Y_{t+1} | Y_t, X_t)
        h_y1_given_yx = conditional_entropy(p_y1yx, p_yx)
    
        # Transfer Entropy: TE = H(Y_{t+1}|Y_t) - H(Y_{t+1}|Y_t, X_t)
        te = h_y1_given_y - h_y1_given_yx
        return float(te)

    

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

    def compute_cross_correlation(x, y):
        x = (x - np.mean(x)) / np.std(x)
        y = (y - np.mean(y)) / np.std(y)
        corr = np.correlate(x, y, mode='full') / len(x)
        max_corr = np.max(np.abs(corr))
        return float(max_corr)


    
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
            ax_main.set_ylim(0, 200)
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


            cross_corr = compute_cross_correlation(mossy_mean, purkinje_mean)
            te_forward = compute_transfer_entropy_binned(mossy_mean, purkinje_mean)
            te_backward = compute_transfer_entropy_binned(purkinje_mean, mossy_mean)


            lag_ms, lagged_corr = compute_optimal_lag(mossy_mean, purkinje_mean, dt=dt)

            # Sincronizzazione
            plv = compute_plv(mossy_mean, purkinje_mean)
            intra_corr_mossy = mean_pairwise_corr(mossy_array)
            intra_corr_purk = mean_pairwise_corr(purkinje_array)
            
            # fooof
            fs_hz = 1000.0 / dt
            fooof_out_mossy = analyze_and_plot_fooof(
            mossy_mean, fs=fs_hz,
            title=f"FOOOF – Mossy (k={kvalue:.3f})",
            out_path=output_dir / f"fooof_mossy_k_{kvalue:.3f}"
            )
            fooof_out_purkinje = analyze_and_plot_fooof(
            purkinje_mean, fs=fs_hz,
            title=f"FOOOF – Purkinje (k={kvalue:.3f})",
            out_path=output_dir / f"fooof_purkinje_k_{kvalue:.3f}"
            )

            # PSD media (usata anche da FOOOF ma qui salvata raw)
            f_mossy, psd_mossy = compute_population_psd(mossy_array, fs=fs_hz)
            f_purkinje, psd_purkinje = compute_population_psd(purkinje_array, fs=fs_hz)

            # Plot PSD overlapped: mossy + purkinje
            plt.figure(figsize=(10, 5))
            plt.plot(f_mossy, psd_mossy, label='Mossy fibers', color='dodgerblue', linewidth=2)
            plt.plot(f_purkinje, psd_purkinje, label='Purkinje cells', color='darkorange', linewidth=2)
            plt.xlim(0,100)
            plt.xlabel("Frequency (Hz)")
            plt.ylabel("Power")
            plt.title(f"PSD – Mossy vs Purkinje (k={kvalue:.3f})")
            plt.grid(True)
            plt.legend(loc='upper right')
            plt.tight_layout()
            plt.savefig(output_dir / f"psd_overlap_k_{kvalue:.3f}.png", dpi=300)
            plt.close()

            summary_data = {
                "k": float(f"{kvalue:.3f}"),
                "atrophy_percent": float(f"{atrophy:.1f}"),
                "mossy_mean_rate": float(np.mean(mossy_mean)),
                "purkinje_mean_rate": float(np.mean(purkinje_mean)),
                "cosine_similarity": float(f"{similarity:.4f}"),
                "cross_correlation": cross_corr,

                "cross_correlation_lag": {
                        "optimal_lag_ms": lag_ms,
                        "correlation_at_lag": lagged_corr
                    },

                "transfer_entropy": {
                    "mossy_to_purkinje": te_forward,
                    "purkinje_to_mossy": te_backward
                },
                "fooof": {
                    "mossy": fooof_out_mossy,
                    "purkinje": fooof_out_purkinje
                },
                "psd": {
                    "mossy": {
                        "frequencies_hz": f_mossy,
                        "power": psd_mossy
                    },
                    "purkinje": {
                        "frequencies_hz": f_purkinje,
                        "power": psd_purkinje
                    }
                },
                "plv": plv,
                "intra_group_corr": {
                    "mossy": intra_corr_mossy,
                    "purkinje": intra_corr_purk
                },
            }
            
            with open(output_dir / f"summary_k_{kvalue:.3f}.yaml", "w") as f:
                yaml.dump(summary_data, f, sort_keys=False)
            print(f"🧾 Saved YAML summary for k={kvalue:.3f}")



if __name__ == "__main__":
    main()
