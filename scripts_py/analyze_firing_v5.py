#!/usr/bin/env python3

def main():
    import gc
    import neo
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from pathlib import Path
    from collections import defaultdict, Counter
    from scipy.signal import convolve, welch, hilbert
    from scipy.signal.windows import triang
    from scipy.stats import pearsonr
    from math import log2
    import argparse
    import yaml
    import re
    from fooof import FOOOF

    def convert_for_yaml(obj):
        if isinstance(obj, dict):
            return {convert_for_yaml(k): convert_for_yaml(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_for_yaml(i) for i in obj]
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.generic,)):
            return obj.item()
        else:
            return obj

    def compute_population_psd(array, fs):
        if array.size == 0:
            return [], []
        psds = []
        for unit in array:
            f, p = welch(unit, fs=fs, nperseg=1024)
            psds.append(p)
        psd_mean = np.mean(psds, axis=0)
        return f.tolist(), psd_mean.tolist()

    def compute_plv(sig1, sig2):
        phase1 = np.angle(hilbert(sig1))
        phase2 = np.angle(hilbert(sig2))
        phase_diff = phase1 - phase2
        plv = np.abs(np.mean(np.exp(1j * phase_diff)))
        return float(plv)

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
            if len(joint_key) == 2:
                y1, y0 = joint_key
                cond_key = y0
            elif len(joint_key) == 3:
                y1, y0, x0 = joint_key
                cond_key = (y0, x0)
            else:
                continue
            if cond_key in p_conditioned and p_joint_val > 0:
                p_cond = p_conditioned[cond_key]
                total_ce += p_joint_val * log2(p_cond / p_joint_val)
        return total_ce

    def compute_transfer_entropy_binned(source, target, bins=16):
        source_d = np.digitize(source, np.histogram(source, bins=bins)[1]) - 1
        target_d = np.digitize(target, np.histogram(target, bins=bins)[1]) - 1
        x_t = source_d[:-1]
        y_t = target_d[:-1]
        y_t1 = target_d[1:]
        p_y1y = Counter(zip(y_t1, y_t))
        p_y = Counter(y_t)
        p_y1yx = Counter(zip(y_t1, y_t, x_t))
        p_yx = Counter(zip(y_t, x_t))
        total = len(y_t1)
        for d in [p_y1y, p_y, p_y1yx, p_yx]:
            for k in d:
                d[k] /= total
        h_y1_given_y = conditional_entropy(p_y1y, p_y)
        h_y1_given_yx = conditional_entropy(p_y1yx, p_yx)
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

    parser = argparse.ArgumentParser(description="Process neuronal firing rates.")
    #parser.add_argument("--path", nargs="?", default=".", help="Main directory")
    parser.add_argument('--path', type=str, required=True, help='Percorso alla cartella dati')
    parser.add_argument("--sigma", type=float, default=100.0, help="Sigma for kernel smoothing in firing rate estimation")
    args = parser.parse_args()
    #root_path = Path(args.root)
    root_path = Path(args.path)  
    filter_mode = True
    trim_ratio = 0.005
    dt = 0.1
    sigma = args.sigma

    for result_dir in root_path.glob('*_deg_*/3.function'):
        deg_folder = result_dir.parent.name
        match = re.match(r'(\d{3})_deg_([0-9]+\.[0-9]+)', deg_folder)
        if not match:
            print(f"⚠️ Skipping unrecognized folder: {deg_folder}")
            continue
        int_part, atrophy_str = match.groups()
        atrophy_percent = float(atrophy_str)
        kvalue = atrophy_percent #round((1 - atrophy_percent) / 100, 3)

        output_dir_k = result_dir / f"plots_firing_rate"
        output_dir_k.mkdir(parents=True, exist_ok=True)

        spike_rows = []
        for file in result_dir.glob("*.nio"):
            reader = neo.io.NixIO(filename=str(file), mode="ro")
            block = reader.read_block()
            for seg_idx, segment in enumerate(block.segments):
                for st_idx, spiketrain in enumerate(segment.spiketrains):
                    spike_times = spiketrain.times.rescale("ms").magnitude
                    #device = spiketrain.annotations.get("device", "unknown")
                    device = spiketrain.annotations.get("device", "unknown").removesuffix("_record")

                    senders = spiketrain.annotations.get("senders", None)
                    senders = np.atleast_1d(senders) if senders is not None else None
                    for idx, time in enumerate(spike_times):
                        sender_id = senders[idx] if senders is not None and idx < len(senders) else f"unit{st_idx}"
                        spike_rows.append({
                            "time_ms": time,
                            "sender_id": sender_id,
                            "pop": device, #device.split('_')[0],
                            "segment": seg_idx,
                            "file": file.name
                        })
            del reader, block
            gc.collect()

        df_spikes = pd.DataFrame(spike_rows)
        if df_spikes.empty:
            continue
        if filter_mode:
            df_spikes = df_spikes[df_spikes['pop'].isin(['mossy_fibers', 'purkinje'])]

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

        mossy_fibers_rates = [rate for pop, rate in zip(firing_data["pop"], firing_data["rate"]) if pop == "mossy_fibers"]
        purkinje_rates = [rate for pop, rate in zip(firing_data["pop"], firing_data["rate"]) if pop == "purkinje"]

        if not mossy_fibers_rates or not purkinje_rates:
            print(f"⚠️ Skipping k={kvalue}: insufficient data.")
            continue

        t = firing_data["time"][0]
        n = len(t)
        start_idx = int(n * trim_ratio)
        end_idx = int(n * (1 - trim_ratio))
        t = t[start_idx:end_idx]

        mossy_fibers_array = np.vstack(mossy_fibers_rates)[:, start_idx:end_idx]
        purkinje_array = np.vstack(purkinje_rates)[:, start_idx:end_idx]
        mossy_fibers_mean = np.mean(mossy_fibers_array, axis=0)
        purkinje_mean = np.mean(purkinje_array, axis=0)

        fig = plt.figure(figsize=(10, 4.5))
        ax = fig.add_subplot(111)
        ax.plot(t, mossy_fibers_mean, label="mossy_fibers fibers", linewidth=1.5)
        ax.plot(t, purkinje_mean, label="Purkinje cells", linewidth=1.5)
        ax.set_ylim(0, 200)
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Firing rate (Hz)")
        ax.set_title(f"Firing Rate – k={kvalue:.3f}")
        ax.grid(True)
        ax.legend()
        fig.savefig(output_dir_k / f"fr_plot_{kvalue:.3f}.png", dpi=300)
        plt.close(fig)

        similarity = np.dot(mossy_fibers_mean, purkinje_mean) / (np.linalg.norm(mossy_fibers_mean) * np.linalg.norm(purkinje_mean))
        cross_corr = compute_cross_correlation(mossy_fibers_mean, purkinje_mean)
        lag_ms, lagged_corr = compute_optimal_lag(mossy_fibers_mean, purkinje_mean, dt=dt)
        plv = compute_plv(mossy_fibers_mean, purkinje_mean)
        te_forward = compute_transfer_entropy_binned(mossy_fibers_mean, purkinje_mean)
        te_backward = compute_transfer_entropy_binned(purkinje_mean, mossy_fibers_mean)

        f_mossy_fibers, psd_mossy_fibers = compute_population_psd(mossy_fibers_array, fs=1000/dt)
        f_purk, psd_purk = compute_population_psd(purkinje_array, fs=1000/dt)

        fooof_out_mossy_fibers = analyze_and_plot_fooof(mossy_fibers_mean, 1000/dt, f"FOOOF – mossy_fibers (k={kvalue:.3f})", output_dir_k / f"fooof_mossy_fibers_{kvalue:.3f}")
        fooof_out_purkinje = analyze_and_plot_fooof(purkinje_mean, 1000/dt, f"FOOOF – Purkinje (k={kvalue:.3f})", output_dir_k / f"fooof_purkinje_{kvalue:.3f}")

        summary = {
            "k": kvalue,
            "time_window" : "Overall Time",
            "atrophy_percent": 100*(1-kvalue),
            "cosine_similarity": similarity,
            "cross_correlation": cross_corr,
            "lag_ms": lag_ms,
            "lagged_correlation": lagged_corr,
            "plv": plv,
            "transfer_entropy": {
                "mossy_fibers_to_purkinje": te_forward,
                "purkinje_to_mossy_fibers": te_backward,
            },
            "fooof": {"mossy_fibers": fooof_out_mossy_fibers, "purkinje": fooof_out_purkinje},
            #"psd": {"mossy_fibers": {"f": f_mossy_fibers, "psd": psd_mossy_fibers}, "purkinje": {"f": f_purk, "psd": psd_purk}},
        }

        with open(output_dir_k / f"summary_{kvalue:.3f}.yaml", 'w') as f:
            yaml.safe_dump(convert_for_yaml(summary), f)

if __name__ == "__main__":
    main()