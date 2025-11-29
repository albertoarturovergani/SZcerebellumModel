#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime
from scipy.signal import welch
from fooof import FOOOF
import argparse
import yaml


def load_resolution_sec(k_suffix):
    yaml_path = f'./configurations/mouse/nest/basal_vitro_{k_suffix}.yaml'
    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f)
    resolution_ms = config['simulations']['basal_activity']['resolution']
    return resolution_ms / 1000.0  # convert ms → sec


def generate_and_analyze_spike_train_extended(T, dt, aperiodic_exponent, mean_rate, seed, custom_datetime=None, k_suffix=None):
    if seed is not None:
        np.random.seed(seed)

    now_str = custom_datetime if custom_datetime else datetime.now().strftime("%Y%m%d_%H%M%S")
    folder_name = f"{now_str}_spikeSlope_{aperiodic_exponent:.3f}"
    save_dir = os.path.join("tempFolders", folder_name)
    os.makedirs(save_dir, exist_ok=True)

    def generate_rate(exponent):
        n = int(T / dt)
        time = np.arange(n) * dt
        freqs = np.fft.rfftfreq(n, d=dt)
        amplitudes = 1.0 / np.maximum(freqs, 1.0 / T)**(exponent / 2)
        random_phases = np.exp(1j * 2 * np.pi * np.random.rand(len(freqs)))
        spectrum = amplitudes * random_phases
        rate_raw = np.fft.irfft(spectrum, n=n)
        rate = rate_raw - np.min(rate_raw)
        rate = rate / np.mean(rate) * mean_rate
        return rate, time

    rate, time = generate_rate(aperiodic_exponent)
    fs = 1 / dt
    f, Pxx = welch(rate, fs=fs, nperseg=min(256, len(rate)))
    fm = FOOOF(peak_width_limits=[1, 6], aperiodic_mode='fixed', verbose=False)
    fm.fit(f, Pxx)

    p_spike = rate * dt
    spikes = np.random.rand(len(rate)) < p_spike
    spike_times = np.round(time[spikes] / dt) * dt * 1e3
    
    if k_suffix:
        resolution = load_resolution_sec(k_suffix)
        spike_times = np.round(spike_times / resolution) * resolution

    out_txt = os.path.join(save_dir, f"spike_times_original_{aperiodic_exponent:.3f}.txt")
    np.savetxt(out_txt, spike_times, fmt="%.6f")

    plt.figure()
    plt.plot(time, rate, label='Original')
    plt.title("Rate Function")
    plt.xlabel("Time (s)")
    plt.ylabel("Rate (Hz)")
    plt.legend()
    plt.savefig(os.path.join(save_dir, f"rate_comparison_extended_{aperiodic_exponent:.3f}.png"))
    plt.close()

    fig, ax = plt.subplots()
    fm.plot(ax=ax)
    fig.suptitle("FOOOF Fit")
    fig.tight_layout()
    plt.savefig(os.path.join(save_dir, f"fooof_plot_{aperiodic_exponent:.3f}.png"))
    plt.close(fig)

    print(f"✔ Spike times salvati in: {out_txt}")
    print(f"✔ Output folder: {save_dir}")
    return save_dir


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Genera spike train aperiodico e salva spike_times.")
    parser.add_argument('--T', type=float, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--slope', type=float, required=True)
    parser.add_argument('--rate', type=float, required=True)
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--datetime', type=str, default=None, help="Timestamp custom in formato YYYYMMDD_HHMMSS")
    parser.add_argument('--k_suffix', type=str, default=None, help="Usato per estrarre la resolution corretta dal file YAML")

    args = parser.parse_args()

    generate_and_analyze_spike_train_extended(
        T=args.T,
        dt=args.dt,
        aperiodic_exponent=round(args.slope, 3),
        mean_rate=args.rate,
        seed=args.seed,
        custom_datetime=args.datetime,
        k_suffix=args.k_suffix
    )
