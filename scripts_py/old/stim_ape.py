import numpy as np
import yaml
from scipy.signal import welch
from fooof import FOOOF
import matplotlib.pyplot as plt

def generate_colored_noise(exponent, n_samples, fs):
    freqs = np.fft.rfftfreq(n_samples, d=1/fs)
    spectrum = np.random.normal(0, 1, len(freqs)) + 1j * np.random.normal(0, 1, len(freqs))
    spectrum[1:] /= freqs[1:] ** (exponent / 2)
    spectrum[0] = 0
    signal = np.fft.irfft(spectrum, n=n_samples)
    return (signal - np.mean(signal)) / np.std(signal)

def generate_colored_spikes(duration_ms, fs=10000, exponent=2.0, threshold=2.5, target_rate=4):
    n_samples = int(duration_ms / 1000 * fs)
    colored_noise = generate_colored_noise(exponent, n_samples, fs)
    spike_train = (colored_noise > threshold).astype(int)
    spike_indices = np.where(spike_train)[0]
    expected_spikes = int((duration_ms / 1000) * target_rate)
    if len(spike_indices) > expected_spikes:
        spike_indices = np.sort(np.random.choice(spike_indices, size=expected_spikes, replace=False))
    spike_times = spike_indices * (1000 / fs)
    return spike_times, colored_noise

def compute_aperiodic_exponent(signal, fs=10000, plot=True):
    f, pxx = welch(signal, fs=fs, nperseg=fs*2)
    f = f[f > 0]
    pxx = pxx[:len(f)]
    if len(f) == 0 or len(pxx) == 0:
        raise ValueError("Power spectrum is empty.")
    fm = FOOOF(peak_width_limits=[1, 6], max_n_peaks=4)
    fm.fit(f, pxx)
    if plot:
        fm.plot()
        plt.show()
    return fm.get_params('aperiodic_params', 'exponent')

def update_yaml_with_spikes(yaml_path, new_spike_times, output_path):
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    spike_times_clean = [float(round(t, 4)) for t in new_spike_times]
    data['simulations']['basal_activity']['devices']['spike_generator']['spike_times'] = spike_times_clean
    with open(output_path, 'w') as f:
        yaml.dump(data, f, sort_keys=False)

def plot_rate(spikes_ms, fs=10000, duration_ms=None):
    if duration_ms is None:
        duration_ms = int(np.max(spikes_ms)) + 1
    time = np.arange(0, duration_ms, 1)
    bin_size = 100
    bins = np.arange(0, duration_ms + bin_size, bin_size)
    hist, _ = np.histogram(spikes_ms, bins=bins)
    psth = hist / (bin_size / 1000)
    spike_train = np.zeros_like(time)
    spike_indices = np.round(spikes_ms).astype(int)
    spike_indices = spike_indices[spike_indices < len(spike_train)]
    spike_train[spike_indices] = 1
    width = 100
    hw = width // 2
    kernel = np.concatenate([np.linspace(0, 1, hw, endpoint=False), np.linspace(1, 0, hw)])
    kernel /= kernel.sum()
    rate_smoothed = np.convolve(spike_train, kernel, mode='same') * fs
    total_spikes = len(spikes_ms)
    mean_rate_hz = total_spikes / (duration_ms / 1000)

    fig, ax = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    ax[0].bar(bins[:-1], psth, width=bin_size, align='edge', color='skyblue', edgecolor='black')
    ax[0].axhline(mean_rate_hz, color='gray', linestyle='--', label=f"Mean rate = {mean_rate_hz:.2f} Hz")
    ax[0].set_title("PSTH binned (100 ms)")
    ax[0].set_ylabel("Rate (Hz)")
    ax[0].legend()
    ax[1].plot(time, rate_smoothed, color='tomato')
    ax[1].axhline(mean_rate_hz, color='gray', linestyle='--', label=f"Mean rate = {mean_rate_hz:.2f} Hz")
    ax[1].set_title("Smoothed rate (triangular kernel 100 ms)")
    ax[1].set_ylabel("Rate (Hz)")
    ax[1].set_xlabel("Time (ms)")
    ax[1].legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    np.random.seed(42)
    duration_ms = 5000
    fs = 10000
    exponent_target = 0 #1.2  # ⬅️ da 0 (white) a 2 (brown)
    threshold = 2.0
    target_rate = 4
    yaml_input = "./configurations/mouse/nest/basal_vitro_withSpikes_2.yaml"
    yaml_output = "./configurations/mouse/nest/basal_vitro_withSpikes_personalized.yaml"

    for attempt in range(10):
        spikes_ms, signal = generate_colored_spikes(
            duration_ms, fs=fs,
            exponent=exponent_target,
            threshold=threshold,
            target_rate=target_rate
        )
        try:
            exponent = compute_aperiodic_exponent(signal, fs=fs, plot=True)
            print(f"Attempt {attempt+1}: Target exponent = {exponent_target} — Estimated = {exponent:.2f}")
            if abs(exponent - exponent_target) < 0.3:
                break
        except Exception as e:
            print(f"Attempt {attempt+1} failed: {e}")

    update_yaml_with_spikes(yaml_input, spikes_ms, yaml_output)
    print(f"Updated YAML saved to {yaml_output}")

    plot_rate(spikes_ms, fs=fs)
