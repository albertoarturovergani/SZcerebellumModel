import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from collections import defaultdict, Counter
from pathlib import Path
from fooof import FOOOF
import matplotlib.pyplot as plt
from math import log2
import pickle
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
import seaborn as sns
import argparse

#base_dir = Path("vbt_temp_results_2025-04-21_21-00-17-exp/")

import matplotlib.pyplot as plt

plt.rcParams.update({
    # Figure
    'figure.figsize': (8, 4),         # Figure size in inches
    'figure.dpi': 150,                # Lower DPI for faster rendering, still sharp
    'savefig.dpi': 300,               # Save at high resolution
    'savefig.format': 'png',
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,

    # Axes
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'axes.labelweight': 'bold',
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
    'grid.color': '#cccccc',
    'axes.spines.top': False,
    'axes.spines.right': False,

    # Lines
    'lines.linewidth': 1.5,
    'lines.markersize': 4,
    'lines.marker': '',  # Vuoto di default, lo specifichi nei plot se serve

    # Font
    'font.size': 12,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],

    # Legend
    'legend.fontsize': 10,
    'legend.title_fontsize': 12,
    'legend.frameon': True,
    'legend.loc': 'best',

    # Ticks
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'xtick.major.pad': 4,
    'ytick.major.pad': 4,

    # Colors
    'axes.prop_cycle': plt.cycler(color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']),
    'image.cmap': 'bwr'
})

import os
import yaml
import matplotlib.pyplot as plt
from pathlib import Path

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
                    "start": float(stim["start"]),
                    "stop": float(stim["stop"]),
                    "rate": float(stim["rate"]),
                    "type": "burst" if "burst" in label else "stim"
                })
    except Exception as e:
        print(f"⚠️ Errore parsing {yaml_path}: {e}")
        return False

    if not stimuli:
        print(f"⚠️ Nessuno stimolo trovato in {yaml_path}")
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
    ax2.scatter(midpoints, rates, marker="o", color="darkgreen")
    ax2.set_title("Rate vs Time")
    ax2.set_xlabel("Time (ms)")
    ax2.set_ylabel("Rate (Hz)")
    ax2.grid(True)

    plt.tight_layout()
    out_path = yaml_path.parent / "stimulus_timeline.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"✅ Salvato: {out_path}")
    return True


def process_reconstruction_yamls(base_dir: Path):
    yaml_paths = list(base_dir.rglob("*_deg_*/2.stimulus/*.yaml"))
    if not yaml_paths:
        print(f"⚠️ Nessun file YAML trovato sotto {base_dir}")
        return

    count = 0
    for yaml_path in yaml_paths:
        if generate_stimulus_timeline(yaml_path):
            count += 1

    print(f"\n🎯 Generati {count} plot di timeline da {len(yaml_paths)} file YAML.")

import matplotlib.pyplot as plt
import numpy as np

def plot_firing_rate_panels(m_mean, p_mean, mossy_fibers, purk, stim, dt, start, end, atrophy, stim_to_mossy_fibers_cosine, stim_to_purk_cosine, k, base_dir):
    t = np.arange(len(m_mean)) * dt
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 9), sharex=True, sharey=True)

    q_low, q_high = 0.25, 0.75

    # --- mossy_fibers ---
    m_min = np.quantile(mossy_fibers[:, start:end], q_low, axis=0)[:len(m_mean)]
    m_max = np.quantile(mossy_fibers[:, start:end], q_high, axis=0)[:len(m_mean)]
    ax1.plot(t, m_mean, color="black", label="Mean firing rate", linestyle='--', linewidth=0.75)
    ax1.fill_between(t, m_min, m_max, color="gray", alpha=0.2, label=f"{int(q_low*100)}–{int(q_high*100)}% quantile range")
    if stim is not None:
        ax1.plot(t, stim, color="red", alpha=0.5, linestyle='--', linewidth=0.75, label="Stimulus profile")
    ax1.set_title(f"mossy_fibers firing rate – Atrophy {atrophy:.1f}% – cos(stim,m)={stim_to_mossy_fibers_cosine:.2f}")
    ax1.legend(loc="upper right")
    ax1.grid(True)

    # --- Purkinje ---
    p_min = np.quantile(purk[:, start:end], q_low, axis=0)[:len(p_mean)]
    p_max = np.quantile(purk[:, start:end], q_high, axis=0)[:len(p_mean)]
    ax2.plot(t, p_mean, color="black", label="Mean firing rate", linestyle='--', linewidth=0.75)
    ax2.fill_between(t, p_min, p_max, color="gray", alpha=0.2, label=f"{int(q_low*100)}–{int(q_high*100)}% quantile range")
    if stim is not None:
        ax2.plot(t, stim, color="red", alpha=0.5, linestyle='--', linewidth=0.75, label="Stimulus profile")
    ax2.set_title(f"Purkinje firing rate – Atrophy {atrophy:.1f}% – cos(stim,p)={stim_to_purk_cosine:.2f}")
    ax2.set_xlabel("Time (ms)")
    ax2.set_ylabel("Firing Rate (Hz)")
    ax2.legend(loc="upper right")
    ax2.grid(True)

    # --- mossy_fibers vs Purkinje ---
    ax3.plot(t, m_mean, label="mossy_fibers", linestyle='-', linewidth=0.75, color="orange")
    ax3.plot(t, p_mean, label="Purkinje", linestyle='--', linewidth=0.75, color="blue")
    if stim is not None:
        ax3.plot(t, stim, color="red", alpha=0.5, linestyle='--', linewidth=0.75, label="Stimulus profile")
    from numpy import dot, linalg
    cos_mp = float(dot(m_mean, p_mean) / (linalg.norm(m_mean) * linalg.norm(p_mean)))
    ax3.set_title(f"mossy_fibers vs Purkinje – cos(m,p)={cos_mp:.2f}")
    ax3.set_xlabel("Time (ms)")
    ax3.legend(loc="upper right")
    ax3.grid(True)

    total_time = len(m_mean) * dt
    ax1.set_xlim(500, total_time - 500)

    plt.tight_layout()
    out_path = f"./{base_dir}/fr_vs_time_k{k:.3f}.png"
    fig.savefig(out_path, dpi=300)
    plt.close()
    print(f"[INFO] Salvato: {out_path}")

def printSynTime():
    import yaml
    from pathlib import Path
    
    # Percorso al file YAML
    yaml_path = Path("/home/alberto/cerebellum/cerebellum_APRIL25/saccades1_results_2025-04-17_19-10-58-exp/results_1.000/reconstruction/basal_vitro_1.000.yaml")
    
    with open(yaml_path, "r") as f:
        data = yaml.safe_load(f)
    
    # Estrai la sezione delle connessioni
    conn_models = data["simulations"]["basal_activity"]["connection_models"]
    
    # Definizione del cammino MF → PC
    path = [
        "mossy_fibers_fibers_to_glomerulus",
        "glomerulus_to_granule",
        "parallel_fiber_to_purkinje"
    ]
    
    # Estrazione dei delay lungo il cammino
    total_delay = 0.0
    detailed_path = []
    
    for conn in path:
        if conn not in conn_models:
            raise ValueError(f"❌ Connessione mancante nel file: {conn}")
        delay = conn_models[conn]["synapse"].get("delay", 0.0)
        total_delay += delay
        detailed_path.append((conn, delay))
    
    # Output
    print("Cammino MF → PC funzionale:")
    for conn, d in detailed_path:
        print(f" - {conn} → {d} ms")
    
    print(f"Delay totale MF → PC = {total_delay:.2f} ms")



def compute_firing_rate(spike_times, duration, dt=1.0, sigma=100):
    time = np.arange(0, duration, dt)
    spike_train = np.zeros_like(time)
    indices = (np.array(spike_times) / dt).astype(int)
    indices = indices[indices < len(spike_train)]
    spike_train[indices] = 1
    kernel_size = max(int(2 * sigma / dt), 3)
    kernel = triang(kernel_size)
    kernel /= kernel.sum()
    rate = convolve(spike_train, kernel, mode='same') * (1000.0 / dt)
    return time, rate

def compute_population_psd(array, fs):
    if array.size == 0: return [], []
    psds = [welch(unit, fs=fs, nperseg=1024)[1] for unit in array]
    return welch(array[0], fs=fs, nperseg=1024)[0].tolist(), np.mean(psds, axis=0).tolist()

def analyze_fooof(signal, fs, freq_range=(1, 100)):
    freqs, psd = welch(signal, fs=fs, nperseg=1024)
    fm = FOOOF(peak_width_limits=[1, 12], max_n_peaks=5, verbose=False)
    fm.fit(freqs, psd, freq_range)
    return {
        "aperiodic_params": fm.aperiodic_params_.tolist(),
        "peak_params": [p.tolist() for p in fm.peak_params_],
        "r_squared": float(fm.r_squared_),
        "error": float(fm.error_),
        "n_peaks": int(fm.n_peaks_)
    }

def compute_cross_correlation_old(x, y):
    x = (x - np.mean(x)) / np.std(x)
    y = (y - np.mean(y)) / np.std(y)
    return float(np.max(np.abs(np.correlate(x, y, mode='full')) / len(x)))

def compute_cross_corr_and_lag(x, y, dt=1.0):
    from scipy.signal import correlate
    from scipy.stats import zscore

    x, y = zscore(x), zscore(y)
    corr = correlate(y, x, mode='full') / len(x)  # Normalizzazione qui
    lags = np.arange(-len(x) + 1, len(x)) * dt
    idx_max = np.argmax(np.abs(corr))
    return float(corr[idx_max]), float(lags[idx_max])


def compute_optimal_lag(x, y, dt=1.0, max_lag_ms=200):
    max_lag_steps = int(max_lag_ms / dt)
    x, y = (x - np.mean(x)) / np.std(x), (y - np.mean(y)) / np.std(y)
    corr = np.correlate(y, x, mode='full') / len(x)
    lags = np.arange(-len(x)+1, len(x)) * dt
    center = len(corr) // 2
    window = slice(center - max_lag_steps, center + max_lag_steps + 1)
    best_idx = np.argmax(np.abs(corr[window]))
    return float(lags[window][best_idx]), float(corr[window][best_idx])

def compute_plv(sig1, sig2):
    phase1, phase2 = np.angle(hilbert(sig1)), np.angle(hilbert(sig2))
    return float(np.abs(np.mean(np.exp(1j * (phase1 - phase2)))))

def mean_pairwise_corr(array):
    n = array.shape[0]
    return float(np.nanmean([pearsonr(array[i], array[j])[0]
                             for i in range(n) for j in range(i+1, n)])) if n > 1 else np.nan

def entropy(probs):
    return -sum(p * log2(p) for p in probs if p > 0)

def conditional_entropy(p_joint, p_conditioned):
    total_ce = 0.0
    for joint_key, p_joint_val in p_joint.items():
        if not isinstance(joint_key, tuple): continue
        cond_key = joint_key[1] if len(joint_key) == 2 else (joint_key[1], joint_key[2]) if len(joint_key) == 3 else None
        if cond_key in p_conditioned and p_joint_val > 0:
            total_ce += p_joint_val * log2(p_conditioned[cond_key] / p_joint_val)
    return total_ce

def compute_transfer_entropy_binned_old(source, target, bins=16):
    s_d = np.digitize(source, np.histogram(source, bins=bins)[1]) - 1
    t_d = np.digitize(target, np.histogram(target, bins=bins)[1]) - 1
    x_t, y_t, y_t1 = s_d[:-1], t_d[:-1], t_d[1:]
    p_y1y, p_y = Counter(zip(y_t1, y_t)), Counter(y_t)
    p_y1yx, p_yx = Counter(zip(y_t1, y_t, x_t)), Counter(zip(y_t, x_t))
    total = len(y_t1)
    for d in [p_y1y, p_y, p_y1yx, p_yx]: d.update({k: v/total for k, v in d.items()})
    return conditional_entropy(p_y1y, p_y) - conditional_entropy(p_y1yx, p_yx)


from collections import Counter
import numpy as np
#bins = int(np.sqrt(len(source)))  # regola di Sturges modificata
def compute_transfer_entropy_binned(source, target, bins=None, normalize=True):
    
    if bins is None: bins = int(np.sqrt(len(source)))  # Regola empirica adattiva (sqrt(N))
        
    s_d = np.digitize(source, np.histogram(source, bins=bins)[1]) - 1
    t_d = np.digitize(target, np.histogram(target, bins=bins)[1]) - 1
    x_t, y_t, y_t1 = s_d[:-1], t_d[:-1], t_d[1:]

    p_y1 = Counter(y_t1)
    p_y1y = Counter(zip(y_t1, y_t))
    p_y1yx = Counter(zip(y_t1, y_t, x_t))
    p_y = Counter(y_t)
    p_yx = Counter(zip(y_t, x_t))

    total = len(y_t1)
    for d in [p_y1, p_y1y, p_y1yx, p_y, p_yx]:
        for k in d:
            d[k] /= total

    def conditional_entropy(joint, marginal):
        return -sum(p * np.log2(p / marginal[k[1:]]) for k, p in joint.items() if p > 0 and marginal[k[1:]] > 0)

    TE = conditional_entropy(p_y1y, p_y) - conditional_entropy(p_y1yx, p_yx)

    if normalize:
        H_y1 = -sum(p * np.log2(p) for p in p_y1.values() if p > 0)
        return TE / H_y1 if H_y1 > 0 else 0.0
    else:
        return TE

def plot_tf(df_all):

    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import pearsonr, wilcoxon
    from pathlib import Path
    
    # Calcolo Atrophy Factor se non presente
    if "Atrophy Factor (%)" not in df_all.columns:
        df_all["Atrophy Factor (%)"] = (1 - df_all["k"]) * 100
    
    # Colori distinti per le due direzioni
    color_m2p = "#1f77b4"  # blu
    color_p2m = "#ff7f0e"  # arancio
    color_diff = "#2ca02c"  # verde per Δ
    
    sns.set(style="whitegrid", context="paper", font_scale=1.2)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), sharex=True)
    
    # ========= LEFT PANEL: MF→PC vs PC→MF ========= #
    for cond, group in df_all.groupby("condition"):
        # MF → PC
        r1, p1 = pearsonr(group["Atrophy Factor (%)"], group["te_mossy_fibers_to_purk"])
        label1 = f"{cond} MF→PC (r={r1:.2f}, p={p1:.1e})"
        sns.lineplot(
            data=group,
            x="Atrophy Factor (%)",
            y="te_mossy_fibers_to_purk",
            label=label1,
            marker="o",
            linewidth=2,
            color=color_m2p,
            ax=ax1
        )
        sns.regplot(
            data=group,
            x="Atrophy Factor (%)",
            y="te_mossy_fibers_to_purk",
            scatter=False,
            order=1,
            color=color_m2p,
            line_kws={"linestyle": "--", "alpha": 0.5},
            ax=ax1
        )
    
        # PC → MF
        r2, p2 = pearsonr(group["Atrophy Factor (%)"], group["te_purk_to_mossy_fibers"])
        label2 = f"{cond} PC→MF (r={r2:.2f}, p={p2:.1e})"
        sns.lineplot(
            data=group,
            x="Atrophy Factor (%)",
            y="te_purk_to_mossy_fibers",
            label=label2,
            marker="s",
            linewidth=2,
            color=color_p2m,
            ax=ax1
        )
        sns.regplot(
            data=group,
            x="Atrophy Factor (%)",
            y="te_purk_to_mossy_fibers",
            scatter=False,
            order=1,
            color=color_p2m,
            line_kws={"linestyle": "--", "alpha": 0.5},
            ax=ax1
        )
    
    # Test statistico tra le due distribuzioni (paired)
    try:
        stat_w, pval_w = wilcoxon(
            df_all["te_mossy_fibers_to_purk"],
            df_all["te_purk_to_mossy_fibers"]
        )
        print(stat_w)
        ax1.text(
            0.05, 0.95,
            f"Wilcoxon p = {pval_w:.2e}",
            transform=ax1.transAxes,
            ha="left",
            va="top",
            fontsize=10,
            bbox=dict(facecolor="white", edgecolor="gray", boxstyle="round,pad=0.2")
        )
    except Exception as e:
        print("⚠️ Errore test statistico:", e)
    
    ax1.set_title("Transfer Entropy MF↔PC", fontsize=13)
    ax1.set_xlabel("Atrophy Factor (%)")
    ax1.set_ylabel("Transfer Entropy (bits)")
    ax1.legend(title="Condition | Direction", fontsize=9)
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # ========= RIGHT PANEL: DIFFERENCE ========= #
    df_all["ΔTE"] = df_all["te_mossy_fibers_to_purk"] - df_all["te_purk_to_mossy_fibers"]
    r_diff, p_diff = pearsonr(df_all["Atrophy Factor (%)"], df_all["ΔTE"])
    
    sns.scatterplot(
        data=df_all,
        x="Atrophy Factor (%)",
        y="ΔTE",
        ax=ax2,
        color=color_diff,
        marker="D",
        label=f"ΔTE (r={r_diff:.2f}, p={p_diff:.1e})"
    )
    sns.regplot(
        data=df_all,
        x="Atrophy Factor (%)",
        y="ΔTE",
        scatter=False,
        color=color_diff,
        line_kws={"linestyle": "--", "alpha": 0.6},
        ax=ax2
    )
    
    ax2.axhline(0, linestyle="--", color="black", linewidth=1, alpha=0.5)
    ax2.set_title("Δ Transfer Entropy (MF→PC − PC→MF)", fontsize=13)
    ax2.set_ylabel("ΔTE (bits)")
    ax2.set_xlabel("Atrophy Factor (%)")
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    
    # Salvataggio
    output_path = Path(f"./{base_dir}/TE_comparison_with_diff.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=300)

def plot_cross_correlation(df_all):

    from pathlib import Path
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import pearsonr
    
    # Calcolo Atrophy Factor
    if "Atrophy Factor (%)" not in df_all.columns:
        df_all["Atrophy Factor (%)"] = (1 - df_all["k"]) * 100
    
    # Setup stile
    sns.set(style="whitegrid", context="paper", font_scale=1.2)
    fig, ax1 = plt.subplots(figsize=(10, 5))
    
    # Secondo asse Y
    ax2 = ax1.twinx()
    
    # Colori
    color1 = "k"
    color2 = "r"
    
    # Cross-correlation (asse sinistro)
    for cond, group in df_all.groupby("condition"):
        # Escludi punto con atrofia massima
        max_atrophy = group["Atrophy Factor (%)"].max()
        group_filtered = group[group["Atrophy Factor (%)"] < max_atrophy]
    
        r1, p1 = pearsonr(group_filtered["Atrophy Factor (%)"], group_filtered["cross_corr"])
        label1 = f"{cond} – Corr (r={r1:.2f}, p={p1:.1e})"
        sns.lineplot(
            data=group,
            x="Atrophy Factor (%)",
            y="cross_corr",
            marker="o",
            linewidth=2,
            color=color1,
            label=label1,
            ax=ax1
        )
    
    # Optimal lag (asse destro)
    for cond, group in df_all.groupby("condition"):
        max_atrophy = group["Atrophy Factor (%)"].max()
        group_filtered = group[group["Atrophy Factor (%)"] < max_atrophy]
    
        r2, p2 = pearsonr(group_filtered["Atrophy Factor (%)"], group_filtered["optimal_lag_ms"])
        label2 = f"{cond} – Lag (r={r2:.2f}, p={p2:.1e})"
        sns.lineplot(
            data=group,
            x="Atrophy Factor (%)",
            y="optimal_lag_ms",
            marker="s",
            linewidth=2,
            color=color2,
            label=label2,
            ax=ax2
        )
    
    # Regressione globale (facoltativa)
    sns.regplot(
        data=df_all,
        x="Atrophy Factor (%)",
        y="cross_corr",
        order=2,
        scatter=False,
        ax=ax1,
        label=None,
        color=color1,
        line_kws={"linestyle": "--", "alpha": 0}
    )
    sns.regplot(
        data=df_all,
        x="Atrophy Factor (%)",
        y="optimal_lag_ms",
        order=2,
        scatter=False,
        ax=ax2,
        label=None,
        color=color2,
        line_kws={"linestyle": "--", "alpha": 0}
    )
    
    # Linee di riferimento
    ax2.axhline(y=7, color='r', linestyle='-', linewidth=1)
    ax2.axhline(y=0, color='k', linestyle='--', linewidth=1)
    # 📌 Minimal MF → PC synaptic delay annotation
    min_atrophy = df_all["Atrophy Factor (%)"].min()
    ax2.text(
        x=min_atrophy - 1,
        y=5,
        s=(
            "7 ms (min MF→PC delay)\n"
            "Pathway:\n"
            "  MF → glomerulus: 1 ms\n"
            "  Glomerulus → granule: 1 ms\n"
            "  Parallel fiber → PC: 5 ms"
        ),
        fontsize=9,
        va="top",
        ha="left",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="gray")
    )

    
    # Assi e titolo
    ax1.set_ylabel("Cross-Correlation (MF → PC)", color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)
    
    ax2.set_ylabel("Optimal Lag (MF → PC) [ms]", color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)
    
    ax1.set_xlabel("Atrophy Factor (%)")
    ax1.set_title("Cross-Correlation & Optimal Lag (MF → PC) vs Atrophy (excl. max atrophy)")
    
    ax1.grid(True)
    plt.tight_layout()
    
    # Salvataggio
    output_path = f"./{base_dir}/twin_corr_lag_plot_excl_max.png"
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.show()

def plot_firing_rate_vs_k(df_all, save_path=None):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from scipy.stats import spearmanr

    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    palette = {"mossy_fibers": "#1f77b4", "purkinje": "#ff7f0e"}
    spearman_vals = {}

    for pop in ["mossy_fibers", "purkinje"]:
        mean_col = f"{pop}_mean_rate"
        sem_col = f"{pop}_sem_rate"

        # Prepara il DataFrame
        sub_df = df_all[["Atrophy Factor (%)", mean_col, sem_col]].copy()
        sub_df.rename(columns={mean_col: "firing_rate", sem_col: "sem"}, inplace=True)

        # Scatter con barre di errore
        ax.errorbar(
            sub_df["Atrophy Factor (%)"], sub_df["firing_rate"],
            yerr=sub_df["sem"],
            fmt='o', capsize=4,
            color=palette[pop],
            label=f"{pop.capitalize()} (mean ± SEM)",
            alpha=0.9,
            markersize=6
        )

        # Fit lineare
        sns.regplot(
            data=sub_df,
            x="Atrophy Factor (%)",
            y="firing_rate",
            scatter=False,
            color=palette[pop],
            line_kws={"label": f"{pop.capitalize()} (fit)"},
            ax=ax
        )

        # Correlazione Spearman
        rho, pval = spearmanr(sub_df["Atrophy Factor (%)"], sub_df["firing_rate"])
        spearman_vals[pop] = (rho, pval)

    # Titoli e assi
    ax.set_title("Firing Rate vs Atrophy Factor")
    ax.set_xlabel("Atrophy Factor (%)")
    ax.set_ylabel("Firing Rate (Hz)")
    ax.grid(True)

    # Legenda combinata
    handles, labels = ax.get_legend_handles_labels()
    extra_labels = [
        f"{pop.capitalize()} ρ = {rho:.2f}, p = {pval:.2e}"
        for pop, (rho, pval) in spearman_vals.items()
    ]
    extra_handles = [plt.Line2D([], [], linestyle='none', label=txt) for txt in extra_labels]

    # Corretto: passiamo sia handles che labels
    ax.legend(handles + extra_handles, labels + extra_labels, loc="upper right", fontsize=10)


    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"✅ Figura salvata in: {save_path}")
    plt.show()

def plot_stimulus_profile_per_k(df_all, base_dir):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from pathlib import Path

    sns.set(style="whitegrid", context="talk", font_scale=1.1)

    for _, row in df_all.iterrows():
        k = row["k"]
        stim = np.array(row["stim_profile"])

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(stim, color="orange", linewidth=2)

        ax.set_title(f"Stimulus Profile – k={k:.3f}", fontsize=14)
        ax.set_xlabel("Time (ms)", fontsize=12)
        ax.set_ylabel("Stimulus Rate (Hz)", fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.6)

        plt.tight_layout()

        output_path = Path(base_dir) / f"stim_profile_k{k:.3f}.png"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")
        plt.close()
        print(f"✅ Salvato: {output_path}")

def plot_stimulus_profile(df_all):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    
    # Setup stile
    sns.set(style="whitegrid", context="talk", font_scale=1.1)
    
    # Estrai stimolo
    stim = np.array(df_all["stim_profile"].iloc[0])
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(stim, color="orange", linewidth=2)
    
    # Titoli e assi
    ax.set_title("Stimulus Profile", fontsize=16)
    ax.set_xlabel("Time (ms)", fontsize=12)
    ax.set_ylabel("Stimulus rate (Hz)", fontsize=12)
    ax.grid(True, linestyle="--", alpha=0.6)
    
    # Togli i margini bianchi inutili
    plt.tight_layout()
    
    output_path = f"./{base_dir}/stim_profile.png"
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.show()

def plot_time_series(df_all, pop="mossy_fibers", dt=0.1, alpha=0.15, linewidth=1.0, dpi=300):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    from pathlib import Path

    # Ordinamento per atrofia
    df_sorted = df_all.sort_values("atrophy_percent").reset_index(drop=True)
    atrophy_vals = df_sorted["atrophy_percent"].values
    norm = mcolors.Normalize(vmin=atrophy_vals.min(), vmax=atrophy_vals.max())

    # Colormap: arancio per mossy_fibers, blu per purkinje
    cmap = cm.Oranges.reversed() if pop == "mossy_fibers" else cm.Blues.reversed()

    fig, ax = plt.subplots(figsize=(10, 5))

    for _, row in df_sorted.iterrows():
        y = np.array(row[f"{pop}_timeseries"])
        t = np.arange(len(y)) * dt
        atrophy = row["atrophy_percent"]
        color = cmap(norm(atrophy))
        ax.plot(t, y, color=color, alpha=alpha, linewidth=linewidth)

    ax.set_title(f"{pop.capitalize()} – Firing rate medio nel tempo", fontsize=12)
    ax.set_xlabel("Tempo (ms)")
    ax.set_ylabel("Firing rate (Hz)")
    ax.grid(True)

    total_time = len(df_sorted.iloc[0][f"{pop}_timeseries"]) * dt
    ax.set_xlim(750, total_time - 750)

    # Colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.01)
    cbar.set_label("Atrophy (%)", fontsize=10)

    output_path = f"./{base_dir}/{pop}_timeseries.png"
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=dpi)
    plt.show()

from collections import defaultdict
import neo
import pandas as pd
import numpy as np
import gc
from pathlib import Path

def extract_spike_df(path, downsample_target_rate_hz=None, duration_ms=None):
    """
    Estrae tutti gli spike da una cartella `results_*/basal_activity` e restituisce un DataFrame.
    Se `downsample_target_rate_hz` è specificato, effettua un sottocampionamento degli spike.
    """
    spikes = []

    for f in sorted(Path(path).glob("*.nio")):
        reader = neo.io.NixIO(filename=str(f), mode="ro")
        block = reader.read_block()
        for seg_idx, segment in enumerate(block.segments):
            for st_idx, st in enumerate(segment.spiketrains):
                times = st.times.rescale("ms").magnitude
                senders = st.annotations.get("senders", [f"unit{st_idx}"])
                if not isinstance(senders, (list, np.ndarray)):
                    senders = [senders]
                pop = st.annotations.get("device", "unknown").split("_")[0]
                for time, sid in zip(times, senders):
                    spikes.append({
                        "time_ms": time,
                        "sender_id": sid,
                        "pop": pop,
                        "file": f.name
                    })

        del reader, block
        gc.collect()

    df_spikes = pd.DataFrame(spikes)

    if downsample_target_rate_hz and not df_spikes.empty:
        if duration_ms is None:
            duration_ms = df_spikes["time_ms"].max()

        total_time_s = duration_ms / 1000.0
        max_spikes = downsample_target_rate_hz * total_time_s

        # Downsample: prendi una frazione casuale degli spike
        keep_fraction = min(1.0, max_spikes / len(df_spikes))
        df_spikes = df_spikes.sample(frac=keep_fraction, random_state=42).reset_index(drop=True)

        print(f"🔻 Downsampled: conservati {len(df_spikes)} spike su {len(spikes)} "
              f"({keep_fraction:.3f} = ~{downsample_target_rate_hz:.1f} Hz)")

    return df_spikes

def plot_classical_raster(df_spikes, df_row_with_stim, pop="mossy_fibers", window=None,
                          save_path=None, event_times=None,
                          title="Classical Raster Plot with Stimulus"):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    df_plot = df_spikes[df_spikes["pop"] == pop].copy()

    if window and event_times:
        spikes = []
        for t_ref in event_times:
            aligned = df_plot["time_ms"] - t_ref
            mask = (aligned >= window[0]) & (aligned <= window[1])
            spikes.append(df_plot.loc[mask].assign(time_aligned=aligned[mask]))
        df_plot = pd.concat(spikes)
        time_col = "time_aligned"
        stim_x = np.arange(len(df_row_with_stim["stim_profile"])) - event_times[0]
    else:
        time_col = "time_ms"
        stim_x = np.arange(len(df_row_with_stim["stim_profile"]))

    stim_y = np.array(df_row_with_stim["stim_profile"])

    # ⬅️ Estrai atrofia
    atrophy = df_row_with_stim.get("atrophy_percent", None)
    if atrophy is not None:
        title = f"{pop.capitalize()} – Raster aligned to event – Atrophy {atrophy:.1f}%"

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.scatter(df_plot[time_col], df_plot["sender_id"], s=1, color="black", label="Spikes")
    ax.set_xlabel("Time (ms)" if not window else "Time from event (ms)")
    ax.set_ylabel("Neuron ID")
    ax.set_title(title)

    # Onset line and stim overlay
    if window:
        ax.axvline(0, color="red", linestyle="--", label="Event onset")
        durations = df_row_with_stim.get("burst_durations_ms", [])
        if durations:
            stim_dur = np.mean(durations)
            ymin, ymax = ax.get_ylim()
            ax.fill_betweenx([ymin, ymax], 0, stim_dur, color="red", alpha=0.3, label="Stim duration")
        ax.legend()

    # Stim inset
    ax2 = ax.inset_axes([0, -0.15, 1, 0.1], transform=ax.transAxes)
    ax2.plot(stim_x, stim_y, color="orange")
    ax2.set_xlim(ax.get_xlim())
    ax2.axis("off")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"✅ Figure saved to: {save_path}")
    plt.show()

def average_psth(df, pop="mossy_fibers", event_times=None, window=(-200, 200), bin_size=5.0,
                 dt=0.1, n_min_events=10, save_path=None, condition_name="cond",
                 plot_raster=False, df_row_with_stim=None):
    import numpy as np
    import matplotlib.pyplot as plt

    if event_times is None or len(event_times) < n_min_events:
        print("⚠️ Non ci sono eventi sufficienti per la PSTH.")
        return

    bins = np.arange(window[0], window[1] + bin_size, bin_size)
    all_counts = []
    raster_data = []

    for idx, t_ref in enumerate(event_times):
        aligned = df[df["pop"] == pop]["time_ms"].values - t_ref
        counts, _ = np.histogram(aligned, bins=bins)
        all_counts.append(counts)

        if plot_raster:
            aligned_in_window = [t for t in aligned if window[0] <= t <= window[1]]
            for i, t in enumerate(aligned_in_window):
                if i % 10 == 0:
                    raster_data.append((idx, t))

    all_counts = np.array(all_counts)
    n_trials = len(event_times)
    n_neurons = df[df["pop"] == pop]["sender_id"].nunique()
    mean_rate = all_counts.mean(axis=0) / (bin_size*2 / 1000) / n_trials
    sem_rate = all_counts.std(axis=0) / np.sqrt(n_trials) / (bin_size*2 / 1000) / n_trials
    centers = bins[:-1] + bin_size / 2

    # 🎨 Colori: nero se il centro bin è < 0, bianco se ≥ 0
    colors = ["black" if c < 0 else "white" for c in centers]

    fig = plt.figure(figsize=(8, 6 if plot_raster else 4))
    fig.suptitle(pop)

    if plot_raster:
        ax1 = fig.add_subplot(211)
        ax1.scatter([t for _, t in raster_data], [i for i, _ in raster_data], s=2, color="black")
        ax1.axvline(0, color="red", linestyle="--")
        ax1.set_ylabel("Trial")
        atrophy = df_row_with_stim.get("atrophy_percent", "N/A")
        ax1.set_title(f"Atrophy {atrophy:.1f}% – Raster & PSTH ({pop})")
        # ➕ Stimulus duration as red shaded area also in raster
        if df_row_with_stim is not None:
            durations = df_row_with_stim.get("burst_durations_ms", [])
            if durations:
                stim_dur = np.mean(durations)
                ax1.fill_betweenx([0, len(event_times)], 0, stim_dur, color="red", alpha=0.3)

        ax2 = fig.add_subplot(212, sharex=ax1)
    else:
        ax2 = fig.add_subplot(111)
        atrophy = df_row_with_stim.get("atrophy_percent", "N/A")
        ax2.set_title(f"Atrophy {atrophy:.1f}% – PSTH ({pop}, n={n_trials} eventi)")

    # ➕ Colore per bin singolo
    ax2.bar(centers, mean_rate, width=bin_size, color=colors, edgecolor="black", alpha=0.9, label=" ")

    # ➕ Onset e durata stimolo
    ax2.axvline(0, color="red", linestyle="--", label="Onset")
    if df_row_with_stim is not None:
        durations = df_row_with_stim.get("burst_durations_ms", [])
        if durations:
            stim_dur = np.mean(durations)
            ax2.fill_betweenx([0, mean_rate.max()*1.1], 0, stim_dur, color="red", alpha=0.3, label="Stim duration")

    ax2.set_xlabel("Tempo (ms)")
    ax2.set_ylabel("Discharge rate (spikes/s)")
    ax2.legend()
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"✅ Figura salvata in: {save_path}")

    plt.show()

def plot_epoched_firing_rate(df, trace_column, pop_name="Population", stim_column=None,
                              onset_column="burst_onsets_ms", group_column=None, 
                              window_ms=(-200, 200), dt=0.1, cmap="viridis",
                              operator=np.mean,
                              sort_by_group=True, save_path=None):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    from pathlib import Path
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    win_pre = int(abs(window_ms[0]) / dt)
    win_post = int(window_ms[1] / dt)
    time_axis = np.arange(-win_pre, win_post) * dt

    if group_column:
        df = df.sort_values(group_column).reset_index(drop=True)
        group_vals = df[group_column].values
        norm = mcolors.Normalize(vmin=group_vals.min(), vmax=group_vals.max())
        cmap = cm.get_cmap(cmap)
    else:
        norm = None
        cmap = None

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    for _, row in df.iterrows():
        trace = np.array(row[trace_column])
        onsets = row[onset_column]
        stim = np.array(row[stim_column]) if stim_column else None
        group_val = row[group_column] if group_column else None

        epochs, stim_epochs = [], []
        for onset in onsets:
            center = int(onset / dt)
            start = center - win_pre
            end = center + win_post
            if start >= 0 and end <= len(trace):
                epochs.append(trace[start:end])
                if stim_column:
                    stim_epochs.append(stim[start:end])

        if epochs:
            epochs = np.array(epochs)
            mean_trace = np.mean(epochs, axis=0)
            lower = np.percentile(epochs, 25, axis=0)
            upper = np.percentile(epochs, 75, axis=0)
            
            color = cmap(norm(group_val)) if cmap else "black"

            ax.plot(time_axis, mean_trace, linewidth=2, color=color, alpha=0.9)
            ax.fill_between(time_axis, lower, upper, color=color, alpha=0.3, label=None)

            if stim_column and stim_epochs:
                stim_epochs = np.array(stim_epochs)
                stim_mean = operator(stim_epochs, axis=0)
                ax.plot(time_axis, stim_mean, color="gray", alpha=0.3, linewidth=1.5)


    # ➕ Onset line and stimulus duration overlay
    ax.axvline(0, color='red', linestyle='--', label="Stim onset")
    if "burst_durations_ms" in df.columns:
        all_durations = [d for durations in df["burst_durations_ms"] for d in durations if d > 0]
        if all_durations:
            stim_duration = operator(all_durations)
            ax.fill_betweenx(
                [ax.get_ylim()[0], ax.get_ylim()[1]],
                0, stim_duration,
                color='red', alpha=0.3, label="Stim duration"
            )

    ax.set_title(f"{pop_name} – Epoched Firing Rate Aligned to Burst Onset")
    ax.set_xlabel("Time from Burst Onset (ms)")
    ax.set_ylabel("Firing Rate (Hz)")

    # Optional colorbar
    if group_column:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.1)
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label("Atrophy (%)")

    ax.legend(loc="upper right", fontsize=9)
    plt.tight_layout()
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300)
        print(f"✅ Saved figure to: {save_path}")
    plt.show()

import numpy as np
import pandas as pd
from collections import defaultdict
from pathlib import Path
import yaml

def generate_stim_profile(devices: dict, sim_duration: float, dt: float = 0.1) -> np.ndarray:
    """
    Genera un profilo dello stimolo temporale (stim_profile) da una lista di devices Poisson nel file YAML BSB.

    Args:
        devices (dict): Dizionario dei devices dalla sezione YAML.
        sim_duration (float): Durata della simulazione in ms.
        dt (float): Risoluzione temporale in ms.

    Returns:
        np.ndarray: array 1D con valori del rate a ogni time step.
    """
    stim_profile = np.zeros(int(sim_duration // dt))
    for name, props in devices.items():
        if props.get("device") != "poisson_generator":
            continue
        start = float(props.get("start", 0.0))
        stop = float(props.get("stop", sim_duration))
        rate = float(props.get("rate", 0.0))
        start_idx = int(start // dt)
        stop_idx = int(stop // dt)
        stim_profile[start_idx:stop_idx] += rate
    return stim_profile

from collections import defaultdict
import neo
import pandas as pd
import numpy as np
import gc
from pathlib import Path
import matplotlib.pyplot as plt

def extract_spike_df(path, downsample_target_rate_hz=None, duration_ms=None):
    spikes = []
    for f in sorted(Path(path).glob("*.nio")):
        reader = neo.io.NixIO(filename=str(f), mode="ro")
        block = reader.read_block()
        for seg_idx, segment in enumerate(block.segments):
            for st_idx, st in enumerate(segment.spiketrains):
                times = st.times.rescale("ms").magnitude
                senders = st.annotations.get("senders", [f"unit{st_idx}"])
                if not isinstance(senders, (list, np.ndarray)):
                    senders = [senders]
                #pop = st.annotations.get("device", "unknown").split("_")[0]
                pop = st.annotations.get("device", "unknown").removesuffix("_record")

                for time, sid in zip(times, senders):
                    spikes.append({
                        "time_ms": time,
                        "sender_id": sid,
                        "pop": pop,
                        "file": f.name
                    })
        del reader, block
        gc.collect()

    df_spikes = pd.DataFrame(spikes)
    if downsample_target_rate_hz and not df_spikes.empty:
        if duration_ms is None:
            duration_ms = df_spikes["time_ms"].max()
        total_time_s = duration_ms / 1000.0
        max_spikes = downsample_target_rate_hz * total_time_s
        keep_fraction = min(1.0, max_spikes / len(df_spikes))
        df_spikes = df_spikes.sample(frac=keep_fraction, random_state=42).reset_index(drop=True)
        print(f"🔻 Downsampled: conservati {len(df_spikes)} spike su {len(spikes)} ({keep_fraction:.3f} = ~{downsample_target_rate_hz:.1f} Hz)")
    return df_spikes


def plot_raster(df_spikes, pop="mossy_fibers", time_window=None, max_units=100):
    df = df_spikes[df_spikes["pop"] == pop]
    unique_units = df["sender_id"].unique()
    if len(unique_units) > max_units:
        unique_units = np.random.choice(unique_units, size=max_units, replace=False)
    df = df[df["sender_id"].isin(unique_units)]
    df = df.sort_values("sender_id")
    unit_map = {unit: idx for idx, unit in enumerate(df["sender_id"].unique())}
    df["unit_idx"] = df["sender_id"].map(unit_map)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(df["time_ms"], df["unit_idx"], s=2, color="black")
    if time_window:
        ax.set_xlim(time_window)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Unit index")
    ax.set_title(f"Raster Plot – {pop} neurons")
    plt.tight_layout()
    plt.show()

def plot_cosine_similarity_vs_atrophy(df_all, column_name="cosine_similarity", save_path=None):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import spearmanr

    plt.figure(figsize=(10, 6))
    ax = plt.gca()

    sns.scatterplot(
        data=df_all,
        x="Atrophy Factor (%)",
        y=column_name,
        s=100,
        color="purple",
        edgecolor="black",
        label=column_name.replace("_", " ").capitalize()
    )

    sns.regplot(
        data=df_all,
        x="Atrophy Factor (%)",
        y=column_name,
        scatter=False,
        color="purple",
        ax=ax
    )

    # Calcolo Spearman
    rho, pval = spearmanr(df_all["Atrophy Factor (%)"], df_all[column_name])
    ax.legend(title=f"ρ = {rho:.2f}, p = {pval:.2e}", loc="upper right")

    # Titoli e assi
    ax.set_title(f"{column_name.replace('_', ' ').capitalize()} vs Atrophy", fontsize=14)
    ax.set_xlabel("Atrophy Factor (%)", fontsize=12)
    ax.set_ylabel("Cosine Similarity", fontsize=12)
    ax.set_ylim(0, 1)  # fissa y tra 0 e 1
    ax.grid(True, linestyle='--', alpha=0.4)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"✅ Cosine similarity plot salvato in: {save_path}")
    plt.show()

def compute_firing_rate(spike_times, duration, dt=1.0, sigma=100):
    time = np.arange(0, duration, dt)
    spike_train = np.zeros_like(time)
    indices = (np.array(spike_times) / dt).astype(int)
    indices = indices[indices < len(spike_train)]
    spike_train[indices] = 1
    kernel_size = max(int(2 * sigma / dt), 3)
    kernel = triang(kernel_size)
    kernel /= kernel.sum()
    rate = convolve(spike_train, kernel, mode='same') * (1000.0 / dt)
    return time, rate

def plot_stim_vs_mossy_fibers_rate(df_all, stim_column="stim_profile", mossy_fibers_column="mossy_fibers_mean_rate", save_path=None):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import spearmanr
    import numpy as np

    # Calcolo del rate medio dello stimolo per ogni riga
    df_all["stim_mean_rate"] = df_all[stim_column].apply(lambda x: np.mean(x) if isinstance(x, list) else np.nan)

    plt.figure(figsize=(7, 5))
    ax = plt.gca()

    sns.scatterplot(
        data=df_all,
        x="stim_mean_rate",
        y=mossy_fibers_column,
        s=80,
        color="teal",
        edgecolor="black",
        label="mossy_fibers vs Stim"
    )

    sns.regplot(
        data=df_all,
        x="stim_mean_rate",
        y=mossy_fibers_column,
        scatter=False,
        color="teal",
        ax=ax
    )

    # Correlazione Spearman
    rho, pval = spearmanr(df_all["stim_mean_rate"], df_all[mossy_fibers_column])
    ax.legend(title=f"ρ = {rho:.2f}, p = {pval:.2e}", loc="upper left")

    ax.set_title("Stimulus Mean Rate vs mossy_fibers Firing Rate")
    ax.set_xlabel("Stimulus Mean Rate (Hz)")
    ax.set_ylabel("mossy_fibers Mean Firing Rate (Hz)")
    ax.grid(True)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"✅ Stimulus vs mossy_fibers rate plot salvato in: {save_path}")
    plt.show()

def cosine_similarity(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))

import pandas as pd
import numpy as np
from collections import defaultdict
from pathlib import Path
from scipy.stats import sem

from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import sem

from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import sem

from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import sem
import yaml  # 🆕 necessario per salvataggio YAML

from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import sem
from scipy.signal import welch
import yaml
from fooof import FOOOF

def plotFRglobal(data_path: Path, df_summary):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from scipy.stats import spearmanr, pearsonr
    
    # Mappa dei colori (non usata perché tutti i plot sono in nero)
    cell_color_map = {
    "mossy_fibers": "purple",
    "granule": "gold",
    "golgi": "green",
    "basket": "blue",
    "stellate": "pink",
    "purkinje": "brown",
    }
    
    # Ordine desiderato
    ordered_cells = ["mossy_fibers", "granule", "golgi", "glomerulus", "basket", "stellate", "purkinje"]
    # add dcn-io
    cell_types = [c for c in ordered_cells if c in df_summary["pop"].unique()]
    #cell_types = df_summary["pop"].unique()
    # Categorie
    excitatory = ["granule"]
    inhibitory = ["golgi", "basket", "stellate", "purkinje"]
    
    # Colonna coerente per asse X
    #df_summary["Atrophy Level (%)"] = 100 - df_summary["k"] * 100
    
    # Funzione per aggregare per sottogruppo
    def compute_summary(df, subset):
        filtered = df[df["pop"].isin(subset)]
        grp = filtered.groupby("Atrophy Level (%)").agg(
            fr_mean=("firing_rate_mean", "mean"),
            fr_sem=("firing_rate_sem", "mean")
        ).reset_index()
        return grp
    
    df_all = compute_summary(df_summary, cell_types)
    df_exc = compute_summary(df_summary, excitatory)
    df_inh = compute_summary(df_summary, inhibitory)
    df_ratio = pd.merge(df_exc, df_inh, on="Atrophy Level (%)", suffixes=("_exc", "_inh"))
    df_ratio["ratio_EI"] = df_ratio["fr_mean_exc"] / df_ratio["fr_mean_inh"]
    
    # Layout: 4 metriche principali + 1 per tipo cellulare
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(18, 12))
    axes = axes.flatten()

    """    
    # Funzione di plot
    def plot_with_corr(ax, df, x_col, y_col, yerr_col=None, label=""):
        x = df[x_col].values
        y = df[y_col].values
        yerr = df[yerr_col].values if yerr_col else None
        ax.errorbar(x, y, yerr=yerr, fmt='o', color="black", ecolor="black", capsize=4, alpha=1)
        sns.regplot(x=x, y=y, ax=ax, scatter=False, order=2,
                    line_kws={"color": "black", "linewidth": 0, 'alpha': 0.15}, ci=None)
        rho, p_rho = spearmanr(x, y)
        r, p_r = pearsonr(x, y)
        ax.set_title(f"{label}\nρ={rho:.2f}, p={p_rho:.1e} | r={r:.2f}, p={p_r:.1e}")
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col.replace("_", " ").capitalize())
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if 0 in x:
            baseline_y = y[x == 0][0]
            ax.axhline(y=baseline_y, color="black", linestyle="--", linewidth=1.5, alpha=0.5)
    """
    from scipy.stats import pearsonr, spearmanr
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    def plot_with_corr(ax, df, x_col, y_col, yerr_col=None, label=""):
        x = df[x_col].values
        y = df[y_col].values
        yerr = df[yerr_col].values if yerr_col and yerr_col in df.columns else None
    
        # Plot base con barre di errore
        ax.errorbar(x, y, yerr=yerr, fmt='o', color="black", ecolor="black", capsize=4, alpha=1)
    
        # Regressione di supporto (linea trasparente)
        sns.regplot(x=x, y=y, ax=ax, scatter=False, order=2,
                    line_kws={"color": "black", "linewidth": 0, 'alpha': 0.15}, ci=None)
    
        # Calcolo delle correlazioni solo se ci sono almeno due punti
        if len(x) >= 2 and len(y) >= 2:
            rho, p_rho = spearmanr(x, y)
            r, p_r = pearsonr(x, y)
            title_corr = f"\nρ={rho:.2f}, p={p_rho:.1e} | r={r:.2f}, p={p_r:.1e}"
        else:
            title_corr = "\n⚠️ not enough data for correlation"
    
        # Titolo e labels
        ax.set_title(f"{label}{title_corr}")
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col.replace("_", " ").capitalize())
    
        # Styling
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
        # Linea di riferimento per baseline (x == 0)
        if 0 in x:
            baseline_y = y[x == 0][0]
            ax.axhline(y=baseline_y, color="black", linestyle="--", linewidth=1.5, alpha=0.5)

    
    # Plot aggregati
    plot_with_corr(axes[0], df_all, "Atrophy Level (%)", "fr_mean", "fr_sem", "All cells")
    #plot_with_corr(axes[1], df_exc, "Atrophy Level (%)", "fr_mean", "fr_sem", "Excitatory")
    #plot_with_corr(axes[2], df_inh, "Atrophy Level (%)", "fr_mean", "fr_sem", "Inhibitory")
    plot_with_corr(axes[1], df_ratio, "Atrophy Level (%)", "ratio_EI", label="E/I ratio")
    
    # Indici dei plot effettivamente usati
    used_slots = [0, 1]  # all cells + E/I ratio
    
    # Plot per tipo cellulare
    for i, cell in enumerate(cell_types):
        cell_df = df_summary[df_summary["pop"] == cell]
        slot_index = 4 + i
        used_slots.append(slot_index)
        plot_with_corr(axes[slot_index], cell_df, "Atrophy Level (%)", "firing_rate_mean", "firing_rate_sem", cell.capitalize())

    # Rendi invisibili gli altri pannelli
    for j, ax in enumerate(axes):
        if j not in used_slots:
            ax.set_visible(False)

    # Titolo principale
    fig.suptitle("Cerebellar Network Activity", fontsize=18, fontweight='bold')

    fig.tight_layout()
    plt.subplots_adjust(hspace=0.55, wspace=0.3, top=0.90)

    # Salvataggio
    out_dir = Path(data_path) / "B_collective_results/2.functions"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "fr_summary_plot.png"    
    fig.savefig(out_path, dpi=300)
    plt.close()
    print(f"💾 Plot salvato in: {out_path}")
    
def plotFOOOFglobal(data_path: Path,df_summary):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from scipy.stats import spearmanr, pearsonr
    from pathlib import Path

    # Ordine desiderato
    ordered_cells = ["mossy_fibers", "granule", "golgi", "basket", "stellate", "purkinje"]
    cell_types = [c for c in ordered_cells if c in df_summary["pop"].unique()]

    # Estrai parametri FOOOF
    def extract_fooof_features(df):
        return pd.DataFrame({
            "Atrophy Level (%)": df["Atrophy Level (%)"],
            "pop": df["pop"],
            "offset": df.get("FOOOF_aperiodic", np.nan).apply(lambda x: x[0] if isinstance(x, list) else np.nan),
            "slope": df.get("FOOOF_aperiodic", np.nan).apply(lambda x: x[1] if isinstance(x, list) else np.nan),
            "n_peaks": df.get("FOOOF_n_peaks", np.nan)
        })

    df_fooof = extract_fooof_features(df_summary)

    features = ["offset", "slope"]
    n_features = len(features)
    fig, axes = plt.subplots(nrows=n_features, ncols=len(cell_types)+1, figsize=(22, 8), sharey='row')
    axes = np.array(axes)

    def plot_feature(ax, df, feature, label):
        df_clean = df[[feature, "Atrophy Level (%)"]].dropna()
        if df_clean.empty:
            ax.set_title(f"{label}\n(no data)")
            return

        x = df_clean["Atrophy Level (%)"]
        y = df_clean[feature]

        sns.scatterplot(x=x, y=y, ax=ax, color="black")
        sns.regplot(x=x, y=y, ax=ax, scatter=False, order=1,
                    line_kws={"color": "black", "alpha": 0}, ci=None)

        # Linea orizzontale baseline (AF = 0)
        if 0 in x.values:
            y0 = y[x == 0].values
            if len(y0) > 0 and not np.isnan(y0[0]):
                ax.axhline(y=y0[0], linestyle="--", linewidth=1.5, alpha=0.5, c='k')

        try:
            rho, p_rho = spearmanr(x, y)
            r, p_r = pearsonr(x, y)
            ax.set_title(f"{label}\nρ={rho:.2f}, p={p_rho:.1e}\nr={r:.2f}, p={p_r:.1e}")
        except Exception:
            ax.set_title(f"{label}\n(no corr)")
        ax.set_xlabel("Atrophy Level (%)")
        ax.set_ylabel(feature.capitalize() if label == "All" else "")

    for i, feature in enumerate(features):
        # Global
        plot_feature(axes[i, 0], df_fooof, feature, "All")
        # Per popolazione
        for j, cell in enumerate(cell_types):
            df_cell = df_fooof[df_fooof["pop"] == cell]
            plot_feature(axes[i, j+1], df_cell, feature, cell.capitalize())

    # Titolo globale
    fig.suptitle("Cerebellar Network Activity — FOOOF Features vs Atrophy Level",
                 fontsize=20, fontweight='bold', y=1.05)
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.4, top=0.88)

    # Salvataggio
    out_dir = Path(data_path) / "B_collective_results/2.functions"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "fooof_summary_plot.png"    
    fig.savefig(out_path, dpi=300)
    plt.close()
    print(f"💾 Plot salvato in: {out_path}")


def loadData(data_path: Path, sigma: float = 10.0, downsample_target_rate_hz: float = None) -> pd.DataFrame:
    dt = 0.1
    sigma = sigma
    freq_range = (1, 100)
    fs = 1000 / dt
    nperseg = 2048

    summary_rows = []
    spike_dfs = {}

    for path in data_path.rglob("*_deg_*/3.function"):
        match = [part for part in path.parts if "_deg_" in part]
        if not match:
            print(f"⚠️ No k found in: {path}")
            continue
        folder_name = match[0]
        try:
            k_str = folder_name.split("_")[2]
            k = float(k_str)
        except Exception as e:
            print(f"⚠️ Errore parsing k da '{folder_name}': {e}")
            continue

        print(f"📁 Analisi per k={k:.3f}")
        df = extract_spike_df(path, downsample_target_rate_hz=downsample_target_rate_hz)
        spike_dfs[k] = df.copy()
        if df.empty:
            print("⚠️ No data after filtering")
            continue

        duration = df["time_ms"].max()
        grouped = df.groupby(["pop", "sender_id"])
        fr_data = defaultdict(list)

        for (pop, sid), g in grouped:
            spike_times = g["time_ms"].values
            t, r = compute_firing_rate(spike_times, duration, dt=dt, sigma=sigma)
            fr_data["pop"].append(pop)
            fr_data["sender_id"].append(sid)
            fr_data["firing_rate"].append(np.mean(r))
            if len(spike_times) >= 2:
                isi_values = np.diff(np.sort(spike_times))
                fr_data["isi"].append(np.mean(isi_values))
            else:
                fr_data["isi"].append(np.nan)

        df_metrics = pd.DataFrame(fr_data)
        grouped_pop = df_metrics.groupby("pop")

        stats_dict = {}

        for pop_name, g in grouped_pop:
            stats_row = {
                "k": k,
                "Atrophy Level (%)": 100*(1-k), 
                "firing_rate_mean": float(g["firing_rate"].mean()),
                "firing_rate_std": float(g["firing_rate"].std()),
                "firing_rate_sem": float(sem(g["firing_rate"], nan_policy="omit")) if len(g) > 1 else 0.0,
                "isi": float(g["isi"].mean()),
                "n_cells": int(len(g))
            }
        
            try:
                sel = df["pop"] == pop_name
                df_pop = df[sel]
                if not df_pop.empty:
                    spike_mat = defaultdict(list)
                    for sid, g_unit in df_pop.groupby("sender_id"):
                        spikes = g_unit["time_ms"].values
                        t, r = compute_firing_rate(spikes, duration, dt=dt, sigma=sigma)
                        spike_mat["rates"].append(r)
                """
                if spike_mat["rates"]:
                        rate_array = np.vstack(spike_mat["rates"])
                        mean_rate = np.mean(rate_array, axis=0)
                        fvals, psd = welch(mean_rate, fs=fs, nperseg=nperseg)
                        resolution = fs / nperseg
                        min_width = max(2 * resolution, 1.0)  # raccomandato >= 2x risoluzione
                        fm = FOOOF(peak_width_limits=[min_width, 12], max_n_peaks=6)
                        fm.fit(fvals, psd, freq_range)
                        
                        # 🔽 Salva il report PDF FOOOF
                        report_name = f"fooof_report_k{k:.3f}_{pop_name}.pdf"
                        report_path = path / report_name
                        fm.save_report(str(report_path), save_format='pdf')
                    
                        # 🔽 Salva i parametri nel dizionario
                        stats_row.update({
                            "FOOOF_aperiodic": list(map(float, fm.aperiodic_params_)),
                            "FOOOF_n_peaks": int(len(fm.gaussian_params_)),
                            "FOOOF_gaussian_params": [list(map(float, p)) for p in fm.gaussian_params_]
                        })
                """
            except Exception as e:
                print(f"⚠️ Errore FOOOF per {pop_name} (k={k:.3f}): {e}")
            
            summary_rows.append({"k": k, "pop": pop_name, **stats_row})
            stats_dict[pop_name] = stats_row

        # Salva YAML locale nella subfolder
        yaml_path = path / "stats.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(stats_dict, f, sort_keys=False)

    df_summary = pd.DataFrame(summary_rows)
    out_dir = data_path / "B_collective_results/2.functions"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "df_summary.pkl"
    df_summary.to_pickle(out_path)
    print(f"💾 Salvato df_summary in: {out_path}")

    return df_summary, spike_dfs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, required=True, help='Percorso alla cartella dati')
    parser.add_argument("--sigma", type=float, default=100.0, help="Sigma for kernel smoothing in firing rate estimation")
    args = parser.parse_args()
    base_dir = Path(args.path)
    #process_reconstruction_yamls(base_dir)
    df_summary, spike_dfs = loadData(base_dir, sigma=args.sigma)
    plotFRglobal(base_dir, df_summary)
    #plotFOOOFglobal(base_dir, df_summary)
    #df_summary, spike_dfs = loadData(base_dir, mainEventKey="burst_burst", condition_name="saccades", use_top_quartile=True, selection_mode="mean")

if __name__ == "__main__":
    main()
