import numpy as np
import pandas as pd
from typing import Union
import shutil
from pathlib import Path
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
import matplotlib.pyplot as plt
import os
import yaml
import matplotlib.pyplot as plt
from pathlib import Path

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


def move_png_to_sigma_folder(path: Path, sigma: float):
    dest_dir = path / "B_collective_results" / "2.functions" / f"sigma_{sigma}"
    dest_dir.mkdir(parents=True, exist_ok=True)

    png_files = list(path.glob("*.png"))
    for f in png_files:
        shutil.move(str(f), dest_dir / f.name)

    print(f"📦 Spostati {len(png_files)} PNG in: {dest_dir}")



import matplotlib.pyplot as plt
import numpy as np



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

def plot_tf(df_all, base_dir):

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
        if len(group) >= 2:
            r1, p1 = pearsonr(group["Atrophy Factor (%)"], group["te_mossy_fibers_to_purk"])
            label1 = f"{cond} MF→PC (r={r1:.2f}, p={p1:.1e})"
        else:
            label1 = f"{cond}"
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
        if len(group) >= 2:
            r2, p2 = pearsonr(group["Atrophy Factor (%)"], group["te_purk_to_mossy_fibers"])
            label2 = f"{cond} PC→MF (r={r2:.2f}, p={p2:.1e})"
        else:
            label2 = f"{cond}"
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
    if len(df_all) >= 2:
        r_diff, p_diff = pearsonr(df_all["Atrophy Factor (%)"], df_all["ΔTE"])
        label=f"ΔTE (r={r_diff:.2f}, p={p_diff:.1e})"
    else:
        label=f"ΔTE"
    sns.scatterplot(
        data=df_all,
        x="Atrophy Factor (%)",
        y="ΔTE",
        ax=ax2,
        color=color_diff,
        marker="D",
        label=label,
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

#def plot_cross_correlation(df_all):
def plot_cross_correlation(df_all, base_dir: Union[str, Path]):

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
        if len(group_filtered) >= 2:
            r1, p1 = pearsonr(group_filtered["Atrophy Factor (%)"], group_filtered["cross_corr"])
            label1 = f"{cond} – Corr (r={r1:.2f}, p={p1:.1e})"
        else:
            label1= f"{cond}"
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
        if len(group_filtered) >= 2:
            r2, p2 = pearsonr(group_filtered["Atrophy Factor (%)"], group_filtered["optimal_lag_ms"])
            label2 = f"{cond} – Lag (r={r2:.2f}, p={p2:.1e})"
        else:
            label2= f"{cond}"
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
    plt.close()

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
    plt.close()

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
    plt.close()

def plot_time_series(df_all, pop="mossy_fibers", base_dir=Union[str, Path], dt=0.1, alpha=0.15, linewidth=1.0, dpi=300):
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
    plt.close()

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
    plt.close()

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

    plt.close()

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
    plt.close()

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
    plt.close()

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
    plt.close()

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
    plt.close()

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import random

def exp_decay(t, A, tau, C):
    return A * np.exp(-t / tau) + C

def plot_tau_fit_examples(trace, onsets, dt=0.1, window=(0, 300), title="", savepath="tau_fits_example.png", color='blue'):
    """
    Plotta esempi di fit esponenziale per un segnale 1D intorno a stimoli.

    Args:
        trace (np.ndarray): segnale 1D.
        onsets (list of float): lista di onset (in ms).
        dt (float): risoluzione temporale in ms.
        window (tuple): finestra post-onset (in ms).
        title (str): titolo del plot.
        savepath (str): path del file .png da salvare.
        color (str): colore del segnale.
    """
    win_start, win_end = int(window[0] / dt), int(window[1] / dt)
    selected_onsets = random.sample(onsets, min(3, len(onsets)))

    plt.figure(figsize=(10, 4))
    for i, onset in enumerate(selected_onsets):
        center = int(onset / dt)
        start = center + win_start
        end = center + win_end
        if start < 0 or end > len(trace):
            continue
        segment = trace[start:end]
        t = np.arange(len(segment)) * dt
        y = segment - np.min(segment)
        try:
            popt, _ = curve_fit(exp_decay, t, y, p0=(np.max(y), 100.0, 0.0), maxfev=2000)
            fit = exp_decay(t, *popt)
            tau = popt[1]
            plt.plot(t, y, label=f"Epoch {i+1} τ={tau:.1f} ms", color=color, alpha=0.6)
            plt.plot(t, fit, linestyle='--', color=color, alpha=0.6)
        except Exception:
            continue

    plt.xlabel("Time (ms)")
    plt.ylabel("Signal")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import random

def exp_decay(t, A, tau, C):
    return A * np.exp(-t / tau) + C

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import random

def exp_decay(t, A, tau, C):
    return A * np.exp(-t / tau) + C

def plot_tau_fit_examples_v2(trace, onsets, dt=0.1, window=(0, 100), title="", savepath="tau_fits_example.png", color='blue'):
    """
    Plotta esempi di fit esponenziale a partire dal picco post-stimolo,
    calcola la media e std di τ solo sul decadimento.

    Args:
        trace (np.ndarray): segnale 1D.
        onsets (list of float): lista di onset (in ms).
        dt (float): risoluzione temporale in ms.
        window (tuple): finestra post-onset (in ms).
        title (str): titolo base del plot.
        savepath (str): path del file .png da salvare.
        color (str): colore del segnale.
    """
    win_start, win_end = int(window[0] / dt), int(window[1] / dt)
    selected_onsets = random.sample(onsets, len(onsets))
    tau_values = []

    plt.figure(figsize=(10, 4))
    for i, onset in enumerate(selected_onsets):
        center = int(onset / dt)
        start = center + win_start
        end = center + win_end
        if start < 0 or end > len(trace):
            continue
        segment = trace[start:end]
        t_segment = np.arange(len(segment)) * dt

        # Trova indice del picco
        peak_idx = np.argmax(segment)
        if peak_idx >= len(segment) - 5:  # troppo vicino alla fine
            continue

        # Seleziona solo il tratto in decadimento
        decay = segment[peak_idx:]
        t_decay = np.arange(len(decay)) * dt
        decay = decay - np.min(decay)  # normalizzazione a zero

        try:
            popt, _ = curve_fit(exp_decay, t_decay, decay, p0=(np.max(decay), 50.0, 0.0), maxfev=2000)
            fit = exp_decay(t_decay, *popt)
            tau = popt[1]
            tau_values.append(tau)

            # Plot completo
            plt.plot(t_segment, segment - np.min(segment), alpha=0.1, color=color)
            plt.plot(t_decay + t_segment[peak_idx], decay, label=f"Epoch {i+1} τ={tau:.1f} ms", color=color, alpha=0.8)
            plt.plot(t_decay + t_segment[peak_idx], fit, linestyle='--', color=color, alpha=0.8)
        except Exception:
            continue

    if tau_values:
        tau_mean = np.mean(tau_values)
        tau_std = np.std(tau_values)
        title += f" | τ = {tau_mean:.1f} ± {tau_std:.1f} ms"

    plt.xlabel("Time (ms)")
    plt.ylabel("Signal (norm.)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def exp_decay(t, A, tau, C):
    return A * np.exp(-t / tau) + C

def exp_rise(t, A, tau, C):
    return A * (1 - np.exp(-t / tau)) + C

def plot_tau_fits_mean_rise_decay(trace, onsets, dt=0.1, window=(-50, 150), title="", savepath="tau_fit_mean.png", color='blue'):
    """
    Plotta la media di epoche allineate all'onset con fit esponenziale di salita e decadimento.
    
    Args:
        trace (np.ndarray): segnale 1D.
        onsets (list): lista di onset (in ms).
        dt (float): tempo tra campioni, in ms.
        window (tuple): finestra attorno all'onset, in ms.
        title (str): titolo base del grafico.
        savepath (str): path del file .png.
        color (str): colore principale.
    """
    win_start, win_end = int(window[0] / dt), int(window[1] / dt)
    epoch_len = win_end - win_start
    epochs = []

    for onset in onsets:
        center = int(onset / dt)
        start = center + win_start
        end = center + win_end
        if start < 0 or end > len(trace):
            continue
        epoch = trace[start:end]
        epochs.append(epoch)

    if not epochs:
        print("⚠️ Nessuna epoca valida trovata.")
        return

    epochs = np.array(epochs)
    t = np.arange(epoch_len) * dt + window[0]
    mean_epoch = epochs.mean(axis=0)
    std_epoch = epochs.std(axis=0)
    norm_epoch = mean_epoch - np.min(mean_epoch)  # normalizza base

    # Trova picco
    peak_idx = np.argmax(norm_epoch)
    t_peak = t[peak_idx]

    # Fit di salita (pre-picco)
    t_rise = t[:peak_idx] - t[0]
    y_rise = norm_epoch[:peak_idx]
    try:
        popt_rise, _ = curve_fit(exp_rise, t_rise, y_rise, p0=(np.max(y_rise), 30.0, 0.0), maxfev=2000)
        fit_rise = exp_rise(t_rise, *popt_rise)
        tau_rise = popt_rise[1]
    except Exception:
        fit_rise = None
        tau_rise = np.nan

    # Fit di decadimento (post-picco)
    t_decay = t[peak_idx:] - t[peak_idx]
    y_decay = norm_epoch[peak_idx:]
    try:
        popt_decay, _ = curve_fit(exp_decay, t_decay, y_decay, p0=(np.max(y_decay), 50.0, 0.0), maxfev=2000)
        fit_decay = exp_decay(t_decay, *popt_decay)
        tau_decay = popt_decay[1]
    except Exception:
        fit_decay = None
        tau_decay = np.nan

    # Plot
    plt.figure(figsize=(10, 4))
    plt.fill_between(t, mean_epoch - std_epoch, mean_epoch + std_epoch, color=color, alpha=0.3, label="±1 std")
    plt.plot(t, mean_epoch, color=color, label="mean")

    if fit_rise is not None:
        plt.plot(t[:peak_idx], fit_rise + np.min(mean_epoch), linestyle='--', color='black',
                 label=f"Rise fit τ={tau_rise:.1f} ms")
    if fit_decay is not None:
        plt.plot(t[peak_idx:], fit_decay + np.min(mean_epoch), linestyle='--', color='gray',
                 label=f"Decay fit τ={tau_decay:.1f} ms")

    title += f" | τ_rise = {tau_rise:.1f} ms | τ_decay = {tau_decay:.1f} ms"
    plt.title(title)
    plt.xlabel("Time (ms)")
    plt.ylabel("Signal")
    plt.legend()
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()



def loadData(data_path: Path, mainEventKey: str, 
             stimName: str, 
             condition_name: str, use_top_quartile: bool = False, THR=50,
             selection_mode: str = "percentile", selection_percentile: float = 0.99, method: str = "mean",
             sigma: float = 10.0, downsample_target_rate_hz: float = None) -> pd.DataFrame:
    dt = 0.1
    #sigma = 10 #50 #10#25#15#5#15#dt*250 #100
    freq_range = (1, 100)
    trim_start = 0 #0.10
    trim_end = 0 #0.10
    base_dir = str(data_path)
    
    def cosine_similarity(a, b):
        return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    
    summary_rows = []
    spike_dfs = {}  # oppure assicurati che sia già stato eseguito

    for path in data_path.rglob("*_deg_*/3.function"):
        match = re.search(r"deg_(\d+\.\d+)", str(path))
        if not match:
            print(f"⚠️ No k found in: {path}")
            continue
        k = float(match.group(1))
        """
        match = [part for part in path.parts if part.startswith("*_deg_")]
        if not match:
            print(f"⚠️ No k found in: {path}")
            continue
        k = float(match[0].replace("*_deg_", ""))
        """
        atrophy = 100 * (1 - k)
        print(f"📁 {condition_name} – Analisi per k={k:.3f}")
    
        burst_onsets_ms, burst_durations_ms = [], []
        stim_profile = []
        yaml_file = path.parent / "2.stimulus" / f"{stimName}_{k:.3f}.yaml"
        if yaml_file.exists():
            import yaml
            with open(yaml_file, "r") as f:
                yaml_data = yaml.safe_load(f)
            #devices = yaml_data.get("simulations", {}).get("basal_activity", {}).get("devices", {})
            devices = yaml_data.get("simulations", {}).get(stimName, {}).get("devices", {})
            
            burst_rates_hz = []
            for name, props in devices.items():
                if name.startswith(mainEventKey) and props.get("device") == "poisson_generator":
                    start = float(props["start"])
                    stop = float(props["stop"])
                    rate = float(props.get("rate", 0.0))
                    burst_onsets_ms.append(start)
                    burst_durations_ms.append(stop - start)
                    burst_rates_hz.append(rate)

            duration = max([float(p.get("stop", 0.0)) for p in devices.values() if p.get("device") == "poisson_generator"] + [0.0])
            n_timepoints = int(duration // dt)
            stim_profile_arr = np.zeros(n_timepoints)
            for props in devices.values():
                if props.get("device") == "poisson_generator":
                    rate = float(props.get("rate", 0.0))
                    start_idx = int(float(props["start"]) // dt)
                    stop_idx = int(float(props.get("stop", duration)) // dt)
                    stim_profile_arr[start_idx:stop_idx] += rate
            stim_profile = stim_profile_arr.tolist()

        df = extract_spike_df(path, downsample_target_rate_hz=downsample_target_rate_hz)
        df = df[df['pop'].isin(['mossy_fibers', 'purkinje'])]
        spike_dfs[k] = df.copy()

        if df.empty:
            print("⚠️ No data after filtering")
            continue

        duration = df["time_ms"].max()
        grouped = df.groupby(["pop", "sender_id"])

        data = defaultdict(list)
        for (pop, sid), g in grouped:
            t, r = compute_firing_rate(g["time_ms"].values, duration, dt=dt, sigma=sigma)
            data["pop"].append(pop)
            data["sender_id"].append(sid)
            data["rate"].append(r)
            data["time"].append(t)

        def select_quartile(data, pop_name, mode="percentile", percentile=0.99):
            scores = []
            for i, r in enumerate(data["rate"]):
                if data["pop"][i] != pop_name:
                    continue
                if mode == "mean":
                    score = np.mean(r)
                elif mode == "max":
                    score = np.max(r)
                elif mode == "auc":
                    score = np.sum(r)
                elif mode == "percentile":
                    score = np.percentile(r, percentile * 100)
                elif mode == "timepoint":
                    score = r[-1]
                else:
                    raise ValueError("Invalid mode")
                scores.append((score, r))
            scores.sort(key=lambda x: x[0], reverse=True)
            top_n = int(len(scores) * 0.25)
            return np.vstack([r for _, r in scores[:top_n]])

        mossy_fibers = select_quartile(data, 'mossy_fibers', mode=selection_mode, percentile=selection_percentile) if use_top_quartile else np.vstack([
            r for i, r in enumerate(data["rate"]) if data["pop"][i] == "mossy_fibers"])
        purk = select_quartile(data, 'purkinje', mode=selection_mode, percentile=selection_percentile) if use_top_quartile else np.vstack([
            r for i, r in enumerate(data["rate"]) if data["pop"][i] == "purkinje"])

        if mossy_fibers.size == 0 or purk.size == 0:
            print("⚠️ Dati insufficienti")
            continue

        n = len(data["time"][0])
        start = int(n * trim_start)
        end = int(n * (1 - trim_end))
        
        # Applica lo stesso trimming anche allo stimolo, se esiste
        if stim_profile:
            stim_profile = stim_profile[start:end]

        if method == "mean":
            m_mean = mossy_fibers[:, start:end].mean(0)
            p_mean = purk[:, start:end].mean(0)
        elif method == "median":
            m_mean = np.median(mossy_fibers[:, start:end], axis=0)
            p_mean = np.median(purk[:, start:end], axis=0)
        elif method.startswith("quantile_"):
            q = float(method.split("_")[1])
            m_mean = np.quantile(mossy_fibers[:, start:end], q, axis=0)
            p_mean = np.quantile(purk[:, start:end], q, axis=0)
        else:
            raise ValueError("Metodo non valido per la stima temporale")

        m_sem, p_sem = mossy_fibers[:, start:end].std(0)/np.sqrt(mossy_fibers.shape[0]), purk[:, start:end].std(0)/np.sqrt(purk.shape[0])

        mossy_fibers_quantile_value = float(np.quantile(m_mean, selection_percentile))
        purkinje_quantile_value = float(np.quantile(p_mean, selection_percentile))
        
        if stim_profile:
            stim_profile = np.array(stim_profile)
            stim = stim_profile[start:end]
        
            # Allinea le lunghezze con m_mean
            if method == "mean":
                m_mean = mossy_fibers[:, start:end].mean(0)
                p_mean = purk[:, start:end].mean(0)
            elif method == "median":
                m_mean = np.median(mossy_fibers[:, start:end], axis=0)
                p_mean = np.median(purk[:, start:end], axis=0)
            elif method.startswith("quantile_"):
                q = float(method.split("_")[1])
                m_mean = np.quantile(mossy_fibers[:, start:end], q, axis=0)
                p_mean = np.quantile(purk[:, start:end], q, axis=0)
            else:
                raise ValueError("Metodo non valido per la stima temporale")
        
            # Allinea tutte le serie alla lunghezza minima
            min_len = min(len(stim), len(m_mean))
            stim = stim[:min_len]
            m_mean = m_mean[:min_len]
            p_mean = p_mean[:min_len]
        
            stim_to_mossy_fibers_cosine = cosine_similarity(m_mean, stim)
            stim_to_purk_cosine = cosine_similarity(p_mean, stim)
        else:
            stim = None
            stim_to_mossy_fibers_cosine = np.nan
            stim_to_purk_cosine = np.nan

        row = {
            "condition": condition_name,
            "k": k,
            'Atrophy Factor (%)': atrophy,
            "atrophy_percent": atrophy,
            "mossy_fibers_mean_rate": float(m_mean.max()), #float(m_mean.mean()),
            "mossy_fibers_quantile_rate": mossy_fibers_quantile_value,
            "mossy_fibers_sem_rate": float(m_sem.mean()),
            "purkinje_mean_rate": float(p_mean.max()), #float(p_mean.mean()),
            "purkinje_quantile_rate": purkinje_quantile_value,
            "purkinje_sem_rate": float(p_sem.mean()),
            "mossy_fibers_timeseries": m_mean.tolist(),
            "purkinje_timeseries": p_mean.tolist(),
            "cosine_similarity_mossy_fibers_purk": float(cosine_similarity(m_mean, p_mean)),
            "cosine_similarity_stim_mossy_fibers": float(stim_to_mossy_fibers_cosine),
            "cosine_similarity_stim_purk": float(stim_to_purk_cosine),
            "cross_corr": compute_cross_corr_and_lag(m_mean, p_mean, dt)[0],
            "optimal_lag_ms": compute_cross_corr_and_lag(m_mean, p_mean, dt)[1],
            "te_mossy_fibers_to_purk": compute_transfer_entropy_binned_old(m_mean, p_mean),
            "te_purk_to_mossy_fibers": compute_transfer_entropy_binned_old(p_mean, m_mean),
            "plv": compute_plv(m_mean, p_mean),
            "intra_corr_mossy_fibers": mean_pairwise_corr(mossy_fibers[:, start:end]),
            "intra_corr_purk": mean_pairwise_corr(purk[:, start:end]),
            "fooof_mossy_fibers": analyze_fooof(m_mean, fs=1000/dt),
            "fooof_purk": analyze_fooof(p_mean, fs=1000/dt),
            "burst_onsets_ms": burst_onsets_ms,
            "burst_durations_ms": burst_durations_ms,
            "burst_rates_hz": burst_rates_hz,
            "stim_profile": stim_profile,
        }

        summary_rows.append(row)
        df_summary = pd.DataFrame(summary_rows)
        df_summary = df_summary[df_summary["atrophy_percent"] < THR].reset_index(drop=True)
        output_dir = Path(base_dir) / "B_collective_results" / "2.functions" / f"sigma_{sigma}"
        output_dir.mkdir(parents=True, exist_ok=True)
        df_summary.to_pickle(output_dir / f"{sigma}_df.pkl")

    
    return df_summary, spike_dfs


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, required=True)
    parser.add_argument("--sigma", type=float, default=100.0, help="Sigma for kernel smoothing in firing rate estimation")
    parser.add_argument("--stimName", type=str, default="basal_activity", help="Stimulation name")
    args = parser.parse_args()
    base_dir = Path(args.path)
    df_summary, spike_dfs = loadData(base_dir, mainEventKey="burst", 
                                     sigma=args.sigma,
                                     stimName=args.stimName,
                                     condition_name="saccades", 
                                     use_top_quartile=True, selection_mode="percentile", method='mean')

if __name__ == "__main__":
    main()
