#!/usr/bin/env python3

import gc
import neo
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import yaml
from collections import defaultdict
import argparse
from scipy.stats import sem
import re
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.integrate import simpson as simps

lowcut = 0.1  # Hz
highcut = 300  # Hz
order = 4  # ordine del filtro

def bandpass_filter(data, fs, lowcut, highcut, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return filtfilt(b, a, data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Cartella principale della simulazione")
    args = parser.parse_args()

    # === Parametri epoca ===
    t_pre = 0.1    # ms prima dell’evento
    t_post = 15  # ms dopo l’evento
    Z_Ohm = 10.0/2 #*1e1  # Impedenza per stima LFP in MOhm

    # === Cartella principale ===
    root_dir = Path(args.root)

    # === Cerca nelle sottocartelle *_deg_*/3.function ===
    for result_dir in sorted(root_dir.glob('*_deg_*/3.function')):
        print(f"\n📁 Entrata in {result_dir}")
        output_dir = result_dir / "plots_multimeter"
        output_dir.mkdir(exist_ok=True)

        # === Trova l’unico YAML nella cartella 2.stimulus ===
        stim_dir = result_dir.parent / "2.stimulus"
        yaml_files = list(stim_dir.glob("*.yaml"))
        assert len(yaml_files) == 1, f"❌ Attesi 1 solo file YAML, trovati {len(yaml_files)} in {stim_dir}"
        stim_yaml_path = yaml_files[0]

        # === Estrai i tempi di inizio dai burst_* ===
        with open(stim_yaml_path, 'r') as f:
            stim_dict = yaml.safe_load(f)

        devices = stim_dict.get("simulations", {}).get("pawan_stim_vitro", {}).get("devices", {})
        burst_events = []

        for name, dev in devices.items():
            if name.startswith("burst_") and isinstance(dev, dict):
                start = dev.get("start")
                if start is not None:
                    burst_events.append(start)

        burst_events = np.array(sorted(burst_events))
        print(f"✅ Trovati {len(burst_events)} eventi burst da {stim_yaml_path.name}")

        # === Accumulatore epoche ===
        epoched_signals = defaultdict(list)
        epoched_lfp = defaultdict(list)

        fs = None  # sarà calcolata al primo segnale valido

        for file in sorted(result_dir.glob("*.nio")):
            print(f"\n📂 Processing {file.name}")
            reader = neo.io.NixIO(filename=str(file), mode="ro")
            block = reader.read_block()

            for segment in block.segments:
                print(f"  ↪️ Segment: {segment.name}")
                all_signals = [a.name for a in segment.analogsignals]
                print(f"    📊 Segnali presenti: {len(all_signals)}")

                # Seleziona segnali i_syn* in modo flessibile
                analogs = {a.name: a for a in segment.analogsignals if a.name.lower().startswith("i_syn")}

                if not analogs:
                    print(f"⚠️ Nessun segnale i_syn nel segmento {segment.name} di {file.name}")
                    continue

                if fs is None:
                    fs = 1.0 / list(analogs.values())[0].sampling_period.rescale("s").magnitude
                    t = list(analogs.values())[0].times.rescale("ms").magnitude
                    n_pre = int(t_pre / 1000 * fs)
                    n_post = int(t_post / 1000 * fs)

                for k, a in analogs.items():
                    s = a.rescale("pA").magnitude.squeeze()
                    #lfp = -s * Z_MOhm #/ 1000.0  # Converti in µV (con segno negativo)
                    """
                    for ev in burst_events:
                        if t[0] <= ev <= t[-1]:
                            idx = np.searchsorted(t, ev)
                            if idx - n_pre >= 0 and idx + n_post < len(s):
                                epoch = s[idx - n_pre: idx + n_post]
                                #epoch_lfp = lfp[idx - n_pre: idx + n_post]
                                epoched_signals[k].append(epoch)
                                #epoched_lfp[k].append(epoch_lfp)
                    """

                    for ev in burst_events:
                        if t[0] <= ev <= t[-1]:
                            idx = np.searchsorted(t, ev)
                            if idx - n_pre >= 0 and idx + n_post < len(s):
                                epoch = s[idx - n_pre: idx + n_post]
                                baseline = np.mean(epoch[:n_pre])
                                epoch = epoch - baseline  # 🟢 qui correggi il riferimento
                                epoched_signals[k].append(epoch)


            reader.close()
            del reader, block
            gc.collect()

        if fs is None:
            print("❌ Nessun segnale i_syn trovato in nessun file .nio. Salto questa cartella.")
            continue

        #time_epoch = np.linspace(-t_pre, t_post, n_pre + n_post)
        time_epoch = np.linspace(-t_pre, t_post, n_pre + n_post, endpoint=False)

        summary_signal = []
        #summary_weights = {"I_syn1": 1.0, "I_syn2": 1, "I_syn3": 1, "I_syn4": 1}
        summary_weights = {
          "I_syn1": 1.0,   # AMPA (eccitatoria fast)
          "I_syn2": 3.0,   # GABA_A (inibitoria)
          "I_syn3": 0.1,   # NMDA (eccitatoria slow)
          "I_syn4": 0.0    # GABA_A extrasinaptica
        }

        np.save(output_dir / "epoched_signals.npy", epoched_signals)

        for k in sorted(epoched_signals):
            stack = np.stack(epoched_signals[k])
            if stack.size == 0:
                continue
            mean = np.mean(stack, axis=0)
            plt.figure(figsize=(6, 4))
            for trial in stack:
                plt.plot(time_epoch, trial, color='gray', alpha=0.3)
            plt.plot(time_epoch, mean, label=f"{k} (mean)", color='black', lw=2)
            plt.title(f"Epoched Current - {k}")
            plt.xlabel("Time (ms)")
            plt.ylabel("Current (pA)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig(output_dir / f"{k}_epoched_butterfly_plot.png", dpi=150)
            plt.close()
            np.save(output_dir / f"{k}_epoched_butterfly_mean.npy", mean)
            w = summary_weights.get(k, 0)
            summary_signal.append(w * mean)

            # === Plot: tutte le correnti sovrapposte ===
            plt.figure(figsize=(7, 5))
            for k_sig in sorted(epoched_signals):
                stack = np.stack(epoched_signals[k_sig])
                if stack.size == 0:
                    continue
                mean = np.mean(stack, axis=0)
                plt.plot(time_epoch, mean, label=k_sig)
            plt.title("Correnti sinaptiche sovrapposte")
            plt.xlabel("Tempo (ms)")
            plt.ylabel("Corrente (pA)")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(output_dir / "all_currents_overlaid.png", dpi=150)
            plt.close()
    
            # === Plot: trasformazioni in Volt delle correnti sovrapposte ===
            plt.figure(figsize=(7, 5))
            for k_sig in sorted(epoched_signals):
                stack = np.stack(epoched_signals[k_sig])
                if stack.size == 0:
                    continue
                mean = np.mean(stack, axis=0)
                weight = summary_weights.get(k_sig, 0.0)
                signal_volt = -mean * weight / (4 * np.pi * Z_Ohm) / 1e6  # in Volt
                plt.plot(time_epoch, signal_volt, label=f"{k_sig} (in V)")
            plt.title("Correnti trasformate in LFP (V)")
            plt.xlabel("Tempo (ms)")
            plt.ylabel("Tensione (V)")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(output_dir / "all_currents_transformed_in_V.png", dpi=150)
            plt.close()


        """
        for k in sorted(epoched_signals):
            stack = np.stack(epoched_signals[k])
            if stack.size == 0:
                continue
            mean = np.mean(stack, axis=0)
            sem_ = sem(stack, axis=0)
            plt.figure(figsize=(6, 4))
            plt.plot(time_epoch, mean, label=f"{k} (mean)", lw=2)
            plt.fill_between(time_epoch, mean - sem_, mean + sem_, alpha=0.3)
            plt.title(f"Epoched Current - {k}")
            plt.xlabel("Time (ms)")
            plt.ylabel("Current (pA)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig(output_dir / f"{k}_epoched_plot.png", dpi=150)
            plt.close()
            np.save(output_dir / f"{k}_epoched_mean.npy", mean)
            w = summary_weights.get(k, 0)
            summary_signal.append(w * mean)
        """
        if summary_signal:
            total_signal = -np.sum(summary_signal, axis=0)/(4*np.pi*Z_Ohm)
            baseline = np.mean(total_signal[:n_pre])
            total_signal = total_signal - baseline
            plt.figure(figsize=(6, 4))
            plt.plot(time_epoch, total_signal, color='k', lw=2)
            plt.title(f"Epoched LFP")
            plt.xlabel("Time (ms)")
            plt.ylabel("LFP (uV)")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(output_dir / "weighted_sum_I_syn.png", dpi=150)
            plt.close()
            np.save(output_dir / "weighted_sum_I_syn.npy", total_signal)


        # === Calcolo N2a, N2b e integrale ===
        from scipy.signal import find_peaks
        
        lfp_path = output_dir / "weighted_sum_I_syn.npy"
        if lfp_path.exists():
            lfp = np.load(lfp_path)
            
            # Inverti per cercare minimi
            inv_lfp = -lfp
            peaks, _ = find_peaks(inv_lfp, distance=5)
            
            if len(peaks) >= 2:
                n2a_idx = peaks[0]
                n2b_idx = peaks[1]
                n2a_amp = lfp[n2a_idx]
                n2b_amp = lfp[n2b_idx]
                n2a_lat = float(time_epoch[n2a_idx])
                n2b_lat = float(time_epoch[n2b_idx])
            else:
                n2a_amp = n2b_amp = n2a_lat = n2b_lat = None
        
            integral = float(np.trapz(lfp, time_epoch))
        
            # Estrai valore di k dal path
            k_match = re.search(r"deg_(\d+\.\d+)", str(result_dir))
            k_val = float(k_match.group(1)) if k_match else None
        
            # Costruisci dizionario dei risultati
            out_dict = {
                "k": k_val,
                "N2a_amp": float(n2a_amp) if n2a_amp is not None else None,
                "N2a_lat": n2a_lat,
                "N2b_amp": float(n2b_amp) if n2b_amp is not None else None,
                "N2b_lat": n2b_lat,
                "LFP_integral": integral
            }
        
            # Salva YAML locale
            with open(output_dir / "lfp_metrics.yaml", "w") as f:
                yaml.dump(out_dict, f, sort_keys=False)


    collective_dir = root_dir / "B_collective_results/2.functions"
    collective_dir.mkdir(parents=True, exist_ok=True)
    deg_dirs = sorted(root_dir.glob("*_deg_*/3.function"))

    all_k = []
    all_time = time_epoch
    all_currents = {}
    all_sums = {}

   
    for d in deg_dirs:
        k_match = re.search(r"deg_(\d+\.\d+)", str(d))
        if not k_match:
            continue
        k_val = float(k_match.group(1))
        out_dir = d / "plots_multimeter"
        if not out_dir.exists():
            continue

        for f in sorted(out_dir.glob("I_syn*.npy")):
            name = f.stem.replace("_epoched_mean", "")
            data = np.load(f)
            if all_time is None:
                all_time = np.linspace(0, len(data) - 1, len(data))
            if name not in all_currents:
                all_currents[name] = {}
            all_currents[name][k_val] = data

        weighted_path = out_dir / "weighted_sum_I_syn.npy"
        if weighted_path.exists():
            all_sums[k_val] = np.load(weighted_path)


    # Colormap e normalizzazione
    cmap = cm.get_cmap('viridis_r')
    norm = mcolors.Normalize(vmin=min(k for kd in all_currents.values() for k in kd), 
                             vmax=max(k for kd in all_currents.values() for k in kd))
    
    for name, kdict in all_currents.items():
        plt.figure(figsize=(6, 4))
        for k, y in sorted(kdict.items()):
            color = cmap(norm(k))
            plt.plot(all_time, y, label=f"k={k:.3f}", color=color)
        plt.title(f"{name} vs k")
        plt.xlabel("Tempo (ms)")
        plt.ylabel("Corrente (pA)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(collective_dir / f"{name}_vs_k.png", dpi=150)
        plt.close()
    
    # Colormap anche per la somma pesata
    plt.figure(figsize=(6, 4))
    for k, y in sorted(all_sums.items()):
        color = cmap(norm(k))
        plt.plot(all_time, y, label=f"k={k:.3f}", color=color)
    plt.ylabel("LFP (uV)")
    plt.xlabel("Time (ms)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(collective_dir / "weighted_sum_vs_k.png", dpi=150)
    plt.close()

    import pandas as pd
    
    # === Raccogli tutti i file lfp_metrics.yaml ===
    yaml_paths = list(root_dir.glob("*_deg_*/3.function/plots_multimeter/lfp_metrics.yaml"))
    metrics = []
    
    for ypath in yaml_paths:
        with open(ypath, "r") as f:
            d = yaml.safe_load(f)
            metrics.append(d)
    
    if metrics:
        df = pd.DataFrame(metrics).sort_values("k")
        df_path = collective_dir / "lfp_metrics_summary.csv"
        df.to_csv(df_path, index=False)
    
        # === Plotta ogni feature vs k ===
        features = ["N2a_amp", "N2b_amp", "N2a_lat", "N2b_lat", "LFP_integral"]
        units = {
            "N2a_amp": "uV", "N2b_amp": "uV",
            "N2a_lat": "ms", "N2b_lat": "ms",
            "LFP_integral": "uV·ms"
        }
    
        for feat in features:
            plt.figure(figsize=(6, 4))
            plt.plot(df["k"], df[feat], 'o-', lw=2)
            plt.xlabel("k (atrofia)")
            plt.ylabel(f"{feat} ({units[feat]})")
            plt.title(f"{feat} vs k")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(collective_dir / f"{feat}_vs_k.png", dpi=150)
            plt.close()
    
      

 

if __name__ == "__main__":
    main()
