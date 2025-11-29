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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Cartella principale della simulazione")
    args = parser.parse_args()

    # === Parametri epoca ===
    t_pre = 0    # ms prima dell’evento
    t_post = 10  # ms dopo l’evento
    Z_MOhm = 1.0  # Impedenza per stima LFP in MOhm
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
                print(f"    📊 Segnali presenti: {all_signals}")

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
                    lfp = -s * Z_MOhm / 1000.0  # Converti in µV (con segno negativo)

                    for ev in burst_events:
                        if t[0] <= ev <= t[-1]:
                            idx = np.searchsorted(t, ev)
                            if idx - n_pre >= 0 and idx + n_post < len(s):
                                epoch = s[idx - n_pre: idx + n_post]
                                epoch_lfp = lfp[idx - n_pre: idx + n_post]
                                epoched_signals[k].append(epoch)
                                epoched_lfp[k].append(epoch_lfp)

            reader.close()
            del reader, block
            gc.collect()

        if fs is None:
            print("❌ Nessun segnale i_syn trovato in nessun file .nio. Salto questa cartella.")
            continue

        time_epoch = np.linspace(-t_pre, t_post, n_pre + n_post)
        summary_signal = []
        summary_weights = {"I_syn1": 1.0, "I_syn2": 1, "I_syn3": 0, "I_syn4": 0}

        for k in sorted(epoched_signals):
            stack = np.stack(epoched_signals[k])
            if stack.size == 0:
                continue
            mean = np.mean(stack, axis=0)
            sem_ = sem(stack, axis=0)

            plt.figure(figsize=(6, 4))
            plt.plot(time_epoch, mean, label=f"{k} (mean)", lw=2)
            plt.fill_between(time_epoch, mean - sem_, mean + sem_, alpha=0.3)
            plt.title(f"Corrente epocata - {k}")
            plt.xlabel("Tempo (ms)")
            plt.ylabel("Corrente (pA)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig(output_dir / f"{k}_epoched_plot.png", dpi=150)
            plt.close()

            w = summary_weights.get(k, 0)
            summary_signal.append(w * mean)

        if summary_signal:
            total_signal = -np.sum(summary_signal, axis=0)
            plt.figure(figsize=(6, 4))
            plt.plot(time_epoch, total_signal, color='k', lw=2)
            plt.title("Weighted Current")
            plt.xlabel("Tempo (ms)")
            plt.ylabel("Corrente combinata (pA)")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(output_dir / "weighted_sum_I_syn.png", dpi=150)
            plt.close()



 

if __name__ == "__main__":
    main()
