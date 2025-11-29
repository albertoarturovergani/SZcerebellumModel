#!/usr/bin/env python3

import yaml
import argparse
import numpy as np

def replace_background_noise(k_suffix, spike_file_path):
    yaml_path = f'./configurations/mouse/nest/basal_vitro_{k_suffix}.yaml'

    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f)

    try:
        sim_config = config['simulations']['basal_activity']
        devices = sim_config.get('devices', {})
    except KeyError:
        raise KeyError("❌ Struttura 'simulations > basal_activity > devices' non trovata.")

    # ✅ Carica i tempi dal file (come lista)
    spike_times = np.loadtxt(spike_file_path).tolist()

    # 🔁 Rimuovi background_noise se presente
    if 'background_noise' in devices:
        print("ℹ️  Rimozione 'background_noise' esistente...")
        devices.pop('background_noise')

    # ➕ Inserisci il nuovo dispositivo 'aperiodic'
    devices['aperiodic'] = {
        'device': 'spike_generator',
        'spike_times': spike_times,
        'targetting': {
            'strategy': 'cell_model',
            'cell_models': ['mossy_fibers']
        },
        'weight': 1,
        'delay': 0.1
    }

    # Salva il file aggiornato
    config['simulations']['basal_activity']['devices'] = devices
    with open(yaml_path, 'w') as f:
        yaml.dump(config, f, sort_keys=False)

    print(f"✅ File aggiornato: {yaml_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Sostituisce 'background_noise' con 'aperiodic' nel file YAML.")
    parser.add_argument('k_suffix', help="Suffisso del file YAML (es. 0.400)")
    parser.add_argument('spike_file_path', help="Percorso al file .txt con spike_times")

    args = parser.parse_args()
    replace_background_noise(args.k_suffix, args.spike_file_path)
