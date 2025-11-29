import os
import yaml
import numpy as np
import argparse

def vbt_stim_protocol(basal_path, seed, output_path, k_suffix):
    np.random.seed(seed)

    # Carica il file di stimolo base
    with open(basal_path, 'r') as file:
        stim_data = yaml.safe_load(file)

    # ✅ Rinomina 'basal_activity' → 'vbt_protocol_{k_suffix}'
    if 'basal_activity' in stim_data['simulations']:
        stim_data['simulations'][f'vbt_protocol_{k_suffix}'] = stim_data['simulations'].pop('basal_activity')
    else:
        raise KeyError("❌ 'basal_activity' not found in the YAML file!")

    # Lavora su vbt_protocol_{k_suffix}
    stim_block = stim_data['simulations'][f'vbt_protocol_{k_suffix}']
    stim_block['duration'] = 25000
    stim_block.setdefault('devices', {})

    # Rimuove background_noise se presente
    stim_block['devices'].pop('background_noise', None)

    devices = {}

    # baseline_1: 0–5000 ms a 4 Hz
    devices['baseline_1'] = {
        'device': 'poisson_generator',
        'rate': 4,
        'start': 0.0,
        'stop': 5000.0,
        'targetting': {
            'strategy': 'cell_model',
            'cell_models': ['mossy_fibers']
        },
        'weight': 1,
        'delay': 0.1
    }

    # stim_step: 8 stimoli da 5000 a 13000 ms
    rates = np.round(np.random.uniform(5, 100, size=8), 3)
    for i, rate in enumerate(rates):
        start = 5000.0 + i * 1000.0
        stop = start + 1000.0
        devices[f'stim_step_{i+1}'] = {
            'device': 'poisson_generator',
            'rate': float(rate),
            'start': float(start),
            'stop': float(stop),
            'targetting': {
                'strategy': 'cell_model',
                'cell_models': ['mossy_fibers']
            },
            'weight': 1,
            'delay': 0.1
        }

    # baseline_2: 13000–25000 ms a 4 Hz
    devices['baseline_2'] = {
        'device': 'poisson_generator',
        'rate': 4,
        'start': 13000.0,
        'stop': 25000.0,
        'targetting': {
            'strategy': 'cell_model',
            'cell_models': ['mossy_fibers']
        },
        'weight': 1,
        'delay': 0.1
    }

    # burst: 8 bursts
    n_bursts = 8
    burst_durations = np.round(np.random.uniform(20.0, 120.0, size=n_bursts), 1)
    burst_rates = np.round(6000.0 / burst_durations, 1)

    for i in range(n_bursts):
        burst_start = 13000.0 + i * 1000.0
        devices[f'burst_{i+1}'] = {
            'device': 'poisson_generator',
            'rate': float(burst_rates[i]),
            'start': float(burst_start),
            'stop': float(burst_start + burst_durations[i]),
            'targetting': {
                'strategy': 'sphere',
                'radius': 90,
                'origin': [150.0, 65.0, 100.0],
                'cell_models': ['mossy_fibers']
            },
            'weight': 1.0,
            'delay': 0.1
        }

    # Assegna i nuovi devices
    stim_block['devices'] = devices

    # Crea output directory
    os.makedirs(output_path, exist_ok=True)

    # ✅ Salva il nuovo stimolo
    stim_file = os.path.join(output_path, f'vbt_stim_protocol_{k_suffix}.yaml')
    with open(stim_file, 'w') as file:
        yaml.safe_dump(stim_data, file, default_flow_style=False, allow_unicode=True)
    print(f"✅ Stimulus YAML saved at {stim_file}")

    # ✅ Crea anche il file di configurazione per il reconfigure
    config_vbt = {
        "name": f"DCN-IO in vitro with VBT stimulus - rel={k_suffix}",
        "$import": {
            "ref": "./dcn_io.yaml#/",
            "values": [
                "packages", "storage", "network", "regions", "partitions",
                "morphologies", "cell_types", "placement", "connectivity",
                "after_connectivity", "components"
            ]
        },
        "simulations": {
            f"vbt_protocol_{k_suffix}": {
                "$import": {
                    "ref": f"./2.stimulus/vbt_stim_protocol_{k_suffix}.yaml#/simulations/vbt_protocol_{k_suffix}"
                }
            }
        }
    }

    config_file = os.path.join(output_path, f'dcn_io_vbt_{k_suffix}.yaml')
    with open(config_file, 'w') as file:
        yaml.safe_dump(config_vbt, file, default_flow_style=False, allow_unicode=True)
    print(f"✅ Configuration YAML saved at {config_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--basal', type=str, required=True, help="Path to basal_vitro.yaml")
    parser.add_argument('--seed', type=int, required=True, help="Seed for random generation")
    parser.add_argument('--output', type=str, required=True, help="Output directory for adapted YAMLs")
    parser.add_argument('--k', type=str, required=True, help="Normalized k value (e.g., 0.850)")
    args = parser.parse_args()

    vbt_stim_protocol(args.basal, args.seed, args.output, args.k)
