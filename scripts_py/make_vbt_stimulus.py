import os
import numpy as np
import yaml

def make_vbt_stimulus(k, seed, stim_name):
    sim_path = "./configurations/mouse/nest/"
    sim_file = os.path.join(sim_path, "basal_vitro.yaml")
    if not os.path.exists(sim_file):
        print(f"❌ ERRORE: File {sim_file} non trovato.")
        return

    with open(sim_file, 'r') as file:
        sim_data = yaml.safe_load(file)

    k_str = f"{k:.3f}"

    # --- Modifica solo simulations ---
    new_sim = {
        stim_name: {
            'simulator': 'nest',
            'resolution': 0.1,
            'duration': 25000,
            'modules': ['cerebmodule'],
            'seed': seed,
            'devices': {}
        }
    }

    # --- Copia i cell_models ---
    if 'simulations' in sim_data and 'basal_activity' in sim_data['simulations']:
        basal = sim_data['simulations']['basal_activity']
        if 'cell_models' in basal:
            new_sim[stim_name]['cell_models'] = basal['cell_models']
        if 'connection_models' in basal:
            new_sim[stim_name]['connection_models'] = basal['connection_models']


    # --- Copia i dispositivi di registrazione ---
    if 'simulations' in sim_data and 'basal_activity' in sim_data['simulations']:
        devices = sim_data['simulations']['basal_activity'].get('devices', {})
        for key, value in devices.items():
            if value.get('device') == 'spike_recorder':
                new_sim[stim_name]['devices'][key] = value


    # --- Aggiunge baseline iniziale ---
    new_sim[stim_name]['devices']['baseline_1'] = {
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

    # --- Aggiunge step 5–13s ---
    rates = np.round(np.random.uniform(5, 100, size=8), 3)
    for i, rate in enumerate(rates):
        start = 5000.0 + i * 1000.0
        stop = start + 1000.0
        stim_key = f"stim_step_{i+1}"
        new_sim[stim_name]['devices'][stim_key] = {
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

    # --- Baseline finale 13–25s ---
    new_sim[stim_name]['devices']['baseline_2'] = {
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

    # --- Stimoli burst 13–21s ---
    n_bursts = 8
    burst_durations = np.round(np.random.uniform(20.0, 120.0, size=n_bursts), 1)
    burst_rates = np.round(6000.0 / burst_durations, 1)
    for i, (dur, rate) in enumerate(zip(burst_durations, burst_rates)):
        burst_start = 13000.0 + i * 1000.0
        label = f"burst_{i+1}"
        new_sim[stim_name]['devices'][label] = {
            'device': 'poisson_generator',
            'rate': float(rate),
            'start': float(burst_start),
            'stop': float(burst_start + dur),
            'targetting': {
                'strategy': 'sphere',
                'radius': 90,
                'origin': [150.0, 65.0, 100.0],
                'cell_models': ['mossy_fibers']
            },
            'weight': 1.0,
            'delay': 0.1
        }

    # --- Assegna la nuova sezione simulations ---
    sim_data['simulations'] = new_sim

    # --- Aggiorna $import e storage ---
    if '$import' in sim_data:
        sim_data['$import']['ref'] = f"../mouse_cerebellar_cortex_{k_str}.yaml#/"
    if 'storage' in sim_data:
        sim_data['storage']['root'] = f"mouse_cerebellar_cortex_{k_str}.hdf5"
        sim_data['storage']['engine'] = "hdf5"

    # --- Salva ---
    stim_outfile = os.path.join(sim_path, f"{stim_name}_{k_str}.yaml")
    with open(stim_outfile, 'w') as file:
        yaml.safe_dump(sim_data, file, default_flow_style=False, allow_unicode=True, sort_keys=False)

    print(f"✅ File stimolo generato: {stim_outfile}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate VBT stimulus YAML.")
    parser.add_argument('--k', type=float, required=True, help='Atrophy factor k (e.g., 0.850)')
    parser.add_argument('--seed', type=int, required=True, help='Random seed')
    parser.add_argument('--stim_name', type=str, required=True, help='Simulation name inside the YAML')

    args = parser.parse_args()
    make_vbt_stimulus(args.k, args.seed, args.stim_name)
