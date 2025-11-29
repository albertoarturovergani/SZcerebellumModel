import os
import numpy as np
import yaml

def make_pawan_stimulus(k, seed, stim_name):
    
    sim_path = "./configurations/mouse/nest/"
    sim_file = os.path.join(sim_path, "basal_vitro.yaml")
    if not os.path.exists(sim_file):
        print(f"❌ ERRORE: File {sim_file} non trovato.")
        return

    with open(sim_file, 'r') as file:
        sim_data = yaml.safe_load(file)

    k_str = f"{k:.3f}"

    resolution = 0.050 #5 instead of 0.055
    T = 1000 # 1000 ms
    simDuration = round(T/resolution)*resolution
    
    # --- Modifica solo simulations ---
    new_sim = {
        stim_name: {
            'simulator': 'nest',
            'resolution': resolution, # -> 18KHz
            'duration': simDuration,
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
    """
        new_sim[stim_name]['devices']['I_syn3'] = {
            'device': 'multimeter',
            'properties': ['I_syn3'],
            'units': ['pA'],
            'delay': 0.1,
            'targetting': {
                'strategy': 'cell_model',
                'cell_models': ['golgi_cell'], #['golgi_cell'] # add golgi_cells
            }
        }
    """
    # Multimetri per I_syn1–I_syn4 su cellule granulari e Golgi, selezionate in una sfera centrata in [150, 65, 100] con raggio 90
    for i in range(1, 5):
        new_sim[stim_name]['devices'][f'I_syn{i}'] = {
            'device': 'multimeter',
            'properties': [f'I_syn{i}'],
            'units': ['pA'],
            'delay': 0.1,
            'targetting': {
                'strategy': 'sphere',
                'radius': 90,  # origin e radius scelti per stimolare mossy_fibers in entrambi i setup (con/senza DCN).
                'origin': [150.0, 65.0, 100.0],  # radius 90 copre z ≈ 95–335 → evita warning 'has no targets'
                'cell_models': ['golgi_cell', 'granule_cell']
            },
        }

   
    # Aggiunta dei n stimoli "burst"
    # Riproduce la relazione durata saccade ↔ durata burst ↔ firing rate (Kase 1980)
    n_bursts = 30  
    burst_durations = [resolution*20] * n_bursts # instead of 0.200 ms   np.round(np.random.uniform(20.0, 120.0, size=n_bursts), 1)
    burst_rates = [10000.0] * n_bursts

    """
    for i, (dur, rate) in enumerate(zip(burst_durations, burst_rates)):
        burst_start = 200.0 + i * 20
        label = f"burst_{i+1}"
        new_sim[stim_name]['devices'][label] = {
            'device': 'poisson_generator',
            'rate': float(rate),
            'start': float(burst_start),
            'stop': float(round((burst_start + dur) / resolution) * resolution), # per avere multipli della resolution
            'targetting': {
                'strategy': 'sphere',
                'radius': 90, #120, #120, # qui cambiato per stimolare tutte le fibre anche con dcn
                'origin': [150.0, 65.0, 100.0], #[150.0, 65.0, 215.0], # qui cambiato per stimolare tutte le fibre anche con dcn
                # origin e radius scelti per stimolare mossy_fibers in entrambi i setup (con/senza DCN).
                # radius 120 copre z ≈ 95–335 → evita warning 'has no targets'.
                'cell_models': ['mossy_fibers']
            },
            'weight': 1.0,
            'delay': 0.1
        }
    """
    I_amp = 50*1e-6*1e3 # from microA to nanoA
    for i, (dur, rate) in enumerate(zip(burst_durations, burst_rates)):
        burst_start = 200.0 + i * 20
        label = f"burst_{i+1}"
        new_sim[stim_name]['devices'][label] = {
            'device': 'poisson_generator', #'dc_generator',
            #'amplitude': I_amp,
            'rate': float(rate),
            'start': float(burst_start),
            'stop': float(round((burst_start + dur) / resolution) * resolution), # per avere multipli della resolution
            'targetting': {
                'strategy': 'sphere',
                'radius': 90, #120, #120, # qui cambiato per stimolare tutte le fibre anche con dcn
                'origin': [150.0, 65.0, 100.0], #[150.0, 65.0, 215.0], # qui cambiato per stimolare tutte le fibre anche con dcn
                # origin e radius scelti per stimolare mossy_fibers in entrambi i setup (con/senza DCN).
                # radius 120 copre z ≈ 95–335 → evita warning 'has no targets'.
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

    parser = argparse.ArgumentParser(description="Generate MEA stimulus YAML.")
    parser.add_argument('--k', type=float, required=True, help='Atrophy factor k (e.g., 0.850)')
    parser.add_argument('--seed', type=int, required=True, help='Random seed')
    parser.add_argument('--stim_name', type=str, required=True, help='Simulation name inside the YAML')

    args = parser.parse_args()
    make_pawan_stimulus(args.k, args.seed, args.stim_name)
