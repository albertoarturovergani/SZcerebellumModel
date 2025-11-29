import numpy as np
import pandas as pd
from tqdm import tqdm
import os
import argparse
from cerebellum.analysis.spiking_results import BasicSimulationReport_v2
from bsb import from_storage
from joblib import Parallel, delayed

def process_simulation(main_dir, file_name, rel, simulation_types):
    """Elabora una simulazione per un valore di rel"""
    rel = f"{rel:.3f}"
    hdf5_file = os.path.join(main_dir, f"results_{rel}", "reconstruction", f"{file_name}_{rel}.hdf5")
    
    if not os.path.exists(hdf5_file):
        print(f"[WARNING] File {hdf5_file} non trovato. Salto iterazione.")
        return []

    print(f"[INFO] Caricamento rete da {hdf5_file}")
    network = from_storage(hdf5_file)
    results = []

    for sim_type in simulation_types:
        sim_dir = os.path.join(main_dir, f"results_{rel}", sim_type)

        if not os.path.exists(sim_dir):
            print(f"[WARNING] Directory {sim_dir} non trovata. Salto simulazione {sim_type}.")
            continue

        print(f"[INFO] Elaborazione simulazione {sim_type} per rel={rel}")
        report = BasicSimulationReport_v2(
            network, 
            simulation_name=sim_type,
            folder_nio=sim_dir
        )
        report.print_report(f"{sim_dir}/{rel}_{sim_type}_dynamics.pdf")

        sim_results_table = report.get_sim_results_table()
        firing_rates = sim_results_table.get_firing_rates()

        for cell_type, values in firing_rates.items():
            mean_fr = np.mean(values) if len(values) > 0 else 0.0
            sem_fr = np.std(values, ddof=1) / np.sqrt(len(values)) if len(values) > 1 else 0.0

            print(f"[INFO] {sim_type} - {cell_type}: Mean FR={mean_fr:.3f}, SEM={sem_fr:.3f}")

            results.append({
                "rel": rel,
                "simulation": sim_type,
                "cell": cell_type,
                "mean_firing_rate": mean_fr,
                "sem_firing_rate": sem_fr
            })
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Processa simulazioni cerebellari.")
    parser.add_argument("main_dir", type=str, help="Directory principale contenente i risultati.")
    parser.add_argument("--file_name", type=str, default="mouse_cereb_dcn_io_nest", help="Nome base del file .hdf5 per ciascun rel.")
    parser.add_argument("--output_name", type=str, default="dynamics_fr_compensation_linear.pkl", help="Nome del file output .pkl")
    args = parser.parse_args()
    
    main_dir = args.main_dir
    simulation_types = ["basal_activity"]
    sweep_range = np.linspace(1, 0.40, 15)
    file_name = args.file_name

    print(f"[INFO] Inizio elaborazione. Directory principale: {main_dir}")
    
    all_results = Parallel(n_jobs=1)(
        delayed(process_simulation)(main_dir, file_name, rel, simulation_types) for rel in tqdm(sweep_range, desc="Processing")
    )
    all_results = [item for sublist in all_results for item in sublist]  # Flatten the list
    
    print("[INFO] Creazione DataFrame e salvataggio dei risultati.")
    df_dynamics = pd.DataFrame(all_results)
    df_dynamics.to_pickle(os.path.join(main_dir, args.output_name))
    
    print("[INFO] Elaborazione completata. File salvato.")

if __name__ == "__main__":
    main()
