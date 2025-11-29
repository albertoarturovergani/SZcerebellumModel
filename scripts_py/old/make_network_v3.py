import os
import sys
import yaml
import argparse
import warnings
from datetime import datetime
import numpy as np

warnings.filterwarnings("ignore", category=DeprecationWarning)

current_time = datetime.now()
formatted_time = current_time.strftime("%Y%m%d%H%M%S")
print(f"Script started at: {formatted_time}")

def convert_np_types(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.float32, np.float64)):
        return float(obj)
    elif isinstance(obj, (np.int32, np.int64)):
        return int(obj)
    elif isinstance(obj, dict):
        return {k: convert_np_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_np_types(v) for v in obj]
    else:
        return obj

def safe_load_yaml(filepath):
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)
    if data is None:
        raise ValueError(f"❌ Error: file {filepath} is empty or invalid.")
    return data

def k_reserve_simple(k, k0=0.75, sharpness=1e-3):
    S = 1 / (1 + np.exp(-sharpness * (k - k0)))
    S_min = 1 / (1 + np.exp(-sharpness * (0.5 - k0)))
    S_max = 1 / (1 + np.exp(-sharpness * (1.0 - k0)))
    return 0.5 + 0.5 * (S - S_min) / (S_max - S_min)

def compute_compensatory_factor(k, R, strategy='linear'):
    if strategy == "linear":
        return R * (1-k)


def writeYAML_cerebellum_cx(k, do_basal=False, reserve=[0.75, 1e-3], compensation=1):
    base_path = "/home/alberto/cerebellum/cerebellum/"
    morpho_path = os.path.join(base_path, "morphologies/")
    config_path = os.path.join(base_path, "configurations/mouse/")
    sim_path = os.path.join(base_path, "configurations/mouse/nest/")
    CELLS = ['GranuleCell', 'GolgiCell', 'PurkinjeCell', 'StellateCell', 'BasketCell']
    k_str = f"{k:.3f}"

    # --- 1. Update morphologies.yaml ---
    morpho_file = os.path.join(morpho_path, "morphologies.yaml")
    morpho_data = safe_load_yaml(morpho_file)

    for i, cell in enumerate(CELLS):
        if i < len(morpho_data['morphologies']):
            morpho_data['morphologies'][i]['file'] = f'/home/alberto/cerebellum/cerebellum/morphologies/{cell}/{cell}_{k_str}_atr.swc'

    morpho_outfile = os.path.join(morpho_path, f"morphologies_{k_str}.yaml")
    with open(morpho_outfile, 'w') as f:
        yaml.safe_dump(convert_np_types(morpho_data), f, default_flow_style=False, allow_unicode=True)

    print(f"✅ Created {morpho_outfile}")

    # --- 2. Update mouse_cerebellar_cortex.yaml ---
    config_file = os.path.join(config_path, "mouse_cerebellar_cortex.yaml")
    config_data = safe_load_yaml(config_file)

    if '$import' in config_data:
        config_data['$import']['ref'] = f"/home/alberto/cerebellum/cerebellum/morphologies/morphologies_{k_str}.yaml#/"
    config_data['storage']['root'] = f"mouse_cerebellar_cortex_{k_str}.hdf5"

    if 'after_connectivity' in config_data:
        config_data['after_connectivity']['print_structure_report']['output_filename'] = f'bsb_report_structure_{k_str}.pdf'

    for i, cell in enumerate(['granule_cell', 'golgi_cell', 'purkinje_cell', 'stellate_cell', 'basket_cell']):
        if 'cell_types' in config_data and cell in config_data['cell_types']:
            config_data['cell_types'][cell]['spatial']['morphologies'] = [f'{CELLS[i]}_{k_str}_atr']
            
    ################
    # NETWORK
    # reserve
    k_modulated = k_reserve_simple(k, k0=reserve[0], sharpness=reserve[1])
    if 'cell_types' in config_data:
        for cell, props in config_data['cell_types'].items():
            spatial = props.get('spatial', {})
            if 'density' in spatial:
                spatial['density'] = np.round(float(spatial['density'] * k_modulated), 6)
            if 'planar_density' in spatial:
                spatial['planar_density'] = np.round(float(spatial['planar_density'] * k_modulated), 6)

    if 'partitions' in config_data:
        for layer in config_data['partitions']:
            for key, value in config_data['partitions'][layer].items():
                if isinstance(value, (int, float)):
                    config_data['partitions'][layer][key] = np.round(float(value * k_modulated), 6)

    # compensation
    if 'connectivity' in config_data:
        for connection in config_data['connectivity']:
            if 'contacts' in config_data['connectivity'][connection]:
                if 'loc' in config_data['connectivity'][connection]['contacts']:
                    original_loc = config_data['connectivity'][connection]['contacts']['loc']
                    factor = 1 + compute_compensatory_factor(k, compensation, strategy='linear')
                    new_loc = original_loc * factor
                    config_data['connectivity'][connection]['contacts']['loc'] = np.round(new_loc, 6)
                if 'scale' in config_data['connectivity'][connection]['contacts']:
                    original_loc = config_data['connectivity'][connection]['contacts']['scale']
                    factor = 1 + compute_compensatory_factor(k, compensation, strategy='linear')
                    new_loc = original_loc * factor
                    config_data['connectivity'][connection]['contacts']['scale'] = np.round(new_loc, 6)

    #########################
    config_outfile = os.path.join(config_path, f"mouse_cerebellar_cortex_{k_str}.yaml")
    with open(config_outfile, 'w') as f:
        yaml.safe_dump(convert_np_types(config_data), f, default_flow_style=False, allow_unicode=True)

    print(f"✅ Created {config_outfile}")

    # --- 3. (Optional) Duplicate basal_vitro.yaml stimulus ---
    if do_basal:
        sim_file = os.path.join(sim_path, "basal_vitro.yaml")
        sim_data = safe_load_yaml(sim_file)
        if '$import' in sim_data:
            sim_data['$import']['ref'] = f"../mouse_cerebellar_cortex_{k_str}.yaml#/"
        sim_data['storage']['root'] = f"mouse_cerebellar_cortex_{k_str}.hdf5"

        sim_outfile = os.path.join(sim_path, f"basal_vitro_{k_str}.yaml")
        with open(sim_outfile, 'w') as f:
            yaml.safe_dump(convert_np_types(sim_data), f, default_flow_style=False, allow_unicode=True)

        print(f"✅ Created cortex stimulus {sim_outfile}")


def main():
    parser = argparse.ArgumentParser(description="Generate network YAMLs with a specific atrophy factor k.")
    parser.add_argument('--k', type=float, required=True, help="atrophy factor k (e.g., 0.85)")
    parser.add_argument('--do_basal', action='store_true', help="also generate basal_vitro yaml")
    parser.add_argument('--reserve', nargs=2, type=float, default=[0.75, 1e-3],
                        metavar=('k0', 'sharpness'), help="reserve function parameters (default: 0.75 1e-3)")
    parser.add_argument('--compensation', type=float, default=1.0,
                        help="compensatory factor R for connectivity scaling (default: 1.0)")
    args = parser.parse_args()

    k = args.k
    do_basal = args.do_basal
    reserve = args.reserve
    compensation = args.compensation

    if not (0 < k <= 1):
        print("❌ ERROR: k must be in (0, 1].")
        sys.exit(1)

    writeYAML_cerebellum_cx(k, do_basal=do_basal, reserve=reserve, compensation=compensation)

if __name__ == "__main__":
    main()
