import os
import sys
import yaml
import argparse
import warnings
from datetime import datetime
import numpy as np
import ast

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
    epsilon = 1e-9
    S = 1 / (1 + np.exp(-sharpness * (k - k0)))
    S_min = 1 / (1 + np.exp(-sharpness * (0.5 - k0)))
    S_max = 1 / (1 + np.exp(-sharpness * (1.0 - k0)))
    return 0.5 + 0.5 * (S - S_min) / (S_max - S_min + epsilon)

def compute_compensatory_factor(k, R, strategy='full'):
    if strategy == "conservative":
        return 1 + (R * (1-k))
    if strategy == "semiexp":
        return R * 1/k
    if strategy == "full":
        return ((1/k)**2)**R

def scale_connectivity_params(connectivity_block, k, R):
    factor = compute_compensatory_factor(k, R)
    if 'contacts' in connectivity_block:
        if 'loc' in connectivity_block['contacts']:
            connectivity_block['contacts']['loc'] *= factor
        if 'scale' in connectivity_block['contacts']:
            connectivity_block['contacts']['scale'] *= factor
    if 'affinity' in connectivity_block:
        connectivity_block['affinity'] = min(connectivity_block['affinity'] * factor, 1.0)

def update_erev_gaba(constants: dict, k: float, key: str = "E_rev2"):
    if key in constants:
        original = constants[key]
        if isinstance(original, (int, float)) and original < 0:
            updated = original + (min(-65 - original, 0) * (1 - k))
            constants[key] = float(np.clip(updated, original, -65))

def scale_weights(connection_models: dict, k: float, listConns: list, effect='upregulation'):
    print("\n🧮 [SCALING PESI SINAPTICI INIBITORI]")
    for conn_name, conn_data in connection_models.items():
        if conn_name in listConns:
            weight = conn_data['synapse'].get('weight', 1.0)
            if effect == 'upregulation':
                new_weight = weight * (1 + k)
            else:
                new_weight = weight * k

            conn_data['synapse']['weight'] = float(np.round(new_weight, 6))
            change = new_weight - weight
            perc_change = (change / weight) * 100
            inverse_factor = weight / new_weight if new_weight != 0 else np.inf

            print(f"🔸 {conn_name}")
            print(f"    Weight:    {weight:.4f} → {new_weight:.4f} ({perc_change:+.2f}%)")
            print(f"    ➕ Compensatory factor needed on excitation: ×{inverse_factor:.3f}\n")

def structural_compensation(connectivity_block, k, threshold=0.999, compensation=0):
    if k >= 1.0 or compensation == 0.0:
        return
    factor = (1 / k) ** compensation
    print(f"\n🛠️  Applying structural compensation for k = {k:.3f} → factor ×{factor:.2f}")
    for key in ['divergence', 'convergence', 'radius', 'voxels_pre', 'voxels_post']:
        if key in connectivity_block:
            old_val = connectivity_block[key]
            new_val = round(old_val * factor, 6)
            connectivity_block[key] = new_val
            print(f"   🔧 {key}: {old_val:.4f} → {new_val:.4f}")
    if 'contacts' in connectivity_block:
        for subkey in ['loc', 'scale']:
            if subkey in connectivity_block['contacts']:
                old_val = connectivity_block['contacts'][subkey]
                new_val = round(old_val * factor, 6)
                connectivity_block['contacts'][subkey] = new_val
                print(f"   🔧 contacts.{subkey}: {old_val:.4f} → {new_val:.4f}")

def filter_connections(hypothesis):
    conn_properties = {
        'golgi_to_glomerulus': ('gaba', 'golgi'),
        'golgi_to_golgi': ('gaba', 'golgi'),
        'stellate_to_stellate': ('gaba', 'stellate'),
        'stellate_to_purkinje': ('gaba', 'stellate'),
        'basket_to_basket': ('gaba', 'basket'),
        'basket_to_purkinje': ('gaba', 'basket'),
        'parallel_fiber_to_purkinje': ('glutamate', 'granule'),
        'ascending_axon_to_purkinje': ('glutamate', 'granule'),
        'parallel_fiber_to_golgi': ('glutamate', 'granule'),
        'parallel_fiber_to_basket': ('glutamate', 'granule'),
        'parallel_fiber_to_stellate': ('glutamate', 'granule'),
        'glomerulus_to_granule': ('glutamate', 'glomerulus'),
        'glomerulus_to_golgi': ('glutamate', 'glomerulus'),
        'mossy_fibers_to_glomerulus': ('glutamate', 'mf'),
        'ascending_axon_to_golgi': ('glutamate', 'granule'),
        'gap_goc': ('electrical', 'golgi'),
    }
    transmitter = hypothesis.get('trasmitter')
    neuron_type = hypothesis.get('neurons')
    filtered = []
    for conn, (tx, src) in conn_properties.items():
        if transmitter and tx != transmitter:
            continue
        if neuron_type != 'all' and neuron_type != src:
            continue
        filtered.append(conn)
    return filtered, conn_properties
    
def writeYAML_cerebellum_cx(k, 
                            do_basal=True, 
                            reserve=[0.75, 1e-3], 
                            compensation=1.0, 
                            sz_hyphothesis={'trasmitter': 'glutamate', 'neurons': 'all'}):

    base_path = "/home/alberto/cerebellum/cerebellum/"
    morpho_path = os.path.join(base_path, "morphologies/")
    config_path = os.path.join(base_path, "configurations/mouse/")
    sim_path = os.path.join(base_path, "configurations/mouse/nest/")
    CELLS = ['GranuleCell', 'GolgiCell', 'PurkinjeCell', 'StellateCell', 'BasketCell']
    k_str = f"{k:.3f}"

    k_modulated = k_reserve_simple(k, k0=reserve[0], sharpness=reserve[1])

    # 1. MORPHOLOGY
    morpho_file = os.path.join(morpho_path, "morphologies.yaml")
    morpho_data = safe_load_yaml(morpho_file)
    for i, cell in enumerate(CELLS):
        if i < len(morpho_data['morphologies']):
            morpho_data['morphologies'][i]['file'] = f'{morpho_path}{cell}/{cell}_{k_str}_atr.swc'
    morpho_outfile = os.path.join(morpho_path, f"morphologies_{k_str}.yaml")
    with open(morpho_outfile, 'w') as f:
        yaml.safe_dump(convert_np_types(morpho_data), f, default_flow_style=False, allow_unicode=True)
    print(f"✅ Created {morpho_outfile}")

    # 2. CONFIGURAZIONE RETE
    config_file = os.path.join(config_path, "mouse_cerebellar_cortex.yaml")
    config_data = safe_load_yaml(config_file)

    # 3. CARICAMENTO PESI e loro SCALING
    weight_file_in = os.path.join(sim_path, "basal_vitro.yaml")
    weight_config_data = safe_load_yaml(weight_file_in)
    A = weight_config_data['simulations']['basal_activity']['connection_models']

    # 4. Modifica riferimenti nella rete
    if '$import' in config_data:
        config_data['$import']['ref'] = f"{morpho_path}morphologies_{k_str}.yaml#/"
    config_data['storage']['root'] = f"mouse_cerebellar_cortex_{k_str}.hdf5"
    if 'after_connectivity' in config_data:
        config_data['after_connectivity']['print_structure_report']['output_filename'] = f'bsb_report_structure_{k_str}.pdf'

    for i, cell in enumerate(['granule_cell', 'golgi_cell', 'purkinje_cell', 'stellate_cell', 'basket_cell']):
        if 'cell_types' in config_data and cell in config_data['cell_types']:
            config_data['cell_types'][cell]['spatial']['morphologies'] = [f'{CELLS[i]}_{k_str}_atr']

    # 5. RIDUZIONE NEURONI (RESERVE FUNCTION)
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

    # 6. FILTRAGGIO CONNESSIONI E SCALING
    listConns, conn_properties = filter_connections(sz_hyphothesis)
    print(f"\n📌 [SELEZIONE CONNESSIONI da SZ_hyphothesis]:\n   → {sz_hyphothesis}")
    print(f"   → Connessioni selezionate: {listConns}\n")
    scale_weights(A, k_modulated, listConns, effect='downregulation')

    for conn_name in config_data['connectivity']:
        if conn_name in listConns:
            conn_block = config_data['connectivity'][conn_name]
            structural_compensation(conn_block, k_modulated, threshold=0.999, compensation=compensation)


    # 7. SALVATAGGIO CONFIG RETE
    config_outfile = os.path.join(config_path, f"mouse_cerebellar_cortex_{k_str}.yaml")
    with open(config_outfile, 'w') as f:
        yaml.safe_dump(convert_np_types(config_data), f, default_flow_style=False, allow_unicode=True)
    print(f"✅ Created {config_outfile}")

    # 8. SALVATAGGIO NUOVO FILE PESI
    weights_outfile = os.path.join(sim_path, f"basal_vitro_{k_str}.yaml")
    with open(weights_outfile, 'w') as f:
        yaml.safe_dump(convert_np_types(weight_config_data), f, default_flow_style=False, allow_unicode=True)
    print(f"✅ Saved scaled weights to {weights_outfile}")

    # 9. SALVATAGGIO STIMOLAZIONE (SE RICHIESTA)
    if do_basal:
        sim_data = safe_load_yaml(weights_outfile)  # usa i pesi scalati salvati sopra
        if '$import' in sim_data:
            sim_data['$import']['ref'] = f"../mouse_cerebellar_cortex_{k_str}.yaml#/"
        sim_data['storage']['root'] = f"mouse_cerebellar_cortex_{k_str}.hdf5"
        sim_outfile = os.path.join(sim_path, f"basal_vitro_{k_str}.yaml")
        with open(sim_outfile, 'w') as f:
            yaml.safe_dump(convert_np_types(sim_data), f, default_flow_style=False, allow_unicode=True)
        print(f"✅ Created cortex stimulus {sim_outfile}")

    return config_data, weight_config_data


def main():
    parser = argparse.ArgumentParser(description="Generate network YAMLs with a specific atrophy factor k.")
    parser.add_argument('--k', type=float, required=True, help="atrophy factor k (e.g., 0.85)")
    parser.add_argument('--do_basal', action='store_true', help="also generate basal_vitro yaml")
    parser.add_argument('--reserve', nargs=2, type=float, default=[0.75, 1e-3], metavar=('k0', 'sharpness'))
    parser.add_argument('--compensation', type=float, default=1.0, help="compensatory factor R for connectivity scaling")
    parser.add_argument('--sz_hyphothesis', type=str,
                        default="{'trasmitter': 'gaba', 'neurons': 'all'}",
                        help="dictionary as string, e.g. \"{'trasmitter': 'gaba', 'neurons': 'basket'}\"")

    args = parser.parse_args()

    if not (0 < args.k <= 1):
        print("❌ ERROR: k must be in (0, 1].")
        sys.exit(1)

    try:
        sz_dict = ast.literal_eval(args.sz_hyphothesis)
        assert isinstance(sz_dict, dict)
    except Exception as e:
        print(f"❌ Invalid --sz_hyphothesis format: {e}")
        sys.exit(1)

    writeYAML_cerebellum_cx(
        k=args.k,
        do_basal=args.do_basal,
        reserve=args.reserve,
        compensation=args.compensation,
        sz_hyphothesis=sz_dict
    )


if __name__ == "__main__":
    main()


