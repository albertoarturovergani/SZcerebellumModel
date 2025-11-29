import argparse
import pandas as pd
import numpy as np
from bsb.core import from_storage
import yaml
import os

def extract_structure_summary(hdf5_path):
    scaffold = from_storage(hdf5_path)

    output_dir = os.path.dirname(hdf5_path)
    basename = os.path.splitext(os.path.basename(hdf5_path))[0]

    ### PLACEMENT ###
    counts = {}
    densities = {}

    for ps in scaffold.get_placement_sets():
        ct = ps.cell_type
        name = ct.name
        n_cells = ps.load_positions().shape[0]
        volume = sum(p.volume() for placement in ct.get_placement() for p in placement.partitions)
        density = (n_cells / volume) if volume > 0 else 0.0
        counts[name] = n_cells
        densities[name] = density

    df_placement = pd.DataFrame({
        "Cell counts": counts,
        "Cell densities [μm⁻³]": densities
    }).sort_index()

    ### CONNECTIVITY ###
    rows = []
    syn_counts = []
    syn_per_pair_mean = []
    conv_mean = []
    div_mean = []

    for cs in scaffold.get_connectivity_sets():
        tag = cs.tag
        conn_set = scaffold.get_connectivity_set(tag).load_connections().as_globals()
        pre_locs, post_locs = conn_set.all()

        if len(pre_locs) == 0 or len(post_locs) == 0:
            continue

        pairs = np.column_stack((pre_locs[:, 0], post_locs[:, 0]))
        unique_pairs, syns_per_pair = np.unique(pairs, axis=0, return_counts=True)

        _, div_counts = np.unique(unique_pairs[:, 0], return_counts=True)
        _, conv_counts = np.unique(unique_pairs[:, 1], return_counts=True)

        rows.append(tag)
        syn_counts.append(len(pre_locs))
        syn_per_pair_mean.append(np.mean(syns_per_pair))
        conv_mean.append(np.mean(conv_counts))
        div_mean.append(np.mean(div_counts))

    df_connectivity = pd.DataFrame({
        "Nb. Synapses": syn_counts,
        "Syn/pair mean": syn_per_pair_mean,
        "Convergence mean": conv_mean,
        "Divergence mean": div_mean,
    }, index=rows)

    ### YAML SUMMARY ###
    summary = {
        "placement": {
            "total_cells": int(df_placement["Cell counts"].sum()),
            "mean_density": float(df_placement["Cell densities [μm⁻³]"].mean()),
            "max_density": float(df_placement["Cell densities [μm⁻³]"].max()),
            "cell_types": {
                k: {
                    "count": int(v),
                    "density": float(df_placement.loc[k, "Cell densities [μm⁻³]"])
                }
                for k, v in counts.items()
            }
        },
        "connectivity": {
            "total_synapses": int(df_connectivity["Nb. Synapses"].sum()),
            "mean_syn_per_pair": float(df_connectivity["Syn/pair mean"].mean()),
            "mean_convergence": float(df_connectivity["Convergence mean"].mean()),
            "mean_divergence": float(df_connectivity["Divergence mean"].mean()),
            "connections": {
                tag: {
                    "synapses": int(df_connectivity.loc[tag, "Nb. Synapses"]),
                    "syn_per_pair": float(df_connectivity.loc[tag, "Syn/pair mean"]),
                    "convergence": float(df_connectivity.loc[tag, "Convergence mean"]),
                    "divergence": float(df_connectivity.loc[tag, "Divergence mean"]),
                }
                for tag in df_connectivity.index
            }
        }
    }

    yaml_path = os.path.join(output_dir, f"{basename}_structure_summary.yaml")
    with open(yaml_path, "w") as f:
        yaml.dump(summary, f, sort_keys=False)

    print(f"✅ YAML salvato in: {yaml_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Estrai riassunto strutturale da un file scaffold HDF5")
    parser.add_argument("--path", required=True, help="Percorso del file .hdf5 dello scaffold")
    args = parser.parse_args()

    extract_structure_summary(args.path)
