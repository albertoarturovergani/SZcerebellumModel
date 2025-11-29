#!/bin/bash

# 📌 Timestamp per organizzare le cartelle
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")
main_dir="results_${timestamp}-exp"
config_path="configurations/mouse/dcn-io/"
backup_dir="${main_dir}/hdf5_backup"

# 📌 Sweep range (Atrophy levels)
rel_values=(0.400 0.443 0.486 0.529 0.571 0.614 0.657 0.700 0.743 0.786 0.829 0.871 0.914 0.957 1.000)
simulation_types=("basal_activity")
fileName="mouse_cereb_dcn_io_nest"

# 📌 Creazione cartelle principali
mkdir -p "${main_dir}" "${backup_dir}"
echo "📁 Created main directory: ${main_dir}"

# 📌 Loop principale su ogni livello di atrofia
for rel in "${rel_values[@]}"
do
    output_dir="${main_dir}/results_${rel}"
    mkdir -p "${output_dir}"
    
    yaml_file="dcn_io_vitro_nest_${rel}.yaml"
    hdf5_file="${fileName}_${rel}.hdf5"

    echo "🔄 Processing ${rel} - Output directory: ${output_dir}"

    # 🔹 Backup HDF5 esistente
    if [ -f "${hdf5_file}" ]; then
        echo "📂 Moving existing ${hdf5_file} to backup directory"
        mv "${hdf5_file}" "${backup_dir}/"
    fi

    # 🔹 Compila la rete
    echo "🚀 Compiling ${yaml_file}"
    mpirun -n 8 bsb compile "${config_path}${yaml_file}" || echo "❌ Compilation failed for ${yaml_file}"

    # 🔹 Sposta il file HDF5 compilato
    if [ -f "${hdf5_file}" ]; then
        mv "${hdf5_file}" "${output_dir}/"
    else
        echo "⚠️ Warning: File ${hdf5_file} not found, skipping move."
    fi

    # 🔹 Esegui le simulazioni
    for sim_type in "${simulation_types[@]}"; do
        sim_dir="${output_dir}/${sim_type}"
        mkdir -p "${sim_dir}"

        echo "⚙️ Reconfiguring & Running ${sim_type} simulation for rel=${rel}"
        mpirun -n 8 bsb reconfigure "${output_dir}/${hdf5_file}" "${config_path}${yaml_file}" || echo "❌ Reconfiguration failed"
        mpirun -n 8 bsb simulate "${output_dir}/${hdf5_file}" ${sim_type} || echo "❌ Simulation failed"

        # 🔹 Sposta i file di output
        for nio_file in *.nio; do
            mv "$nio_file" "${sim_dir}/"
        done

        echo "✅ Completed ${sim_type} simulation for rel = ${rel}"
    done

    echo "--------------------------------------------------"
done

echo "🎯 All results stored in: ${main_dir}"
python saveDynamics_v5.py "${main_dir}"

