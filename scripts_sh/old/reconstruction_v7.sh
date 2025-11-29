#!/bin/bash

timestamp=$(date +"%Y-%m-%d_%H-%M-%S")
main_dir="results_${timestamp}-exp"

mkdir -p "${main_dir}"
echo "📁 Created main directory: ${main_dir}"

backup_dir="${main_dir}/hdf5_backup"
mkdir -p "${backup_dir}"

#rel_values=(1.000 0.957 0.914 0.871 0.829 0.786 0.743 0.700 0.657 0.614 0.571 0.529 0.486 0.443 0.400)
rel_values=(0.400 0.443 0.486 0.529 0.571 0.614 0.657 0.700 0.743 0.786 0.829 0.871 0.914 0.957 1.000)

config_path="configurations/mouse/"

for rel in "${rel_values[@]}"
do
    output_dir="${main_dir}/results_${rel}"
    
    mkdir -p "${output_dir}"
    echo "🔄 Processing ${rel} - Output directory: ${output_dir}"

    yaml_file="mouse_cerebellar_cortex_${rel}.yaml"
    hdf5_file="mouse_cerebellar_cortex_${rel}.hdf5"

    if [ -f "${hdf5_file}" ]; then
        echo "📂 Moving existing ${hdf5_file} to backup directory: ${backup_dir}/"
        mv "${hdf5_file}" "${backup_dir}/"
    fi

    echo "🚀 Running: mpirun -n 8 bsb compile ${config_path}${yaml_file}"
    mpirun -n 8 bsb compile "${config_path}${yaml_file}" || echo "❌ Error compiling ${yaml_file}"

    if [ -f "${hdf5_file}" ]; then
        echo "📁 Moving ${hdf5_file} to ${output_dir}/"
        mv "${hdf5_file}" "${output_dir}/"
    else
        echo "⚠️ Warning: File ${hdf5_file} not found, skipping move."
    fi

    config_used_dir="${main_dir}/configurations_used"
    mkdir -p "${config_used_dir}"
    echo "📄 Copying used YAML files to: ${config_used_dir}/"
    for rel in "${rel_values[@]}"
    do
        yaml_file="mouse_cerebellar_cortex_${rel}.yaml"
        cp "${config_path}${yaml_file}" "${config_used_dir}/" || echo "⚠️ Warning: Could not copy ${yaml_file}"
    done

    echo "✅ Completed processing for rel = ${rel}"
    echo "--------------------------------------------------"
done

echo "🎯 All results stored in: ${main_dir}"

read -p "Press Enter to exit..."
