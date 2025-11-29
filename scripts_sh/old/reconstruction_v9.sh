#!/bin/bash

# Parametri iniziali
hdf5_base_name="mouse_cereb_dcn_io_nest"
yaml_base_name="dcn_io_vitro_nest"
config_path="configurations/mouse/dcn-io/"

rel_values=(0.400 0.443 0.486 0.529 0.571 0.614 0.657 0.700 0.743 0.786 0.829 0.871 0.914 0.957 1.000)

# Directory di output principale con timestamp
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")
main_dir="results_${timestamp}-exp"
mkdir -p "${main_dir}"
echo "📁 Created main directory: ${main_dir}"

# Ricostruzione per ogni rel
for rel in "${rel_values[@]}"; do
    output_dir="${main_dir}/results_${rel}/reconstruction"
    mkdir -p "${output_dir}"
    echo "🔄 Processing rel = ${rel} → ${output_dir}"

    yaml_file="${yaml_base_name}_${rel}.yaml"
    hdf5_file="${hdf5_base_name}_${rel}.hdf5"
    original_pdf="bsb_report_structure.pdf"
    final_pdf="${hdf5_base_name}_${rel}.pdf"

    echo "🚀 Running: mpirun -n 8 bsb compile ${config_path}${yaml_file}"
    mpirun -n 8 bsb compile "${config_path}${yaml_file}" || echo "❌ Error compiling ${yaml_file}"

    if [ -f "${hdf5_file}" ]; then
        echo "📁 Moving ${hdf5_file} to ${output_dir}/"
        mv "${hdf5_file}" "${output_dir}/"
    else
        echo "⚠️ Warning: File ${hdf5_file} not found, skipping move."
    fi

    echo "📄 Copying YAML files ending in _${rel}.yaml from relevant folders..."
    yaml_sources=(
        "${config_path}"
        "configurations/mouse/"
        "configurations/mouse/nest/"
    )

    for src in "${yaml_sources[@]}"; do
        files=("${src}"*_"${rel}".yaml)
        if compgen -G "${src}*_${rel}.yaml" > /dev/null; then
            echo "📂 Found ${#files[@]} file(s) in ${src}, copying..."
            cp "${files[@]}" "${output_dir}/"
        else
            echo "⚠️ No YAML files ending in _${rel}.yaml found in ${src}"
        fi
    done

    if [ -f "${original_pdf}" ]; then
        echo "📄 Moving ${original_pdf} to ${output_dir}/"
        mv "${original_pdf}" "${output_dir}/"
        echo "✏️ Renaming PDF to ${final_pdf}"
        mv "${output_dir}/${original_pdf}" "${output_dir}/${final_pdf}"
    else
        echo "⚠️ Warning: File ${original_pdf} not found, skipping move."
    fi

    echo "✅ Completed processing for rel = ${rel}"
    echo "--------------------------------------------------"
done

echo "🎯 All results stored in: ${main_dir}"
read -p "Press Enter to exit..."
