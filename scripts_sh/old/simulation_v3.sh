#!/bin/bash

timestamp='2025-04-15_10-26-52'  # aggiorna se necessario
main_dir="results_${timestamp}-exp"

# Nome della simulazione (BSB deve trovarlo nello YAML compilato)
sim_name="basal_activity"

# YAML stimolo: parte fissa del nome
stim_prefix="basal_vitro"

hdf5_files=$(find "${main_dir}" -type f -path "*/reconstruction/*.hdf5")

if [ -z "$hdf5_files" ]; then
    echo "❌ No reconstruction HDF5 files found in ${main_dir}."
    exit 1
fi

for hdf5_file in $hdf5_files; do
    recon_dir=$(dirname "$hdf5_file")
    rel_dir=$(dirname "$recon_dir")
    rel=$(basename "$rel_dir" | grep -oP '\d+\.\d+')

    echo "🔹 Running simulation for rel = ${rel}"

    #stim_file="./configurations/mouse/nest/${stim_prefix}_${rel}.yaml"
    stim_file="./configurations/mouse/dcn-io/dcn_io_vitro_nest_${rel}.yaml"
    
    sim_dir="${rel_dir}/${sim_name}"
    mkdir -p "${sim_dir}"

    echo "⚙️ Reconfiguring with stimulus file: ${stim_file}"
    if [ ! -f "${stim_file}" ]; then
        echo "❌ Stimulus YAML ${stim_file} not found. Skipping rel = ${rel}"
        continue
    fi
    mpirun -n 8 bsb reconfigure "${hdf5_file}" "${stim_file}" || echo "❌ Reconfiguration failed"

    echo "🚀 Running simulation: ${sim_name}"
    mpirun -n 8 bsb simulate "${hdf5_file}" "${sim_name}" || echo "❌ Simulation failed"

    shopt -s nullglob
    for nio_file in *.nio; do
        mv "$nio_file" "${sim_dir}/"
    done


    echo "📊 Generating dynamics PDF..."
    python report_dynamics_v3.py "${hdf5_file}" "${sim_name}" "${sim_dir}" || echo "❌ PDF generation failed"
    python plotSimSeq.py "${main_dir}"

    echo "✅ Completed simulation ${sim_name} for rel = ${rel}"
    echo "--------------------------------------------------"
done

echo "🎯 All simulations saved in: ${main_dir}"
read -s -n 1 -p "Press any key to exit..."
