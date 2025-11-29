#!/bin/bash

timestamp='2025-04-08_16-11-43' # Update manually if needed
main_dir="results_${timestamp}-exp"

simulation_types=("basal_activity") #"mf_cf_stimulus")

hdf5_files=$(find "${main_dir}" -type f -name "*.hdf5")

if [ -z "$hdf5_files" ]; then
    echo "❌ No previous HDF5 files found in ${main_dir}. Exiting."
    exit 1
fi

for hdf5_file in $hdf5_files; do
    result_dir=$(dirname "$hdf5_file")
    rel=$(basename "$hdf5_file" | grep -oP '\d+\.\d+')

    echo "🔹 Processing rel = ${rel} (HDF5: $(basename "$hdf5_file"))"

    for sim_type in "${simulation_types[@]}"; do
        sim_dir="${result_dir}/${sim_type}"  
        mkdir -p "${sim_dir}"
        
        # 🔹 Reconfigure simulation
        echo "⚙️  Reconfiguring with bsb..."
        mpirun -n 8 bsb reconfigure "${hdf5_file}" "./configurations/mouse/nest/basal_vitro_${rel}.yaml" || echo "❌ Reconfiguration failed"

        echo "🚀 Running ${sim_type} simulation for rel = ${rel}"
        mpirun -n 8 bsb simulate "${hdf5_file}" ${sim_type} || echo "❌ Simulation failed for ${rel} - ${sim_type}"
        shopt -s nullglob
        for nio_file in *.nio; do
            mv "$nio_file" "${sim_dir}/"
        done

        #echo "📊 Generating dynamics report..."
        #ipython ./report_dynamics_v2.py "${hdf5_file}" ${sim_type} "${sim_dir}" || echo "❌ Dynamics report failed for ${rel} - ${sim_type}"

        echo "✅ Completed ${sim_type} simulation for rel = ${rel}"
        echo "--------------------------------------------------"
    done
done

echo "🎯 All results stored in: ${main_dir}"

# Prevents terminal from closing immediately
read -s -n 1 -p "Press any key to exit..."
