#!/bin/bash

# Parametri iniziali
hdf5_base_name="mouse_cerebellar_cortex" #"mouse_cereb_dcn_io_nest"
yaml_base_name="mouse_cerebellar_cortex" #"dcn_io_vitro_nest"
config_path="configurations/mouse/" #"configurations/mouse/dcn-io/"
stim_fil="basal_vitro"
stimFolder="nest" #dcn-io
sim_name="basal_activity"
#rel_values=(0.400 0.443 0.486 0.529 0.571 0.614 0.657 0.700 0.743 0.786 0.829 0.871 0.914 0.957 1.000)
#rel_values=(3.000 2.733 2.467 2.200 1.933 1.667 1.400 1.133 0.867 0.600)

# Timestamp e directory principale
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")
main_dir="results_${timestamp}-exp"
mkdir -p "${main_dir}"
echo "📁 Created main directory: ${main_dir}"

for rel in "${rel_values[@]}"; do
    echo "🔄 Processing rel = ${rel}"
    
    yaml_file="${yaml_base_name}_${rel}.yaml"
    hdf5_file="${hdf5_base_name}_${rel}.hdf5"
    output_dir="${main_dir}/results_${rel}"
    recon_dir="${output_dir}/reconstruction"
    sim_dir="${output_dir}/${sim_name}"

    mkdir -p "${recon_dir}"
    mkdir -p "${sim_dir}"

    echo "🚧 Compiling network..."
    mpirun -n 8 bsb compile "${config_path}${yaml_file}" || {
        echo "❌ Error compiling ${yaml_file}, skipping rel = ${rel}"
        continue
    }

    if [ -f "${hdf5_file}" ]; then
        mv "${hdf5_file}" "${recon_dir}/"
    else
        echo "⚠️ HDF5 file not found, skipping rel = ${rel}"
        continue
    fi

    echo "📁 Copying YAMLs for rel = ${rel}"
    yaml_sources=(
        "${config_path}"
        "configurations/mouse/"
        "configurations/mouse/nest/"
    )
    for src in "${yaml_sources[@]}"; do
        files=("${src}"*_"${rel}.yaml")
        if compgen -G "${src}*_${rel}.yaml" > /dev/null; then
            cp "${files[@]}" "${recon_dir}/"
        fi
    done

    pdf_orig="bsb_report_structure.pdf"
    pdf_target="${hdf5_base_name}_${rel}.pdf"
    if [ -f "${pdf_orig}" ]; then
        mv "${pdf_orig}" "${recon_dir}/"
        mv "${recon_dir}/${pdf_orig}" "${recon_dir}/${pdf_target}"
    fi

    echo "🧠 Plotting stimulus timeline..."
    #python3 plot_stim_timeline_v5.py "${recon_dir}"
    python3 plot_stim_timeline_v4.py "${recon_dir}"
    #python3 plot_stim_timeline_v3.py "${recon_dir}"

    echo "⚙️ Reconfiguring and simulating..."
    sim_yaml="./${recon_dir}/${yaml_file}"
    full_hdf5_path="${recon_dir}/${hdf5_file}"
    stim_file="./${config_path}/${stimFolder}/${stim_fil}_${rel}.yaml"
    
    if [ -f "${sim_yaml}" ]; then
        mpirun -n 8 bsb reconfigure "${full_hdf5_path}" "${stim_file}" || echo "❌ Reconfigure failed"
        mpirun -n 8 bsb simulate "${full_hdf5_path}" "${sim_name}" || echo "❌ Simulation failed"
    else
        echo "❌ Simulation YAML not found: ${sim_yaml}, skipping"
        continue
    fi

    for nio_file in *.nio; do
        mv "$nio_file" "${sim_dir}/"
    done

    echo "📊 Generating dynamics PDF..."
    python report_dynamics_v3.py "${full_hdf5_path}" "${sim_name}" "${sim_dir}" || echo "❌ PDF generation failed"
    #python plotSimSeq.py "${main_dir}"
    


    echo "✅ Done for rel = ${rel}"
    echo "--------------------------------------------------"
done

echo "🎯 All simulations and results saved in: ${main_dir}"
python analyze_firing_v2.py "${main_dir}"
python overallAnalysis_v2.py "${main_dir}"

read -s -n 1 -p "Press any key to exit..."
