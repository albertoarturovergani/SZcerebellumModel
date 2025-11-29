#!/bin/bash

# === SETTINGS ===
name="vbt_test"
span_mode="linear" # "linear" oppure "log"
start=60
stop=70
n_steps=1
cores_reconstruction=8
cores_simulations=8
network_type="cortex" # "dcn_io" oppure "cortex"

# === FUNCTIONS ===
log_and_print() {
    local datetime=$(date +"%Y%m%d%H%M")
    echo "[$datetime] $1" | tee -a "$log_file"
}

safe_mkdir() {
    mkdir -p "$1" || { echo "[ERROR] Failed to create directory $1"; exit 1; }
}

# === INIT FOLDERS ===
datetime=$(date +"%Y%m%d%H%M%S")
folder_name="${datetime}_${name}"
safe_mkdir "$folder_name"

log_file="${folder_name}/log.txt"
touch "$log_file"
log_and_print "Created folder '$folder_name' and log."

# === GENERATE SPAN VALUES ===
values=()
if [[ "$span_mode" == "linear" ]]; then
    step=$(awk -v s="$start" -v e="$stop" -v n="$n_steps" 'BEGIN{print (n==1)?0:(e-s)/(n-1)}')
    for ((i=0; i<n_steps; i++)); do
        val=$(awk -v s="$start" -v st="$step" -v i="$i" 'BEGIN{printf "%.6f", s+st*i}' | tr ',' '.')
        values+=("$val")
    done
elif [[ "$span_mode" == "log" ]]; then
    log_start=$(awk -v s="$start" 'BEGIN{print log(s+1)}')
    log_stop=$(awk -v s="$stop" 'BEGIN{print log(s+1)}')
    for ((i=0; i<n_steps; i++)); do
        fraction=$(awk -v i="$i" -v n="$n_steps" 'BEGIN{print (n==1)?0:i/(n-1)}')
        log_val=$(awk -v ls="$log_start" -v le="$log_stop" -v f="$fraction" 'BEGIN{print ls+(le-ls)*f}')
        val=$(awk -v lv="$log_val" 'BEGIN{printf "%.6f", exp(lv)-1}' | tr ',' '.')
        values+=("$val")
    done
else
    log_and_print "[ERROR] Invalid span mode: $span_mode"
    exit 1
fi
log_and_print "Generated span values."

# === CREATE BASE FOLDERS ===
global_config_file="${folder_name}/info.yaml"
mkdir -p "$folder_name/C_seeds" "$folder_name/A_replica" "$folder_name/B_collective_results/1.structures" "$folder_name/B_collective_results/2.functions"

echo "experiment:" > "$global_config_file"
echo "  name: ${name}" >> "$global_config_file"
echo "  date: ${datetime}" >> "$global_config_file"
echo "  user: ${USER:-unknown}" >> "$global_config_file"
echo "  span_mode: ${span_mode}" >> "$global_config_file"
echo "  start: ${start}" >> "$global_config_file"
echo "  stop: ${stop}" >> "$global_config_file"
echo "  n_steps: ${n_steps}" >> "$global_config_file"
echo "  cores_reconstruction: ${cores_reconstruction}" >> "$global_config_file"
echo "  cores_simulations: ${cores_simulations}" >> "$global_config_file"
echo "deg_sweep:" >> "$global_config_file"

# === PREPARE DEGENERATION SWEEP ===
index=0
base_seed=$(date +"%s")
for v in "${values[@]}"; do
    v_normalized=$(LC_NUMERIC=C awk -v d="$v" 'BEGIN{printf "%.3f", (100 - d)/100}')
    k_suffix="$v_normalized"  
    log_and_print "Launching morphology generation for degeneration $v% (k=$v_normalized)"

    python ./scripts_py/make_morphologies.py --k "$v_normalized"
    python ./scripts_py/make_network.py --k "$v_normalized"

    echo "🚧 Compiling network..."
    yaml_file="mouse_cerebellar_cortex_${k_suffix}.yaml"
    config_path="/home/alberto/cerebellum/cerebellum/configurations/mouse/"
    hdf5_base_name="mouse_cerebellar_cortex_${k_suffix}"
    hdf5_file="${hdf5_base_name}.hdf5"
    pdf_orig="bsb_report_structure.pdf"

    # Define destination folders
    index_formatted=$(printf "%03d" "$index")
    v_rounded=$(LC_NUMERIC=C printf "%.2f" "$v")
    int_part=$(echo "$v_rounded" | awk -F. '{printf "%02d", $1}')
    dec_part=$(echo "$v_rounded" | awk -F. '{printf "%02d", $2}')
    val_formatted="${int_part}_${dec_part}"

    recon_dir="${folder_name}/${index_formatted}_deg_${val_formatted}/3.structure"
    config_yaml_dir="${folder_name}/${index_formatted}_deg_${val_formatted}/1.config_files/config_yaml"
    stimulus_dir="${folder_name}/${index_formatted}_deg_${val_formatted}/2.stimulus"
    sim_dir="${folder_name}/${index_formatted}_deg_${val_formatted}/4.function"

    mkdir -p "$recon_dir" "$config_yaml_dir" "$stimulus_dir" "$sim_dir"

    # === Compile ===
    if ! mpirun -n 8 bsb compile "${config_path}${yaml_file}"; then
        echo "❌ Error compiling ${yaml_file}, skipping degeneration = ${v}%"
        ((index++))
        continue
    fi

    # === Move compiled HDF5 ===
    if [ -f "${hdf5_file}" ]; then
        mv "${hdf5_file}" "${recon_dir}/"
        echo "✅ Moved ${hdf5_file} to ${recon_dir}/"
    else
        echo "⚠️ No ${hdf5_file} found after compilation."
    fi

    # === Move Structure PDF ===
    if [ -f "${pdf_orig}" ]; then
        pdf_target="${hdf5_base_name}.pdf"
        mv "${pdf_orig}" "${recon_dir}/${pdf_target}"
        echo "✅ Moved and renamed PDF to ${recon_dir}/${pdf_target}"
    else
        echo "⚠️ No ${pdf_orig} found after compilation."
    fi

    # === Move Morphology Folders ===
    morphologies_path="./morphologies/"
    mkdir -p "${recon_dir}/morphologies/"
    for cell_type in GranuleCell GolgiCell PurkinjeCell StellateCell BasketCell; do
        cell_folder="${morphologies_path}${cell_type}"
        if [ -d "$cell_folder" ]; then
            mv "$cell_folder" "${recon_dir}/morphologies/"
            echo "✅ Moved full folder ${cell_type} to ${recon_dir}/morphologies/"
        else
            echo "⚠️ No folder found for ${cell_type}."
        fi
    done

    # === Move Morphometrics CSV ===
    morphometrics_csv="morphometrics_${k_suffix}.csv"
    if [ -f "${morphologies_path}${morphometrics_csv}" ]; then
        mv "${morphologies_path}${morphometrics_csv}" "${recon_dir}/morphologies/"
        echo "✅ Moved ${morphometrics_csv} to ${recon_dir}/morphologies/"
    else
        echo "⚠️ No morphometrics CSV for k=${k_suffix}"
    fi

    # === Move Morphologies YAML ===
    morphologies_yaml="${morphologies_path}morphologies_${k_suffix}.yaml"
    if [ -f "$morphologies_yaml" ]; then
        mv "$morphologies_yaml" "${recon_dir}/morphologies/"
        echo "✅ Moved ${morphologies_yaml} to ${recon_dir}/morphologies/"
    else
        echo "⚠️ No morphologies YAML for k=${k_suffix}"
    fi

    # === Simulazione ===
    sim_yaml="${config_path}${yaml_file}"  # Ancora nel path originale!
    full_hdf5_path="${recon_dir}/${hdf5_base_name}.hdf5"
    stim_file="./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml"  # Ancora nel path originale!
    sim_name="basal_activity"

    if [ -f "$sim_yaml" ]; then
        echo "⚙️ Reconfiguring network..."
        mpirun -n 8 bsb reconfigure "$full_hdf5_path" "$stim_file" || {
            echo "❌ Reconfigure failed"
            ((index++))
            continue
        }

        echo "🚀 Running simulation..."
        mpirun -n 8 bsb simulate "$full_hdf5_path" "$sim_name" || {
            echo "❌ Simulation failed"
            ((index++))
            continue
        }
    else
        echo "❌ Simulation YAML not found: ${sim_yaml}, skipping"
        ((index++))
        continue
    fi

    # === Move NIO Files ===
    for nio_file in *.nio; do
        if [ -f "$nio_file" ]; then
            mv "$nio_file" "${sim_dir}/"
            echo "✅ Moved ${nio_file} to ${sim_dir}/"
        fi
    done

    # === Move YAMLs After Sim ===
    echo "📦 Moving YAML files into correct folders after simulation..."
    for yaml_file in \
        "mouse_cerebellar_cortex_${k_suffix}.yaml" \
        "morphologies_${k_suffix}.yaml" \
        "basal_vitro_${k_suffix}.yaml" \
        "dcn_${k_suffix}.yaml" \
        "dcn_io_${k_suffix}.yaml" \
        "dcn_vitro_nest_${k_suffix}.yaml" \
        "dcn_io_vitro_nest_${k_suffix}.yaml"; do

        if [ -f "./configurations/mouse/${yaml_file}" ]; then
            src="./configurations/mouse/${yaml_file}"
        elif [ -f "./configurations/mouse/dcn-io/${yaml_file}" ]; then
            src="./configurations/mouse/dcn-io/${yaml_file}"
        elif [ -f "./configurations/mouse/nest/${yaml_file}" ]; then
            src="./configurations/mouse/nest/${yaml_file}"
        else
            echo "⚠️ Warning: ${yaml_file} not found!"
            continue
        fi

        if [[ "$yaml_file" == *"vitro"* ]]; then
            dest="$stimulus_dir"
        else
            dest="$config_yaml_dir"
        fi

        mv "$src" "$dest/"
        echo "✅ Moved ${yaml_file} to ${dest}/"
    done

    # === Generate Dynamics PDF ===
    echo "📊 Generating dynamics PDF..."
    python ./scripts_py/report_dynamics_v3.py "$full_hdf5_path" "$sim_name" "$sim_dir" || echo "❌ PDF generation failed"

    # === Update index ===
    ((index++))
done

log_and_print "✅ Finished morphology, network compilation, simulation, YAML organization."

# === FINAL ANALYSIS ===
python ./scripts_py/analyze_firing_v4.py "$folder_name"
python ./scripts_py/overallAnalysis_v3.py "$folder_name"

log_and_print "✅ Finished full pipeline."
read -s -n 1 -p "Press any key to exit..."


