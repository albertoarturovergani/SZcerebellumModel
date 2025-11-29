#!/bin/bash

# === SETTINGS ===
name="vbt_test"
span_mode="linear" # "linear" oppure "log"
start=60
stop=100
n_steps=1
cores_reconstruction=8
cores_simulations=8
network='cereb_cortex'

# === COLORS ===
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# === CHECK DEPENDENCIES ===
for cmd in bsb python mpirun awk date mkdir tr tee; do
    if ! command -v $cmd &>/dev/null; then
        echo -e "${RED}❌ Required command '$cmd' not found, aborting.${NC}"
        exit 1
    fi
done

# === FUNCTIONS ===
log_and_print() {
    local datetime=$(date +"%Y%m%d%H%M")
    echo -e "[$datetime] $1" | tee -a "$log_file"
}

safe_mkdir() {
    mkdir -p "$1" || { echo -e "${RED}[ERROR] Failed to create directory $1${NC}"; exit 1; }
}

timer_start() {
    start_time=$(date +%s)
}

timer_eta() {
    local now=$(date +%s)
    local elapsed=$(( now - start_time ))
    local remaining=$(( (elapsed / current_step) * (n_steps - current_step) ))
    printf "%02d:%02d remaining" $((remaining/60)) $((remaining%60))
}

# === INIT FOLDERS ===
datetime=$(date +"%Y%m%d%H%M%S")
folder_name="${datetime}_${name}"
safe_mkdir "$folder_name"

log_file="${folder_name}/log.txt"
touch "$log_file"
log_and_print "Created folder '$folder_name' and log."

# === CREATE BASE FOLDERS ===
global_config_file="${folder_name}/info.yaml"
mkdir -p "$folder_name/A_replica" "$folder_name/B_collective_results/1.structures" "$folder_name/B_collective_results/2.functions"

cat > "$global_config_file" <<EOF
experiment:
  name: ${name}
  date: ${datetime}
  user: ${USER:-unknown}
  span_mode: ${span_mode}
  start: ${start}
  stop: ${stop}
  n_steps: ${n_steps}
  cores_reconstruction: ${cores_reconstruction}
  cores_simulations: ${cores_simulations}
  network: ${network}
seeds:
EOF

# === GENERATE SPAN VALUES ===
values=()
seeds=()

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
    log_and_print "${RED}[ERROR] Invalid span mode: $span_mode${NC}"
    exit 1
fi

log_and_print "Generated span values: ${values[*]}"

# === DEGENERATION SWEEP ===
index=0
current_step=0
timer_start

for v in "${values[@]}"; do
    ((current_step++))
    echo ""
    echo "=== STEP ${current_step}/${n_steps} ($(timer_eta)) ==="
    v_normalized=$(LC_NUMERIC=C awk -v d="$v" 'BEGIN{printf "%.3f", (100 - d)/100}')
    k_suffix="$v_normalized"  
    log_and_print "Launching morphology generation for degeneration $v% (k=$v_normalized)"

    # === SEED GENERATION ===
    this_seed=$(date +%s%N | cut -b1-13)
    echo "  - $this_seed" >> "$global_config_file"

    python ./scripts_py/make_morphologies.py --k "$v_normalized" 
    python ./scripts_py/make_network_v2.py --k "$v_normalized"

    yaml_file="mouse_cerebellar_cortex_${k_suffix}.yaml"
    config_path="./configurations/mouse/"
    hdf5_base_name="mouse_cerebellar_cortex_${k_suffix}"
    hdf5_file="${hdf5_base_name}.hdf5"
    pdf_orig="bsb_report_structure_${k_suffix}.pdf"

    index_formatted=$(printf "%03d" "$index")
    v_rounded=$(LC_NUMERIC=C printf "%.2f" "$v")
    int_part=$(echo "$v_rounded" | awk -F. '{printf "%02d", $1}')
    dec_part=$(echo "$v_rounded" | awk -F. '{printf "%02d", $2}')
    val_formatted="${int_part}_${dec_part}"

    recon_dir="${folder_name}/${index_formatted}_deg_${val_formatted}/3.structure"
    config_yaml_dir="${folder_name}/${index_formatted}_deg_${val_formatted}/1.config_files/"
    stimulus_dir="${folder_name}/${index_formatted}_deg_${val_formatted}/2.stimulus"
    sim_dir="${folder_name}/${index_formatted}_deg_${val_formatted}/4.function"

    mkdir -p "$recon_dir" "$config_yaml_dir" "$stimulus_dir" "$sim_dir"

    # === Compile ===
    compile_log="${recon_dir}/compile.log"
    if ! mpirun -n "$cores_reconstruction" bsb compile "${config_path}${yaml_file}" 2>&1 | tee "$compile_log"; then
        log_and_print "${RED}❌ Error compiling ${yaml_file}, see ${compile_log}${NC}"
        ((index++))
        continue
    fi

    # === Simulation ===
    sim_yaml="${config_path}${yaml_file}"
    full_hdf5_path="./${hdf5_base_name}.hdf5"
    stim_file="./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml"
    sim_name="basal_activity"
    sim_log="${sim_dir}/simulation.log"

    mpirun -n "$cores_simulations" bsb reconfigure "$full_hdf5_path" "$stim_file" 
    mpirun -n "$cores_simulations" bsb simulate "$full_hdf5_path" "$sim_name" 

    echo -e "${YELLOW}📊 Generating dynamics PDF...${NC}"
    python ./scripts_py/report_dynamics_v4.py "$full_hdf5_path" "$sim_name" "$sim_dir" || echo -e "${RED}❌ PDF generation failed.${NC}"

    # === Move Files ===
    for nio_file in *.nio; do
        [ -f "$nio_file" ] && mv "$nio_file" "$sim_dir/"
    done

    [ -f "${config_path}${yaml_file}" ] && mv "${config_path}${yaml_file}" "$config_yaml_dir/"
    [ -f "${config_path}morphologies_${k_suffix}.yaml" ] && mv "${config_path}morphologies_${k_suffix}.yaml" "$config_yaml_dir/"
    [ -f "./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml" ] && mv "./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml" "$stimulus_dir/"

    [ -f "$hdf5_file" ] && mv "$hdf5_file" "$recon_dir/"
    [ -f "$pdf_orig" ] && mv "$pdf_orig" "${recon_dir}/"

    morphologies_path="./morphologies/"
    mkdir -p "${recon_dir}/morphologies/"
    for cell_type in GranuleCell GolgiCell PurkinjeCell StellateCell BasketCell; do
        [ -d "${morphologies_path}${cell_type}" ] && mv "${morphologies_path}${cell_type}" "${recon_dir}/morphologies/"
    done
    [ -f "${morphologies_path}morphometrics_${k_suffix}.csv" ] && mv "${morphologies_path}morphometrics_${k_suffix}.csv" "${recon_dir}/morphologies/"
    [ -f "${morphologies_path}morphologies_${k_suffix}.yaml" ] && mv "${morphologies_path}morphologies_${k_suffix}.yaml" "${recon_dir}/morphologies/"

    ((index++))
done

log_and_print "${GREEN}✅ Finished all reconstructions and simulations.${NC}"

# === FINAL ANALYSIS ===
python ./scripts_py/analyze_firing_v4.py "$folder_name"
python ./scripts_py/overallAnalysis_v3.py "$folder_name"

log_and_print "${GREEN}✅ Full pipeline completed.${NC}"
read -s -n 1 -p "Press any key to exit..."
