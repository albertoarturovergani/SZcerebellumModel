#!/bin/bash

# comandi di slurm

# === SETTINGS ===
span_mode="linear" # "linear" oppure "log"
start=60
stop=50
n_steps=1
cores_reconstruction=8
cores_simulations=8
network='cereb_cortex'
stim_name='pawan_stim_vitro' #'pawan_stim_vitro' #'basal_activity' #vbt_stim_protocol_vitro' #'basal_activity' #'kase_stim_vitro' #vbt_stim_protocol_vitro' #'basal_activity' 
name=${stim_name}

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
  stim_name: ${stim_name}
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

global_start=$(date +%s)
for v in "${values[@]}"; do
    ((current_step++))
    echo ""
    echo "=== STEP ${current_step}/${n_steps} ($(timer_eta)) ==="
    v_normalized=$(LC_NUMERIC=C awk -v d="$v" 'BEGIN{printf "%.3f", (100 - d)/100}')
    k_suffix="$v_normalized"  
    log_and_print "Launching morphology generation for degeneration $v% (k=$v_normalized)"

    yaml_file="mouse_cerebellar_cortex_${k_suffix}.yaml"
    config_path="./configurations/mouse/"
    hdf5_base_name="mouse_cerebellar_cortex_${k_suffix}"
    hdf5_file="${hdf5_base_name}.hdf5"
    pdf_orig="bsb_report_structure_${k_suffix}.pdf"
    pdf_orig_f="bsb_report_function_${k_suffix}.pdf"

    index_formatted=$(printf "%03d" "$index")
    v_rounded=$(LC_NUMERIC=C printf "%.2f" "$v")
    int_part=$(echo "$v_rounded" | awk -F. '{printf "%02d", $1}')
    dec_part=$(echo "$v_rounded" | awk -F. '{printf "%02d", $2}')

    recon_dir="${folder_name}/${index_formatted}_deg_${k_suffix}/1.structure"
    config_yaml_dir="${folder_name}/${index_formatted}_deg_${k_suffix}/1.structure"
    stimulus_dir="${folder_name}/${index_formatted}_deg_${k_suffix}/2.stimulus"
    sim_dir="${folder_name}/${index_formatted}_deg_${k_suffix}/3.function"

    mkdir -p "$recon_dir" "$config_yaml_dir" "$stimulus_dir" "$sim_dir"

    # === SEED GENERATION ===
    this_seed=$(date +%s)  # Solo secondi
    echo "  - $this_seed" >> "$global_config_file"

    python ./scripts_py/make_morphologies.py --k "$v_normalized" 
    
    if [ "$stim_name" = "basal_activity" ]; then
        python ./scripts_py/make_network_v2.py --k "$v_normalized" --do_basal
    else
        python ./scripts_py/make_network_v2.py --k "$v_normalized"
    fi

    # === Stimulus generation ===    
    if [[ "$stim_name" == "vbt_stim_protocol_vitro" ]]; then
        python ./scripts_py/make_vbt_stimulus.py --k "$v_normalized" --seed "$this_seed" --stim_name "${stim_name}"
        stim_file="./configurations/mouse/nest/${stim_name}_${k_suffix}.yaml"

    elif [[ "$stim_name" == "kase_stim_vitro" ]]; then
        python ./scripts_py/make_kase_stimulus.py --k "$v_normalized" --seed "$this_seed" --stim_name "${stim_name}"
        stim_file="./configurations/mouse/nest/${stim_name}_${k_suffix}.yaml"

    elif [[ "$stim_name" == "pawan_stim_vitro" ]]; then
        python ./scripts_py/make_pawan_stimulus.py --k "$v_normalized" --seed "$this_seed" --stim_name "${stim_name}"
        stim_file="./configurations/mouse/nest/${stim_name}_${k_suffix}.yaml"
        
    else
        stim_file="./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml"
    fi


    # === Parametri stimolo ===
    #slope="0.000"
    #dt=$(awk '/resolution:/ {print $2}' ./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml)
    #python ./scripts_py/generate_spikes.py --T 5 --dt "$dt" --slope $slope --rate 4 --seed 1234 --datetime $datetime --k_suffix ${k_suffix}
    #python ./scripts_py/replace_background_noise.py ${k_suffix} tempFolders/${datetime}_spikeSlope_${slope}/spike_times_original_${slope}.txt
    #if grep -q "aperiodic:" "./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml"; then
    #  echo "✅ YAML aggiornato con stimolo aperiodico."
    #else
    #  echo -e "${RED}❌ Errore: stimolo 'aperiodic' non presente nel YAML.${NC}"
    #  exit 1
    #fi

    sim_log="${sim_dir}/simulation.log"
    
    # === Compile ===
    compile_log="${recon_dir}/compile.log"
    mpirun -n "$cores_reconstruction" bsb compile ${stim_file}

    # === Simulation ===
    sim_yaml="${config_path}${yaml_file}"
    full_hdf5_path="./${hdf5_file}"

    mpirun -n "$cores_simulations" bsb reconfigure "$full_hdf5_path" "$stim_file" 
    mpirun -n "$cores_simulations" bsb simulate "$full_hdf5_path" "$stim_name" 
    
done

