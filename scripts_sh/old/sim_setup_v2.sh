#!/bin/bash

# comandi di slurm

# === SETTINGS ===
span_mode="linear" # "linear" oppure "log"
start=0
stop=60 # no more than 50 not bioplausible
n_steps=10
cores_reconstruction=8
cores_simulations=8
network='cereb_cortex'
stim_name='pawan_stim_vitro' #'vbt_stim_protocol_vitro' #kase_stim_vitro' #'pawan_stim_vitro' #'basal_activity' #'pawan_stim_vitro' #'basal_activity'#'pawan_stim_vitro' #'basal_activity' #vbt_stim_protocol_vitro' #'basal_activity' #'kase_stim_vitro' #vbt_stim_protocol_vitro' #'basal_activity'

# === MODELING PARAMETERS ===
reserve_k0=0.75
compensation=1         # se non specificato, default = 1
reserve_sharpness=50 #0.001 #0.001  # se non specificato, default = 0.001
gaba=0
# === SIGMA (solo per stimoli non-basali) ===
if [[ "$stim_name" == "basal_activity" ]]; then
    sigmafr=100
else
    sigmafr=10
fi
note=${reserve_sharpness}_${compensation}_${gaba}
name=${stim_name}_${note}

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
modeling_parameters:
  reserve_k0: ${reserve_k0}
  reserve_sharpness: ${reserve_sharpness}
  compensation: ${compensation}
  gaba: ${gaba}
  sigmafr: ${sigmafr}

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
        #python ./scripts_py/make_network_v3.py --k "$v_normalized" --do_basal --reserve "$reserve_k0" "$reserve_sharpness" --compensation "$compensation"
        python ./scripts_py/make_network_v4.py --k "$v_normalized" --do_basal --reserve "$reserve_k0" "$reserve_sharpness" --compensation "$compensation" --gaba "$gaba"

    else
        python ./scripts_py/make_network_v4.py --k "$v_normalized" --reserve "$reserve_k0" "$reserve_sharpness" --compensation "$compensation" --gaba "$gaba"
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
    
    # === Move Files ===
    for nio_file in *.nio; do
        [ -f "$nio_file" ] && mv "$nio_file" "$sim_dir/"
    done

    echo -e "${YELLOW}📊 Generating dynamics PDF...${NC}"
    python ./scripts_py/report_dynamics_v4.py "$full_hdf5_path" "$stim_name" "$sim_dir" || echo -e "${RED}❌ PDF generation failed.${NC}"

    sleep 5
    
    [ -f "${config_path}${yaml_file}" ] && mv "${config_path}${yaml_file}" "$config_yaml_dir/"
    [ -f "${config_path}morphologies_${k_suffix}.yaml" ] && mv "${config_path}morphologies_${k_suffix}.yaml" "$config_yaml_dir/"
    # === Move stimulus YAML into stimulus_dir ===
    if [[ "$stim_name" != "basal_activity" ]]; then
        if [ -f "./configurations/mouse/nest/${stim_name}_${k_suffix}.yaml" ]; then
            mv "./configurations/mouse/nest/${stim_name}_${k_suffix}.yaml" "$stimulus_dir/"
        fi
    else
        if [ -f "./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml" ]; then
            mv "./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml" "$stimulus_dir/"
        fi
    fi

    #python ./scripts_py/plot_stim_profile.py "${stim_file}" --savePath "${stimulus_dir}"
    python ./scripts_py/generate_stimulus_timeline.py \
      --stimDir "$stimulus_dir/*.yaml" \
      --stimName "${stim_file}" \
      --outPath "${stimulus_dir}"
    
    [ -f "$hdf5_file" ] && mv "$hdf5_file" "$recon_dir/"
    [ -f "$pdf_orig" ] && mv "$pdf_orig" "${recon_dir}/"
    [ -f "$pdf_orig_f" ] && mv "$pdf_orig" "${recon_dir}/"

    echo -e "${YELLOW}🧠 Generating structure YAML summary...${NC}"
    python ./scripts_py/do_structure_yaml.py --path "${recon_dir}/${hdf5_file}"

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
# STRUCTURE
echo "➡️  Running morphology_analysis.py ..."
python ./scripts_py/do_morphology_analysis.py --path "$folder_name"
echo "➡️  Running do_network_analysis.py ..."
python ./scripts_py/do_network_analysis.py --path "$folder_name"
echo "➡️  Running analyze_firing_v4.py ..."
# FUNCTIONS
python ./scripts_py/analyze_firing_v5.py --path "$folder_name" --sigma "$sigmafr"
python ./scripts_py/do_completeAnalysis_v2.py --path "$folder_name" --sigma "$sigmafr"
python ./scripts_py/overallAnalysis_v3.py "$folder_name"


if [[ "$stim_name" == "pawan_stim_vitro" ]]; then
  echo "➡️  Running analyze_multimeter_v4.py ..."
  python ./scripts_py/analyze_multimeter_v4.py --path "$folder_name"
fi

if [[ "$stim_name" == "vbt_stim_protocol_vitro" || "$stim_name" == "kase_stim_vitro" ]]; then
  echo "➡️  Running overallAnalysis_v3.py ..."
  for sigma in 10 100; do
    echo "▶️  Sigma = $sigma"
    python ./scripts_py/do_completeAnalysis_vbtkase.py --path "$folder_name" --sigma "$sigma" --stimName "$stim_name"
  done
fi


global_end=$(date +%s)
global_total_time=$((global_end - global_start))

echo "" >> "$global_config_file"
echo "global_timings:" >> "$global_config_file"
echo "  total_time_seconds: $global_total_time" >> "$global_config_file"

echo "➡️  Merging structural and functional summaries (PKL) ..."

python3 <<EOF
import pandas as pd
from pathlib import Path

folder = Path("$folder_name")
summary_struct = folder / "B_collective_results/1.structures/df_summary.pkl"
summary_func = folder / "B_collective_results/2.functions/df_summary.pkl"
output_merged = folder / "B_collective_results/df_summary_merged.pkl"

if summary_struct.exists() and summary_func.exists():
    df1 = pd.read_pickle(summary_struct)
    df2 = pd.read_pickle(summary_func)

    if 'k' not in df1.columns or 'k' not in df2.columns:
        print("❌ One of the summaries does not contain the 'k' column.")
    else:
        df_merged = pd.merge(df1, df2, on='k', suffixes=('_struct', '_func'))
        df_merged.to_pickle(output_merged)
        print(f"✅ Merged summary saved to: {output_merged}")
else:
    print("❌ One or both df_summary.pkl files not found.")
EOF


log_and_print "${GREEN}✅ Full pipeline completed.${NC}"
read -s -n 1 -p "Press any key to exit..."

