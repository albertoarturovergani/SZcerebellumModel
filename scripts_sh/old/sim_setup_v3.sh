#!/bin/bash

# comandi di slurm
# bash sim_setup_v3.sh "$stim" "$start" "$stop" "$n_steps" "$comp" "$gaba" "$do_dcnio" "$reserve_k0" "$reserve_sharpness"
# === SETTINGS === (parametri dinamici da terminale)

stim_name="${1:-basal_activity}"   # es. basal_activity, kase_stim_vitro, ...
start="${2:-0}"                    # in %
stop="${3:-50}"                    # in %
n_steps="${4:-10}"                 # numero steps
comp_state="${5:-1}"  # 0 o 1
compensation_str="{'state': ${comp_state}, 'trasmitter': ['gaba','glutamate'], 'neurons': 'all'}"
sz_hyph_raw="${6:-gaba:all}"       # es. gaba:all
do_dcnio="${7:-1}"  # 0 o 1
reserve_k0="${8:-0.5}"
reserve_sharpness="${9:-40}"
reserve_str="{'k0': ${reserve_k0}, 'sharpness': ${reserve_sharpness}}"

# Costanti
span_mode="linear"                # oppure "log"
cores_reconstruction=8
cores_simulations=8
network='cereb_cortex'

# Parsing sz_hyph_raw
if [[ "$sz_hyph_raw" == *:* ]]; then
  IFS=':' read -r sz_transmitter sz_neurons <<< "$sz_hyph_raw"
else
  echo "❌ Errore: parametro sz_hyphothesis deve essere nel formato 'trasmettitore:neuroni' (es. gaba:all)"
  exit 1
fi

# === SIGMA (solo per stimoli non-basali) ===
if [[ "$stim_name" == "basal_activity" ]]; then
    sigmafr=100
else
    sigmafr=10
fi
note=${reserve_sharpness}_${comp_state}_${sz_transmitter}_${sz_neurons}
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
  do_dcnio: ${do_dcnio}
modeling_parameters:
  reserve: "${reserve_str}"
  compensation:
    state: "$comp_state"
    general: "${compensation_str}"
  sz_hyphothesis:
    transmitter: "${sz_transmitter}"
    neurons: "${sz_neurons}"
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
    sz_hyphothesis_str="{'trasmitter': '${sz_transmitter}', 'neurons': '${sz_neurons}'}"
    log_and_print "Launching morphology generation for degeneration $v% (k=$v_normalized)"

    yaml_file="mouse_cerebellar_cortex_${k_suffix}.yaml"
    config_path="./configurations/mouse/"
    if [[ "$do_dcnio" == "1" ]]; then
        hdf5_base_name="mouse_cereb_dcn_io_nest_${k_suffix}"
    else
        hdf5_base_name="mouse_cerebellar_cortex_${k_suffix}"
    fi
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
    
    # === Stimulus generation ===    
    if [[ "$stim_name" == "vbt_stim_protocol_vitro" ]]; then
        python ./scripts_py/make_vbt_stimulus.py --k "$v_normalized" --seed "$this_seed" --stim_name "${stim_name}"
        stim_file="./configurations/mouse/nest/${stim_name}_${k_suffix}.yaml"
    
    elif [[ "$stim_name" == "kase_stim_vitro" ]]; then
        python ./scripts_py/make_kase_stimulus_v2.py --k "$v_normalized" --seed "$this_seed" --stim_name "${stim_name}"
        stim_file="./configurations/mouse/nest/${stim_name}_${k_suffix}.yaml"
    
    elif [[ "$stim_name" == "pawan_stim_vitro" ]]; then
        python ./scripts_py/make_pawan_stimulus_v2.py --k "$v_normalized" --seed "$this_seed" --stim_name "${stim_name}"
        stim_file="./configurations/mouse/nest/${stim_name}_${k_suffix}.yaml"
    fi
    
    if [[ "$stim_name" == "basal_activity" ]]; then
        if [[ "$do_dcnio" == "1" ]]; then
            stim_file="./configurations/mouse/dcn-io/dcn_io_vitro_nest_${k_suffix}.yaml"
        else
            stim_file="./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml"
        fi
    fi
    
    # === MORPHOLOGICAL ALTERATIONS ===
    python ./scripts_py/make_morphologies_v2.py \
        --k "$v_normalized" \
        --reserve "$reserve_k0" "$reserve_sharpness"

    # === NETWORK ALTERATIONS ===
    if [ "$stim_name" = "basal_activity" ]; then
        python ./scripts_py/make_network_v7.py \
            --k "$v_normalized" \
            --do_basal \
            --do_dcnio "$do_dcnio" \
            --reserve "{'k0': $reserve_k0, 'sharpness': $reserve_sharpness}" \
            --compensation "$compensation_str" \
            --sz_hyphothesis "$sz_hyphothesis_str"
    else
        python ./scripts_py/make_network_v7.py \
            --k "$v_normalized" \
            --do_dcnio "$do_dcnio" \
            --reserve "{'k0': $reserve_k0, 'sharpness': $reserve_sharpness}" \
            --compensation "$compensation_str" \
            --sz_hyphothesis "$sz_hyphothesis_str"
    fi


    sim_log="${sim_dir}/simulation.log"
    
    # === Compile ===
    compile_log="${recon_dir}/compile.log"
    mpirun -n "$cores_reconstruction" bsb compile ${stim_file}

    # === Simulation ===
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
    
    # === Move Files ===
    # for nio_file in *.nio; do
    #    [ -f "$nio_file" ] && mv "$nio_file" "$sim_dir/"
    # done
    # echo -e "${YELLOW}📊 Generating dynamics PDF...${NC}"
    # python ./scripts_py/report_dynamics_v4.py "$full_hdf5_path" "$stim_name" "$sim_dir" || echo -e "${RED}❌ PDF generation failed.${NC}"
    # sleep 5
    
    # Move YAML config + morphologies YAML
    [ -f "${config_path}${yaml_file}" ] && mv "${config_path}${yaml_file}" "$config_yaml_dir/"
    [ -f "${config_path}morphologies_${k_suffix}.yaml" ] && mv "${config_path}morphologies_${k_suffix}.yaml" "$config_yaml_dir/"
    
    # Move DCN/IO YAMLs (non-vitro) into config_yaml_dir
    if [[ "$do_dcnio" == "1" ]]; then
        for dcnfile in ./configurations/mouse/nest/dcn-io/dcn_${k_suffix}.yaml ./configurations/mouse/nest/dcn-io/dcn_io_${k_suffix}.yaml; do
            [ -f "$dcnfile" ] && mv "$dcnfile" "$config_yaml_dir/"
        done
    fi
    
    # === Move stimulus YAML into stimulus_dir ===
    if [[ "$stim_name" != "basal_activity" ]]; then
        stim_file_path="./configurations/mouse/nest/${stim_name}_${k_suffix}.yaml"
        [ -f "$stim_file_path" ] && mv "$stim_file_path" "$stimulus_dir/"
    else
        basal_file="./configurations/mouse/nest/basal_vitro_${k_suffix}.yaml"
        [ -f "$basal_file" ] && mv "$basal_file" "$stimulus_dir/"
        if [[ "$do_dcnio" == "1" ]]; then
            # Muovi tutti i file YAML con ${k_suffix} da nest/dcn-io a stimulus_dir
            files_to_move=(./configurations/mouse/dcn-io/*${k_suffix}.yaml)
            if compgen -G "${files_to_move[0]}" > /dev/null; then
                echo "📦 Sposto file YAML con suffisso ${k_suffix} in $stimulus_dir"
                mv ./configurations/mouse/dcn-io/*${k_suffix}.yaml "$stimulus_dir/"
            else
                echo "⚠️ Nessun file trovato con suffisso ${k_suffix} da spostare."
            fi
        fi
    fi

    # === Move stimulus YAML into stimulus_dir ===
    if [[ "$stim_name" != "basal_activity" ]]; then
        python scripts_py/generate_stimulus_timeline.py \
          --stimDir "$yaml_file" \
          --stimName "${stim_file}" \
          --outPath "${stimulus_dir}"
    fi
          
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
#if "$stim_name" == "kase_stim_vitro" ]]; then
  echo "➡️  Running do_completeAnalysis_vbtkase_v2.py ..."
  for sigma in 10 25 50 75 100; do
    echo "▶️  Sigma = $sigma"
    python ./scripts_py/do_completeAnalysis_vbtkase_v2.py --path "$folder_name" --sigma "$sigma" --stimName "$stim_name"
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


#log_and_print "${GREEN}✅ Full pipeline completed.${NC}"
#read -s -n 1 -p "Press any key to exit..."

