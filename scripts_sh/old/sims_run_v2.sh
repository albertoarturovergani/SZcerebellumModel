#!/bin/bash

# === SETTINGS ===
name="vbt_test"
span_mode="linear" # oppure "log"
start=0
stop=50
n_steps=1

cores_reconstruction=8
cores_simulations=8

# === CREATE MAIN FOLDER ===
datetime=$(date +"%Y%m%d%H%M%S")
folder_name="${datetime}_${name}"
mkdir -p "$folder_name"

log_file="${folder_name}/log.txt"
touch "$log_file"

log_and_print() {
    echo "$1"
    echo "$1" >> "$log_file"
}

log_and_print "[INFO] Created folder '$folder_name' and log."

# === GENERATE SPAN VALUES ===
values=()

if [[ "$span_mode" == "linear" ]]; then
    step=$(echo "($stop - $start) / ($n_steps - 1)" | bc -l)
    val=$start
    for ((i=0; i<n_steps; i++)); do
        values+=("$val")
        val=$(echo "$val + $step" | bc -l)
    done
elif [[ "$span_mode" == "log" ]]; then
    log_start=$(echo "l($start+1)" | bc -l)
    log_stop=$(echo "l($stop+1)" | bc -l)
    for ((i=0; i<n_steps; i++)); do
        fraction=$(echo "$i / ($n_steps - 1)" | bc -l)
        log_val=$(echo "$log_start + ($log_stop - $log_start) * $fraction" | bc -l)
        val=$(echo "e($log_val) - 1" | bc -l)
        values+=("$val")
    done
else
    log_and_print "[ERROR] Invalid span mode selected: $span_mode"
    exit 1
fi

log_and_print "[INFO] Generated span values."

# === INFO YAML ===
global_config_file="${folder_name}/info.yaml"
{
    echo "experiment:"
    echo "  name: ${name}"
    echo "  date: ${datetime}"
    echo "  user: ${USER}"
    echo "  span_mode: ${span_mode}"
    echo "  start: ${start}"
    echo "  stop: ${stop}"
    echo "  n_steps: ${n_steps}"
    echo "  cores_reconstruction: ${cores_reconstruction}"
    echo "  cores_simulations: ${cores_simulations}"
    echo "deg_sweep:"
    for v in "${values[@]}"; do
        v_rounded=$(LC_NUMERIC=C printf "%.2f" "$v")
        echo "  - ${v_rounded}"
    done
    echo "simulation_timings: []"
    echo "time_reconstruction: []"
    echo "time_simulation: []"
    echo "time_postprocess: []"
    echo "time_total: []"
} > "$global_config_file"

log_and_print "[INFO] Created info.yaml."

# === CREATE SEEDS ===
year=${datetime:0:4}
month=${datetime:4:2}
day=${datetime:6:2}
hour=${datetime:8:2}
minute=${datetime:10:2}
second=${datetime:12:2}

seeds_folder="${folder_name}/C_seeds"
mkdir -p "$seeds_folder"
seeds_file="${seeds_folder}/seeds.txt"
touch "$seeds_file"

seeds_list=()

index=0
for v in "${values[@]}"; do
    index_formatted=$(printf "%03d" "$index")
    rand1=$(printf "%02d" $((RANDOM % 100)))
    rand2=$(printf "%02d" $((RANDOM % 100)))
    rand3=$(printf "%02d" $((RANDOM % 100)))
    rand4=$(printf "%02d" $((RANDOM % 100)))
    rand5=$(printf "%02d" $((RANDOM % 100)))
    seed_full="${index_formatted}${year}${month}${day}${hour}${minute}${second}${rand1}${rand2}${rand3}${rand4}${rand5}"
    seeds_list+=("$seed_full")
    echo "$seed_full" >> "$seeds_file"
    ((index++))
done

log_and_print "[INFO] Created seeds."

# === CREATE SUBFOLDERS FOR EACH DEGENERATION ===
index=0
for seed in "${seeds_list[@]}"; do
    v_rounded=$(LC_NUMERIC=C printf "%.2f" "${values[$index]}")
    int_part=$(LC_NUMERIC=C printf "%02d" $(echo "$v_rounded" | awk -F. '{print $1}'))
    decimal_part=$(LC_NUMERIC=C printf "%02d" $(echo "$v_rounded" | awk -F. '{print $2}'))
    val_formatted="${int_part}_${decimal_part}"
    index_formatted=$(printf "%03d" "$index")
    subfolder="${folder_name}/${index_formatted}_deg_${val_formatted}"
    mkdir -p "$subfolder"/{1.config_files,2.stimulus,3.structure,4.function}
    echo "Seed: ${seed}" > "$subfolder/log_deg_${val_formatted}.txt"
    log_and_print "[INFO] Created '$subfolder' and log file."
    ((index++))
done

# === CREATE B_collective_results ===
collective_folder="${folder_name}/B_collective_results"
mkdir -p "$collective_folder"/{1.structures,2.functions}
log_and_print "[INFO] Created 'B_collective_results'."

# === CREATE A_replica and copy scripts ===
replica_folder="${folder_name}/A_replica"
mkdir -p "$replica_folder"

datetime_now=$(date +"%Y-%m-%d %H:%M:%S")
user_now="$USER"

replica_sims_script="${replica_folder}/${datetime}_run_sims.sh"

{
    echo "# Copied on ${datetime_now} by user '${user_now}'"
    echo "# Seed base: ${datetime}"
    echo "# Seeds list:"
    for seed_line in "${seeds_list[@]}"; do
        echo "# $seed_line"
    done
    echo ""
    echo "SEEDS=("
    for seed in "${seeds_list[@]}"; do
        echo "  \"$seed\""
    done
    echo ")"
    echo ""
    echo "# Example loop:"
    echo "for seed in \"\${SEEDS[@]}\"; do"
    echo "  echo \"Running simulation with seed \$seed\""
    echo "done"
    echo ""
    cat "scripts_sh/sims_run.sh"
} > "$replica_sims_script"

chmod +x "$replica_sims_script"

# Copy and patch juelich_run.sh
if [[ -f "scripts_sh/juelich_run.sh" ]]; then
    {
        echo "# Copied on ${datetime_now} by user '${user_now}'"
        sed "s|source \$XXX/run_sims.sh|source \$XXX/${datetime}_run_sims.sh|" "scripts_sh/juelich_run.sh"
    } > "${replica_folder}/juelich_run.sh"
    chmod +x "${replica_folder}/juelich_run.sh"
    log_and_print "[INFO] Copied and updated 'juelich_run.sh' into 'A_replica'."
else
    log_and_print "[WARNING] 'scripts_sh/juelich_run.sh' not found!"
fi

# === CREATE log_replica.txt ===
{
    echo "To submit the simulation:"
    echo ""
    echo "cd cerebellum"
    echo "source BSB/bin/activate"
    echo "cd ${folder_name}/A_replica"
    echo "sbatch juelich_run.sh"
} > "${replica_folder}/log_replica.txt"

log_and_print "[INFO] Created 'log_replica.txt'."

for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
    if [[ -d "$degfolder" ]]; then
        degeneration_suffix=$(basename "$degfolder" | awk -F'_' '{print $(NF-1)"_"$NF}')

        # COPY MORPHOLOGIES
        if [[ -d "morphologies" && "$(ls -A morphologies)" ]]; then
            suffix=$(basename "$degfolder" | cut -d'_' -f3,4)
            for file in morphologies/*; do
                if [[ -f "$file" ]]; then
                    filename=$(basename "$file")
                    extension="${filename##*.}"
                    filename_base="${filename%.*}"
                    cp "$file" "$degfolder/3.structure/${filename_base}_${suffix}.${extension}"
                fi
            done
            log_and_print "[INFO] Copied and renamed morphologies into '$degfolder/3.structure' with suffix '_${suffix}'."
        else
            log_and_print "[WARNING] 'morphologies/' folder not found or is empty. No structures copied."
        fi

        # === COPY stimulus files containing 'vitro' from /nest and specific files ===
        if [[ -d "configurations/mouse/nest/" ]]; then
            for stim in configurations/mouse/nest/*vitro*.yaml; do
                if [[ -f "$stim" ]]; then
                    base_stim=$(basename "$stim" .yaml)
                    cp "$stim" "$degfolder/2.stimulus/${base_stim}_${degeneration_suffix}.yaml"
                fi
            done
            log_and_print "[INFO] Copied 'vitro' stimulus files from 'nest' into '$degfolder/2.stimulus'."
        else
            log_and_print "[WARNING] 'configurations/mouse/nest/' not found. No stimulus copied."
        fi
        
        # Copia anche specificamente dcn_io_vitro_nest.yaml se esiste
        if [[ -f "configurations/mouse/dcn-io/dcn_io_vitro_nest.yaml" ]]; then
            cp "configurations/mouse/dcn-io/dcn_io_vitro_nest.yaml" "$degfolder/2.stimulus/dcn_io_vitro_nest_${degeneration_suffix}.yaml"
            log_and_print "[INFO] Copied 'dcn_io_vitro_nest.yaml' into '$degfolder/2.stimulus'."
        fi

        # COPY AND PATCH CONFIG FILES
        if [[ -f "configurations/mouse/mouse_cerebellar_cortex.yaml" && -f "configurations/mouse/dcn-io/dcn_io.yaml" && -f "configurations/mouse/dcn-io/dcn.yaml" ]]; then

            # mouse_cerebellar_cortex
            yamlfile="configurations/mouse/mouse_cerebellar_cortex.yaml"
            base_yaml=$(basename "$yamlfile" .yaml)
            cp "$yamlfile" "$degfolder/1.config_files/${base_yaml}_${degeneration_suffix}.yaml"
            
            # Patch inside copied mouse_cerebellar_cortex yaml
            sed -i "s|mouse_cerebellum.hdf5|mouse_cerebellum_${degeneration_suffix}.hdf5|g" "$degfolder/1.config_files/${base_yaml}_${degeneration_suffix}.yaml"
            sed -i "s|bsb_report_structure.pdf|bsb_report_structure_${degeneration_suffix}.pdf|g" "$degfolder/1.config_files/${base_yaml}_${degeneration_suffix}.yaml"
            sed -i "s|\.\./\.\./morphologies/morphologies.yaml#/|../3.structure/morphologies_${degeneration_suffix}.yaml#/|g" "$degfolder/1.config_files/${base_yaml}_${degeneration_suffix}.yaml"


            # dcn_io
            yamlfile="configurations/mouse/dcn-io/dcn_io.yaml"
            base_yaml=$(basename "$yamlfile" .yaml)
            cp "$yamlfile" "$degfolder/1.config_files/${base_yaml}_${degeneration_suffix}.yaml"
            sed -i "s|mouse_cereb_dcn_io.hdf5|mouse_cereb_dcn_io_${degeneration_suffix}.hdf5|g" "$degfolder/1.config_files/${base_yaml}_${degeneration_suffix}.yaml"
            sed -i "s|./dcn.yaml#/|./dcn_${degeneration_suffix}.yaml#/|g" "$degfolder/1.config_files/${base_yaml}_${degeneration_suffix}.yaml"

            # dcn
            yamlfile="configurations/mouse/dcn-io/dcn.yaml"
            base_yaml=$(basename "$yamlfile" .yaml)
            cp "$yamlfile" "$degfolder/1.config_files/${base_yaml}_${degeneration_suffix}.yaml"
            sed -i "s|mouse_cereb_dcn.hdf5|mouse_cereb_dcn_${degeneration_suffix}.hdf5|g" "$degfolder/1.config_files/${base_yaml}_${degeneration_suffix}.yaml"
            sed -i "s|\.\./mouse_cerebellar_cortex.yaml#|./mouse_cerebellar_cortex_${degeneration_suffix}.yaml#|g" "$degfolder/1.config_files/${base_yaml}_${degeneration_suffix}.yaml"

            log_and_print "[INFO] Copied and patched config files into '$degfolder/1.config_files'."
        else
            log_and_print "[WARNING] Missing configuration files for patching!"
        fi
    fi
done

# === PATCH basal_vitro_XX_XX.yaml ===
for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
    degeneration_suffix=$(basename "$degfolder" | awk -F'_' '{print $(NF-1)"_"$NF}')

    stimulus_folder="$degfolder/2.stimulus"

    for basal_file in "$stimulus_folder"/basal_vitro_*"${degeneration_suffix}".yaml; do
        if [[ -f "$basal_file" ]]; then
            sed -i "s|mouse_cereb_nest.hdf5|mouse_cereb_nest_${degeneration_suffix}.hdf5|g" "$basal_file"
            sed -i "s|\.\./mouse_cerebellar_cortex.yaml#|../1.config_files/mouse_cerebellar_cortex_${degeneration_suffix}.yaml#|g" "$basal_file"
            log_and_print "[INFO] Patched '$(basename "$basal_file")' with correct storage and import paths."
        fi
    done
done

# === PATCH stimulus_mossy_XX_XX.yaml ===
for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
    degeneration_suffix=$(basename "$degfolder" | awk -F'_' '{print $(NF-1)"_"$NF}')

    stimulus_folder="$degfolder/2.stimulus"

    for mossy_file in "$stimulus_folder"/stimulus_mossy_*"${degeneration_suffix}".yaml; do
        if [[ -f "$mossy_file" ]]; then
            sed -i "s|\.\/basal_vitro.yaml#/simulations/basal_activity|./basal_vitro_${degeneration_suffix}.yaml#/simulations/basal_activity|g" "$mossy_file"
            sed -i "s|\.\/basal_vitro.yaml#/|./basal_vitro_${degeneration_suffix}.yaml#/|g" "$mossy_file"
            log_and_print "[INFO] Patched '$(basename "$mossy_file")' with correct import references."
        fi
    done
done

# === PATCH dcn_io_vitro_nest_XX_XX.yaml ===
for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
    degeneration_suffix=$(basename "$degfolder" | awk -F'_' '{print $(NF-1)"_"$NF}')
    stimulus_folder="$degfolder/2.stimulus"

    for dcnio_file in "$stimulus_folder"/dcn_io_vitro_nest_*"${degeneration_suffix}".yaml; do
        if [[ -f "$dcnio_file" ]]; then
            sed -i "s|\.\/dcn_vitro_nest.yaml#/simulations/basal_activity/|./dcn_vitro_nest_${degeneration_suffix}.yaml#/simulations/basal_activity/|g" "$dcnio_file"
            sed -i "s|\"#/simulations/basal_activity\"|\"./dcn_io_vitro_nest_${degeneration_suffix}.yaml#/simulations/basal_activity\"|g" "$dcnio_file"
            sed -i "s|\.\/dcn_io.yaml#/|./dcn_io_${degeneration_suffix}.yaml#/|g" "$dcnio_file"
            sed -i "s|mouse_cereb_dcn_io_nest.hdf5|mouse_cereb_dcn_io_nest_${degeneration_suffix}.hdf5|g" "$dcnio_file"
            log_and_print "[INFO] Patched '$(basename "$dcnio_file")' with correct imports and storage."
        fi
    done
done

# === COPY AND GENERATE MORPHOLOGIES YAML ===
if [[ -d "morphologies" && "$(ls -A morphologies)" ]]; then
    for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
        if [[ -d "$degfolder/3.structure" ]]; then
            suffix=$(basename "$degfolder" | cut -d'_' -f3,4)
            morph_yaml="${degfolder}/3.structure/morphologies_${suffix}.yaml"

            echo "morphologies:" > "$morph_yaml"

            for file in morphologies/*.swc; do
                if [[ -f "$file" ]]; then
                    filename=$(basename "$file" .swc)
                    new_filename="${filename}_${suffix}.swc"
                    cp "$file" "$degfolder/3.structure/$new_filename"

                    echo "  - file: ./${new_filename}" >> "$morph_yaml"
                    echo "    pipeline:" >> "$morph_yaml"
                    echo "      - func: swap_axes" >> "$morph_yaml"
                    echo "        axis1: 1" >> "$morph_yaml"
                    echo "        axis2: 2" >> "$morph_yaml"
                    echo "    parser:" >> "$morph_yaml"
                    case "$filename" in
                        "GranuleCell")
                            echo "      tags:" >> "$morph_yaml"
                            echo "        16: [axon, axon_hillock]" >> "$morph_yaml"
                            echo "        17: [axon, axon_initial_segment]" >> "$morph_yaml"
                            echo "        18: [axon, ascending_axon]" >> "$morph_yaml"
                            echo "        19: [axon, parallel_fiber]" >> "$morph_yaml"
                            ;;
                        "GolgiCell")
                            echo "      tags:" >> "$morph_yaml"
                            echo "        16: [dendrites, basal_dendrites]" >> "$morph_yaml"
                            echo "        17: [dendrites, apical_dendrites]" >> "$morph_yaml"
                            echo "        18: [axon, axon_initial_segment]" >> "$morph_yaml"
                            ;;
                        "PurkinjeCell")
                            echo "      tags:" >> "$morph_yaml"
                            echo "        16: [axon, AIS]" >> "$morph_yaml"
                            echo "        17: [axon, AIS_K]" >> "$morph_yaml"
                            echo "        18: [axon, axonmyelin]" >> "$morph_yaml"
                            echo "        19: [axon, nodes]" >> "$morph_yaml"
                            echo "        20: [dendrites, basal]" >> "$morph_yaml"
                            echo "        21: [dendrites, pf_targets, sc_targets]" >> "$morph_yaml"
                            echo "        22: [dendrites, aa_targets, sc_targets]" >> "$morph_yaml"
                            ;;
                        "StellateCell")
                            echo "      tags:" >> "$morph_yaml"
                            echo "        16: [dendrites, proximal_dendrites]" >> "$morph_yaml"
                            echo "        17: [dendrites, distal_dendrites]" >> "$morph_yaml"
                            echo "        18: [axon, axon_initial_segment]" >> "$morph_yaml"
                            ;;
                        "BasketCell")
                            echo "      tags:" >> "$morph_yaml"
                            echo "        16: [axon, axon_initial_segment]" >> "$morph_yaml"
                            ;;
                        *)
                            echo "      tags: {}" >> "$morph_yaml"
                            ;;
                    esac
                fi
            done

            log_and_print "[INFO] Copied and generated 'morphologies_${suffix}.yaml' in '$degfolder/3.structure'."
        fi
    done
else
    log_and_print "[WARNING] 'morphologies/' folder not found or is empty. No structures copied."
fi

# === FIX IMPORT PATH in mouse_cerebellar_cortex_XX_XX.yaml ===
for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
    degeneration_suffix=$(basename "$degfolder" | awk -F'_' '{print $(NF-1)"_"$NF}')
    cortex_file="$degfolder/1.config_files/mouse_cerebellar_cortex_${degeneration_suffix}.yaml"
    if [[ -f "$cortex_file" ]]; then
        sed -i "s|\.\./\.\./morphologies/morphologies.yaml#|../3.structure/morphologies_${degeneration_suffix}.yaml#|g" "$cortex_file"
        log_and_print "[INFO] Fixed \$import in 'mouse_cerebellar_cortex_${degeneration_suffix}.yaml'."
    fi
done


for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
    if [[ -d "$degfolder" ]]; then
        degeneration_suffix=$(basename "$degfolder" | awk -F'_' '{print $(NF-1)"_"$NF}')
        
        config_file="1.config_files/mouse_cerebellar_cortex_${degeneration_suffix}.yaml"
        hdf5_file="4.function/mouse_cereb_nest_${degeneration_suffix}.hdf5"
        stim_file="2.stimulus/basal_vitro_${degeneration_suffix}.yaml"
        sim_name="basal_activity"
        
        if [[ -f "$degfolder/$config_file" && -f "$degfolder/$stim_file" ]]; then
            log_and_print "[INFO] Compiling network for ${degfolder}..."
            cd "$degfolder"
            mpirun -n ${cores_reconstruction} bsb compile "$config_file" || {
                log_and_print "[ERROR] Compile failed on ${degfolder}, skipping."
                cd -
                continue
            }

            log_and_print "[INFO] Running reconfigure on ${degfolder}..."
            mpirun -n ${cores_reconstruction} bsb reconfigure "$hdf5_file" "$stim_file" || {
                log_and_print "[ERROR] Reconfigure failed on ${degfolder}, skipping simulate."
                cd -
                continue
            }
            
            log_and_print "[INFO] Running simulate on ${degfolder}..."
            mpirun -n ${cores_simulations} bsb simulate "$hdf5_file" "$sim_name" || {
                log_and_print "[ERROR] Simulation failed on ${degfolder}"
                cd -
                continue
            }
            cd -
        else
            log_and_print "[WARNING] Missing config or stimulus file in ${degfolder}, skipping."
        fi
    fi
done




# === UPDATE info.yaml WITH FINAL SIZES ===
total_bytes=$(du -sb "$folder_name" | awk '{print $1}')
echo "" >> "$global_config_file"
echo "folder_sizes:" >> "$global_config_file"
echo "  total_bytes: $total_bytes" >> "$global_config_file"
for dir in "$folder_name"/*/; do
    size=$(du -sb "$dir" | awk '{print $1}')
    name=$(basename "$dir")
    echo "  ${name}: $size" >> "$global_config_file"
done

log_and_print "[INFO] Updated 'info.yaml' with folder sizes."

# === DONE ===
log_and_print "[INFO] Finished simulation structure preparation."



