#!/bin/bash

name="vbt_test"
span_mode="linear" # oppure "log"
start=0
stop=50
n_steps=10

cores_reconstruction=8
cores_simulations=8

datetime=$(date +"%Y%m%d%H%M%S")
folder_name="${datetime}_${name}"
mkdir -p "$folder_name"

log_file="${folder_name}/log.txt"
touch "$log_file"

log_and_print() {
    echo "$1"
    echo "$1" >> "$log_file"
}

log_and_print "[INFO] Folder '$folder_name' created."
log_and_print "[INFO] Log file '$log_file' created."

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

log_and_print "[INFO] Generated span values:"
index=0
for v in "${values[@]}"; do
    v_rounded=$(LC_NUMERIC=C printf "%.2f" "$v")
    echo " [$index] $v_rounded" | tee -a "$log_file"
    ((index++))
done

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

index=0
for seed in "${seeds_list[@]}"; do
    v_rounded=$(LC_NUMERIC=C printf "%.2f" "${values[$index]}")
    int_part=$(LC_NUMERIC=C printf "%02d" $(echo "$v_rounded" | awk -F. '{print $1}'))
    decimal_part=$(LC_NUMERIC=C printf "%02d" $(echo "$v_rounded" | awk -F. '{print $2}'))
    val_formatted="${int_part}_${decimal_part}"
    index_formatted=$(printf "%03d" "$index")
    subfolder="${folder_name}/${index_formatted}_deg_${val_formatted}"
    mkdir -p "$subfolder"/{1.config_files,2.stimulus,3.structure,4.function}
    log_and_print "[INFO] Created subfolder '${index_formatted}_deg_${val_formatted}' with internal structure."

    log_deg_file="${subfolder}/log_deg_${val_formatted}.txt"
    echo "Seed: ${seed}" > "$log_deg_file"
    log_and_print "[INFO] Created log file 'log_deg_${val_formatted}.txt' inside '${index_formatted}_deg_${val_formatted}' with seed."
    ((index++))
done

# Copy and rename morphologies inside each structure folder
if [[ -d "morphologies" && "$(ls -A morphologies)" ]]; then
    for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
        if [[ -d "$degfolder/3.structure" ]]; then
            suffix=$(basename "$degfolder" | cut -d'_' -f3,4) # estrae il XX_XX
            for file in morphologies/*; do
                if [[ -f "$file" ]]; then
                    filename=$(basename "$file")
                    extension="${filename##*.}"
                    filename_base="${filename%.*}"
                    cp "$file" "$degfolder/3.structure/${filename_base}_${suffix}.${extension}"
                fi
            done
            log_and_print "[INFO] Copied and renamed morphologies into '$degfolder/3.structure' with suffix '_${suffix}'."
        fi
    done
else
    log_and_print "[WARNING] 'morphologies/' folder not found or is empty. No structures copied."
fi

# Copy only stimulus YAMLs containing "vitro" inside each stimulus folder and rename with degeneration
if [[ -d "configurations/mouse/nest" && -f "configurations/mouse/dcn-io/dcn_io_vitro_nest.yaml" ]]; then
    for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
        if [[ -d "$degfolder/2.stimulus" ]]; then
            degeneration_suffix=$(basename "$degfolder" | awk -F'_' '{print $(NF-1)"_"$NF}')

            # Copy only *vitro*.yaml from nest/
            for yamlfile in configurations/mouse/nest/*vitro*.yaml; do
                if [[ -f "$yamlfile" ]]; then
                    filename=$(basename "$yamlfile")
                    base="${filename%.*}"
                    ext="${filename##*.}"
                    new_filename="${base}_${degeneration_suffix}.${ext}"
                    cp "$yamlfile" "$degfolder/2.stimulus/$new_filename"
                    log_and_print "[INFO] Copied and renamed '$filename' -> '$new_filename' into '$degfolder/2.stimulus'."
                fi
            done

            # Copy the specific dcn-io yaml
            yamlfile="configurations/mouse/dcn-io/dcn_io_vitro_nest.yaml"
            if [[ -f "$yamlfile" ]]; then
                filename=$(basename "$yamlfile")
                base="${filename%.*}"
                ext="${filename##*.}"
                new_filename="${base}_${degeneration_suffix}.${ext}"
                cp "$yamlfile" "$degfolder/2.stimulus/$new_filename"
                log_and_print "[INFO] Copied and renamed 'dcn_io_vitro_nest.yaml' -> '$new_filename' into '$degfolder/2.stimulus'."
            fi

        fi
    done
else
    log_and_print "[WARNING] No stimulus YAML files with 'vitro' found. No stimulus files copied."
fi

# Copy configuration YAMLs into each config_files folder and rename with degeneration
if [[ -f "configurations/mouse/mouse_cerebellar_cortex.yaml" && -f "configurations/mouse/dcn-io/dcn_io.yaml" ]]; then
    for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
        if [[ -d "$degfolder/1.config_files" ]]; then
            degeneration_suffix=$(basename "$degfolder" | awk -F'_' '{print $(NF-1)"_"$NF}')

            # Copy mouse_cerebellar_cortex.yaml
            yamlfile="configurations/mouse/mouse_cerebellar_cortex.yaml"
            filename=$(basename "$yamlfile")
            base="${filename%.*}"
            ext="${filename##*.}"
            new_filename="${base}_${degeneration_suffix}.${ext}"
            cp "$yamlfile" "$degfolder/1.config_files/$new_filename"
            log_and_print "[INFO] Copied and renamed '$filename' -> '$new_filename' into '$degfolder/1.config_files'."

            # Copy dcn_io.yaml
            yamlfile="configurations/mouse/dcn-io/dcn_io.yaml"
            filename=$(basename "$yamlfile")
            base="${filename%.*}"
            ext="${filename##*.}"
            new_filename="${base}_${degeneration_suffix}.${ext}"
            cp "$yamlfile" "$degfolder/1.config_files/$new_filename"
            log_and_print "[INFO] Copied and renamed 'dcn_io.yaml' -> '$new_filename' into '$degfolder/1.config_files'."
        fi
    done
else
    log_and_print "[WARNING] Configuration files not found. No config files copied."
fi


collective_folder="${folder_name}/B_collective_results"
mkdir -p "$collective_folder"/{1.structures,2.functions}
log_and_print "[INFO] Created 'B_collective_results' with '1.structures' and '2.functions'."

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
    echo "# Example of using seeds:"
    echo "for seed in \"\${SEEDS[@]}\"; do"
    echo "  echo \"Running simulation with seed \$seed\""
    echo "  # here your simulation code using \$seed"
    echo "done"
    echo ""
    cat "scripts_sh/sims_run.sh"
} > "$replica_sims_script"

chmod +x "$replica_sims_script"
log_and_print "[INFO] Created '${datetime}_run_sims.sh' with embedded seeds and script copy."

if [[ -f "scripts_sh/juelich_run.sh" ]]; then
    {
        echo "# Copied on ${datetime_now} by user '${user_now}'"
        sed "s|source \$XXX/run_sims.sh|source \$XXX/${datetime}_run_sims.sh|" "scripts_sh/juelich_run.sh"
    } > "${replica_folder}/juelich_run.sh"
    chmod +x "${replica_folder}/juelich_run.sh"
    log_and_print "[INFO] Copied and updated 'juelich_run.sh' with corrected source."
else
    log_and_print "[ERROR] 'scripts_sh/juelich_run.sh' not found!"
fi

log_replica_file="${replica_folder}/log_replica.txt"

{
    echo "To submit the simulation to the cluster:"
    echo ""
    echo "1. Move into your cerebellum working directory:"
    echo "   cd /path/to/your/cerebellum"
    echo ""
    echo "2. Activate your BSB environment:"
    echo "   source BSB/bin/activate"
    echo ""
    echo "3. Move into the replica folder:"
    echo "   cd ${replica_folder}"
    echo ""
    echo "4. Submit the batch job:"
    echo "   sbatch juelich_run.sh"
} > "$log_replica_file"

log_and_print "[INFO] Created 'log_replica.txt' inside 'A_replica'."

# CREAZIONE info.yaml COMPLETO
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
    echo ""
    echo "deg_sweep:"
    for v in "${values[@]}"; do
        v_rounded=$(LC_NUMERIC=C printf "%.2f" "$v")
        echo "  - ${v_rounded}"
    done
    echo ""
    echo "seeds:"
    for seed in "${seeds_list[@]}"; do
        echo "  - ${seed}"
    done
    echo ""
    echo "simulation_timings: []"
    echo "time_reconstruction: []"
    echo "time_simulation: []"
    echo "time_postprocess: []"
    echo "time_total: []"
} > "$global_config_file"

# APPEND FOLDER SIZES
{
    echo ""
    echo "folder_sizes:"
    total_size=$(du -sb "$folder_name" | awk '{print $1}')
    echo "  total_size: ${total_size}"

    for sub in A_replica B_collective_results C_seeds; do
        if [[ -d "${folder_name}/${sub}" ]]; then
            size=$(du -sb "${folder_name}/${sub}" | awk '{print $1}')
            echo "  ${sub}: ${size}"
        fi
    done

    for degfolder in "${folder_name}"/[0-9][0-9][0-9]_deg_*; do
        if [[ -d "$degfolder" ]]; then
            foldername_only=$(basename "$degfolder")
            size=$(du -sb "$degfolder" | awk '{print $1}')
            echo "  ${foldername_only}: ${size}"
        fi
    done
} >> "$global_config_file"

log_and_print "[INFO] Completed 'info.yaml' with folder sizes."
