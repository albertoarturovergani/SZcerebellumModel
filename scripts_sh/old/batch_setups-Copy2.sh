#!/bin/bash

# === CONFIGURAZIONI ===
stim_name="pawan_stim_vitro"
start=40
stop=40
n_steps=1

# === COMBINAZIONI DI PARAMETRI ===
declare -a compensations=(0)
declare -a gabas=("glutamate:all")
declare -a dcnio_flags=("1")
declare -a reserve_k0_list=("0.5") #1 
declare -a reserve_sharpness_list=("1") #30

# === TIMER INIZIALE ===
start_time=$(date +%s)

# === CICLO COMBINATO ===
count=1
total=$(( ${#compensations[@]} * ${#gabas[@]} * ${#dcnio_flags[@]} * ${#reserve_k0_list[@]} * ${#reserve_sharpness_list[@]} ))

for gaba in "${gabas[@]}"; do
    for comp in "${compensations[@]}"; do
        for do_dcnio in "${dcnio_flags[@]}"; do
            for reserve_k0 in "${reserve_k0_list[@]}"; do
                for reserve_sharpness in "${reserve_sharpness_list[@]}"; do

                    # Calcola tempo trascorso
                    current_time=$(date +%s)
                    elapsed=$((current_time - start_time))
                    elapsed_min=$((elapsed / 60))
                    elapsed_sec=$((elapsed % 60))

                    # Stampa variante e tempo trascorso
                    echo -e "\n▶️ Variant $count/$total: $gaba + compensation=$comp + dcnio=$do_dcnio + k0=$reserve_k0 + sharp=$reserve_sharpness"
                    echo "⏱️ Tempo trascorso: ${elapsed_min}m ${elapsed_sec}s"

                    # Lancia lo script principale con TUTTI i parametri
                    bash scripts_sh/sim_setup_v4.sh "$stim_name" "$start" "$stop" "$n_steps" "$comp" "$gaba" "$do_dcnio" "$reserve_k0" "$reserve_sharpness"
                    # usa per basal bash scripts_sh/sim_setup_v3.sh "$stim_name" "$start" "$stop" "$n_steps" "$comp" "$gaba" "$do_dcnio" "$reserve_k0" "$reserve_sharpness"

                    ((count++))
                done
            done
        done
    done
done

echo -e "\n✅ Tutte le simulazioni sono state lanciate."
