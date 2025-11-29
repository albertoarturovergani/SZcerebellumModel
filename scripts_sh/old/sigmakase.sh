#!/usr/bin/env bash


# Elenco delle cartelle da analizzare
folders=(
  "20251113161643_kase_stim_vitro_2_0_glutamate_all"
  "20251113214858_kase_stim_vitro_15_0_glutamate_all"
)

# Sigma che vuoi usare
sigmas=($(seq 1 1 100))

# Stim name (di’ esplicitamente cosa stai usando)
stim_name="kase_stim_vitro"

if [[ "$stim_name" == "vbt_stim_protocol_vitro" || "$stim_name" == "kase_stim_vitro" ]]; then
  echo "➡️  Running overallAnalysis_v3.py ..."
  for folder_name in "${folders[@]}"; do
    for sigma in "${sigmas[@]}"; do
      echo "▶️  Folder = $folder_name | Sigma = $sigma"
      python ./scripts_py/do_completeAnalysis_vbtkase_v2.py \
        --path "$folder_name" \
        --sigma "$sigma" \
        --stimName "$stim_name"
    done
  done
fi

