#!/bin/bash

main_dir="$1"  # prende il primo argomento da terminale

if [ -z "$main_dir" ]; then
    echo "❌ Devi fornire la directory principale (es: ./simulation_dir)"
    exit 1
fi

echo "📂 Directory principale: $main_dir"

# 1. Genera timeline degli stimoli
python plot_stim_timeline_v3.py "${main_dir}"

# 2. Analizza i dati di firing e salva plot/YAML
python analyze_firing_v2.py "${main_dir}"

# 3. Crea analisi globale e grafici riassuntivi
python overallAnalysis_v2.py "${main_dir}"
