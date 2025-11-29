#!/bin/bash

# Requisiti: enscript, ps2pdf, pdfunite
for cmd in pdfunite enscript ps2pdf; do
    command -v $cmd >/dev/null 2>&1 || { echo "❌ $cmd non trovato. Installa con: sudo apt install poppler-utils enscript ghostscript"; exit 1; }
done

results_dir=${1:-"."}
timestamp=$(date +%Y%m%d_%H%M%S)
output_pdf="FullReport_${timestamp}.pdf"
temp_dir="temp_report_merge"
mkdir -p "$temp_dir"

merged_list=()

results_folders=$(find "$results_dir" -maxdepth 2 -type d -name "results_*" | sort -V)

for folder in $results_folders; do
    rel=$(basename "$folder" | grep -oP '\d+\.\d+')
    recon_dir="$folder/reconstruction"

    echo "🔹 Processing rel = ${rel}"

    # 1. YAML → PDF
    yaml_file=$(find "$recon_dir" -maxdepth 1 -name "*.yaml" 2>/dev/null)
    if [ -f "$yaml_file" ]; then
        yaml_pdf="${temp_dir}/01_yaml_${rel}.pdf"
        enscript -B -q -p - "$yaml_file" | ps2pdf - "$yaml_pdf"
        merged_list+=("$yaml_pdf")
    fi

    # 2. Reconstruction PDF
    recon_pdf=$(find "$recon_dir" -maxdepth 1 -name "*.pdf" 2>/dev/null)
    if [ -f "$recon_pdf" ]; then
        recon_pdf_copy="${temp_dir}/02_reconstruction_${rel}.pdf"
        cp "$recon_pdf" "$recon_pdf_copy"
        merged_list+=("$recon_pdf_copy")
    fi

    # 3. Simulation PDFs
    sim_dirs=$(find "$folder" -mindepth 1 -maxdepth 1 -type d ! -name "reconstruction")
    for sim in $sim_dirs; do
        sim_name=$(basename "$sim")
        sim_pdf=$(find "$sim" -maxdepth 1 -name "*.pdf" 2>/dev/null)
        if [ -f "$sim_pdf" ]; then
            sim_pdf_copy="${temp_dir}/03_sim_${rel}_${sim_name}.pdf"
            cp "$sim_pdf" "$sim_pdf_copy"
            merged_list+=("$sim_pdf_copy")
        fi
    done
done

# Unione finale
if [ ${#merged_list[@]} -eq 0 ]; then
    echo "❌ Nessun file trovato da unire."
    exit 1
fi

echo "📄 Generazione PDF finale: $output_pdf"
pdfunite "${merged_list[@]}" "$output_pdf"

# Pulizia
rm -r "$temp_dir"
echo "✅ Report completo creato: $output_pdf"
