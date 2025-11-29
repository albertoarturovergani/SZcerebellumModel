import os
import glob
import subprocess

# Lista dei percorsi base (cartelle esperimenti)
paths = [
    "./20250613154056_basal_activity_nores_nocomp",
    "./20250613103505_basal_activity_highres_nocomp",
    "./20250614064723_basal_activity_lowres_comp",
    "./20250613123509_basal_activity_lowres_nocomp",
    "./20250613154056_basal_activity_nores_nocomp",
    "./20250613235358_basal_activity_highres_comp"
]

# Path allo script che genera il summary
summary_script = "scripts_py/do_structure_yaml.py"

# Loop su ciascuna cartella
for base_path in paths:
    print(f"🔍 Scansione di: {base_path}")

    if not os.path.isdir(base_path):
        print(f"❌ Cartella non trovata: {base_path}")
        continue

    # Trova tutte le sottocartelle tipo: XXX_deg_YYY/1.structure/
    subfolders = glob.glob(os.path.join(base_path, "*_deg_*/1.structure"))

    for folder in subfolders:
        # Cerca l'unico file .hdf5 nella cartella
        hdf5_files = glob.glob(os.path.join(folder, "*.hdf5"))
        if not hdf5_files:
            print(f"⚠️  Nessun file .hdf5 in: {folder}")
            continue

        hdf5_path = hdf5_files[0]
        print(f"🧠 Analisi struttura: {hdf5_path}")

        # Esegui lo script Python
        result = subprocess.run(["python", summary_script, "--path", hdf5_path])

        if result.returncode == 0:
            print("✅ YAML salvato con successo.")
        else:
            print(f"❌ Errore durante analisi di: {hdf5_path}")
