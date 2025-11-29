import os
import sys
from cerebellum.analysis.structure_analysis import StructureReport
from bsb import from_storage

# Verifica se il percorso del file HDF5 è fornito come argomento
if len(sys.argv) < 2:
    print("Usage: python report.py <path-to-hdf5-file>")
    sys.exit(1)

# Ottieni il percorso completo del file HDF5
hdf5_path = sys.argv[1]

# Estrai il nome del file HDF5 senza estensione
hdf5_filename = os.path.splitext(os.path.basename(hdf5_path))[0]

# Crea il percorso del file PDF aggiungendo il prefisso 'report_' e l'estensione '.pdf'
pdf_report_path = os.path.join(os.path.dirname(hdf5_path), f"report_{hdf5_filename}.pdf")

# Carica il modello e crea il report
scaffold = from_storage(hdf5_path)
report = StructureReport(scaffold)

# Stampa il report
report.print_report(pdf_report_path)

