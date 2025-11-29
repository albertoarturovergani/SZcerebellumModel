#!/usr/bin/env python3

import sys
import os
from cerebellum.analysis.spiking_results import BasicSimulationReport
from bsb import from_storage

# Ensure correct number of arguments
if len(sys.argv) != 4:
    print("Usage: python myFile.py <hdf5_file> <simulation_name> <output_folder>")
    print("Example: python myFile.py mouse_cerebellar_cortex.hdf5 basal_activity results_folder")
    sys.exit(1)

# Read input arguments
hdf5_filename = sys.argv[1]  # HDF5 file (e.g., mouse_cerebellar_cortex_0.400.hdf5)
simulation_name = sys.argv[2]  # Simulation type (e.g., basal_activity)
output_folder = sys.argv[3]  # Output directory

# Load the network
network = from_storage(hdf5_filename)
print(f"📂 Network loaded with configuration: {network.configuration}")

# Extract the basename (e.g., mouse_cerebellar_cortex_0.400.hdf5 -> mouse_cerebellar_cortex_0.400)
base_name = os.path.basename(hdf5_filename).replace(".hdf5", "")

# Try to extract the 'k' part if present (after the last "_")
parts = base_name.split("_")
if len(parts) >= 2 and parts[-1].replace(".", "").isdigit():
    k_suffix = parts[-1]  # Example: "0.400"
    pdf_filename = f"bsb_report_function_{k_suffix}.pdf"
else:
    pdf_filename = "bsb_report_function.pdf"

# Full output path
pdfPath = os.path.join(output_folder, pdf_filename)

try:
    # Generate report
    report = BasicSimulationReport(network, simulation_name, output_folder)
    report.print_report(pdfPath)

    print(f"✅ Report generated and saved to {pdfPath}")

except Exception as e:
    print(f"❌ An error occurred: {e}")
