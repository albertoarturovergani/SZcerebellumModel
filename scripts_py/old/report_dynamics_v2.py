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
hdf5_filename = sys.argv[1]  # HDF5 file
simulation_name = sys.argv[2]  # Simulation type (e.g., basal_activity)
output_folder = sys.argv[3]  # Output directory

# Load the network
network = from_storage(hdf5_filename)
print(f"📂 Network loaded with configuration: {network.configuration}")

# Extract filename without path
pdf_filename = os.path.basename(hdf5_filename).replace(".hdf5", ".pdf")

# Create the specific dynamics subfolder
dynamics_folder = os.path.join(output_folder, f"dynamics_{pdf_filename.replace('.pdf', '')}")
os.makedirs(dynamics_folder, exist_ok=True)  # ✅ Ensure the directory exists

# Define final report path
pdfPath = os.path.join(dynamics_folder, pdf_filename)

try:
    # Generate report
    report = BasicSimulationReport(network, simulation_name, output_folder)
    report.print_report(pdfPath)

    print(f"✅ Report generated and saved to {pdfPath}")

except Exception as e:
    print(f"❌ An error occurred: {e}")
