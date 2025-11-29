from cerebellum.analysis.structure_analysis import StructureReport
from bsb import from_storage

scaffold = from_storage('./test1/mouse_cerebellum.hdf5')
report = StructureReport(scaffold)

report.print_report('./test1/report_repo_3.pdf')




