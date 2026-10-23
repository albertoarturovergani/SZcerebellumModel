SZ Numerical Regression v1
==========================

1) Place regression_SZ_smoke_v1.py inside validation/ of SZcerebellumModel-TEST.
2) Capture the previously validated full-run CSV files as provisional baselines.
3) Run the three validators WITHOUT --dry-run into NEW, distinct output directories.
4) Compare the new CSVs with the snapshot using explicit tolerances.

Features: Python standard library only. Checks headers, row counts, nonnumeric
cells, every numeric cell using abs/rel tolerances, CSV/JSON readability, and
SHA256 integrity of baseline snapshots. Exits 0 (PASS), 1 (FAIL), 2 (error).
Checks compare values in the same row order. Row order changes count as a
regression. Check JSON files are parsed and snapshotted for audit, but JSON
numerical fields are not compared in v1 because they may contain metadata.
The script deliberately rejects a comparison with the original capture files.

The baseline captured from existing smoke-test outputs is not independent
validation of the original scientific calculations or manuscripts.

Do not commit the baseline until comparing a fresh run and reviewing values.
No simulations or notebooks are modified by this runner.
