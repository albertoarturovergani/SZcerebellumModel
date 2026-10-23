# Reproduction inputs and outputs

`configurations/SZ_repro_2025/experiments/` preserves each experimental YAML tree,
with local absolute morphology imports rewritten to relative paths and a stale
upstream package dependency corrected to cerebellum 0.5.0. The shared SWC inputs
are under `morphologies/SZ_repro_2025/<step>/`.

Use `scripts/replay.py` to materialize a **writable** per-case/step directory
outside this Git repository. The runner copies input YAML, five ATR SWCs,
`morphologies_<k>.yaml`, and a source snapshot `cerebellum/` to the working
folder. It rewrites only the morphology YAML `$import` path for the new location,
then runs BSB compile, reconfigure and NEST simulate.

Exact bitwise network reproduction is not expected, because the original
structural RNG state was not archived. NEST seeds ARE preserved per simulation.

Output HDF5 and NIO files belong outside the repository; the original 2025
HDF5/NIO are available on Zenodo: https://zenodo.org/records/23391996.
