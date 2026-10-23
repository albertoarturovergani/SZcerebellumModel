# Validation scope

The 2026 replay completed 96/96 historical case/level simulations with expected HDF5
and NIO outputs. Direct HDF5 structural comparison was performed on 64/64 of the
originally selected Auto-Rhythm and Kase reconstructions. Paired raw Basal 4 Hz
firing-rate estimates were available for 31/32 old-vs-new case/level pairs.

A missing original historical NIO group cannot be used as a new result and must be
marked `historical_output_unavailable` in quantitative comparisons. Source scripts
in `scripts_py/` retain their historical analysis assumptions; not all panels
have been independently re-executed from a freshly cloned public repository.

The outputs of the validation scripts are not the original Zenodo data archive.
See https://zenodo.org/records/23391996 for source data and its provenance.
