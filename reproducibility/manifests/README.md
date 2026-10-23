# Provenance manifest

`source_files.json`: for every copied input YAML/SWC, SHA256 of the local
2026 replay input and SHA256 of the portable GitHub input; scientific settings
are not intentionally changed. `changes` records portability substitutions.

`SHA256SUMS.txt`: checksums of the staged repository files, excluding this
checksum file itself. Historical originals on Zenodo may retain the outdated
upstream package-version metadata and local absolute reference paths.
