# Changelog

## v1.0.0 - Proposed

### Added
- Single-source application version in `src/version.py`.
- SparCC-compatible second-omics transformation: `Relative abundance (+1 pseudocount)`.
- OS-aware Joint-RPCA/Gemelli worker dependency guidance.
- Explicit hosted/restricted versus local deployment-mode handling.

### Changed
- SparCC outer iterations are fixed at 1 with deterministic normalization.
- Target and decoy dataframes are displayed in tabs.
- Target and decoy association result tables are displayed in tabs.
- Association result tables use scientific-notation display while preserving numeric data for exports.
- Homepage and README terminology now use association for the general workflow.
- Cytoscape GraphML guidance documents the actual exported node and edge attributes.
- README installation guidance is organized around hosted, Git clone, WSL2, Docker, Windows installer, and package routes.

### Fixed
- Docker-hosted deployments can be explicitly restricted instead of being inferred as local only because they run in Docker.
- `GEMELLI_WORKER_PYTHON` is validated and honored when launching the Joint-RPCA worker.
- Native Windows Joint-RPCA guidance no longer implies Gemelli can be installed into the native executable workflow.

### Documentation
- Added Git installation guidance for macOS, Ubuntu/Debian Linux, and Windows.
- Added WSL repeat-start instructions so users do not rerun setup each time.
- Added `conda: command not found` troubleshooting.
- Added macOS Gemelli worker setup guidance.

### Known Limitations
- Joint-RPCA remains unavailable in the native Windows executable unless a future implementation delegates to a WSL2 or Docker worker.
- SparCC and Joint-RPCA do not provide calibrated p-values or target-decoy FDR filtering in the Streamlit app.
- Streamlit displays numeric result cells in scientific notation, while full precision remains in the underlying dataframe and CSV export.
