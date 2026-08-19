# CorrOmics Repository Audit

This audit supports the restructuring pass described in `CORROMICS_RESTRUCTURE_SPEC.md`.
No files were deleted during this pass.

## GNPS Template Comparison

The upstream Wang Bioinformatics Lab `Template_Streamlit_App` contains the GNPS-relevant layout conventions used here: `Dockerfile`, `Makefile`, `docker-compose.yml`, `docker-compose-dev.yml`, `requirements.txt`, `run_server.sh`, a Streamlit homepage entry point, and `pages/`.

CorrOmics intentionally differs from the template by adding app-specific scientific code under `src/`, vendored SparCC code, example datasets, Docker Gemelli worker setup, WSL Joint-RPCA setup, Windows executable packaging files, documentation, and assets. These differences are functional and should be preserved unless replaced by an explicitly approved deployment/build plan.

## File Categories

| Path | Category | Notes |
|---|---|---|
| `Home.py` | app runtime | Main Streamlit homepage. |
| `pages/` | app runtime | Streamlit workflow pages. |
| `src/correlation.py`, `src/dist_correlation.py`, `src/sparcc.py`, `src/joint_rpca.py`, `src/joint_rpca_worker.py`, `src/fdr.py` | analysis/scientific code | Core association methods and worker bridge. |
| `src/common.py`, `src/fileselection.py`, `src/progress.py`, `src/visualization.py`, `src/molnet.py`, `src/correlation_sections.py` | app runtime | Shared UI, import, export, and visualization helpers. |
| `src/vendor/sparcc3/` | third-party/vendor | Vendored SparCC implementation and upstream README. Preserve notices. |
| `Dockerfile`, `docker-compose.yml`, `docker-compose-dev.yml`, `run_server.sh`, `Makefile`, `requirements.txt`, `environment-gemelli-worker.yml` | hosted GNPS/Docker deployment and local Docker | Template-critical plus Gemelli worker environment. |
| `Docker_User_Guide/` | docs | Local Docker Desktop guidance. |
| `WSL_Joint_RPCA/` | Windows/WSL | WSL setup for Gemelli/Joint-RPCA. |
| `.github/workflows/corromics-build-windows-executable.yaml`, `run.py`, `run_app_temp.spec`, `hooks/hook-streamlit.py`, `assets/corromics_icon.ico`, `assets/corromics_license.rtf` | Windows packaging/build | Preserve until executable build is independently verified without them. |
| `assets/` | example data/assets and packaging assets | Logos/icons/license text used by app and packaging. |
| `example-data/` | example data/assets | Existing example datasets. |
| `.streamlit/config.toml`, `.streamlit/secrets.toml` | deployment/configuration | Streamlit runtime configuration. `secrets.toml` should be reviewed for deployability and sensitivity. |
| `.devcontainer/` | development/configuration | Developer environment support. |
| `pyproject.toml`, `requirements.lock.txt`, `environment.lock.yml`, `install_corromics.sh` | development/build | Packaging metadata, lockfiles, and local install helper. |
| `README.md`, `CHANGELOG.md`, `README_JOINT_RPCA.txt`, `REPOSITORY_AUDIT.md` | docs | User/deployment/release/audit documentation. |
| `AGENTS.md`, `CORROMICS_RESTRUCTURE_SPEC.md` | local Codex/spec files | Untracked local planning files in this clone. |
| `__pycache__/`, `pages/__pycache__/`, `src/__pycache__/`, `src/vendor/sparcc3/__pycache__/`, `.DS_Store`, `build_arm64.log` | generated/temporary | Candidate cleanup items. |

## Candidate Removal Or Ignore List

Do not remove these until approved:

- `.DS_Store`: macOS Finder artifact.
- `build_arm64.log`: generated build log; keep only if needed for debugging a current Apple Silicon build issue.
- `__pycache__/`, `pages/__pycache__/`, `src/__pycache__/`, `src/vendor/sparcc3/__pycache__/`: generated Python bytecode caches.
- `README_JOINT_RPCA.txt`: appears partially redundant with `README.md` and `WSL_Joint_RPCA/README.md`, but may still be a legacy user-facing note.
- `src/helpers.py`: appears unreferenced by current imports and may be obsolete; needs manual review before deletion because it contains older Streamlit workflow code.
- `.streamlit/secrets.toml`: review whether this should remain local-only for deployments. Do not remove if GNPS or Streamlit configuration depends on it.

## Intentionally Kept

- GNPS template-shaped deployment files: required for hosted deployment compatibility.
- Docker Gemelli files: required for local Docker and hosted Docker Joint-RPCA support.
- WSL setup files: required for Windows Joint-RPCA route.
- Windows executable files: required to preserve the native Windows installer workflow.
- Vendored SparCC files: required by current SparCC implementation and license/notice preservation.
- Existing example datasets and assets: required for demo workflows and app branding.

## Licensing Audit

- CorrOmics-owned code appears to be the Streamlit app, local `src/` modules outside `src/vendor/`, documentation, deployment scripts, and packaging files.
- Vendored third-party code is present under `src/vendor/sparcc3/`, including `ReadMe.md` with upstream attribution to SparCC and its Python 3 port.
- The current top-level README states BSD 3-Clause and points to `assets/corromics_license.rtf`; there is no top-level `LICENSE` file in this clone.
- Recommended final layout, pending maintainer approval: add a top-level `LICENSE` for CorrOmics BSD 3-Clause and keep third-party notices/licenses with vendored code or in a dedicated notices file.
- No license text was replaced or removed in this pass.
