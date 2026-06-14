#!/usr/bin/env bash
set -euo pipefail

CORROMICS_HOME="${HOME}/.corromics"
MINIFORGE_DIR="${CORROMICS_HOME}/miniforge"
ENV_NAME="gemelli-standalone"
ENV_FILE="${CORROMICS_HOME}/environment-gemelli-worker.yml"

mkdir -p "${CORROMICS_HOME}"

if ! command -v curl >/dev/null 2>&1; then
  sudo apt-get update
  sudo apt-get install -y curl bzip2 ca-certificates
fi

if [ ! -x "${MINIFORGE_DIR}/bin/conda" ]; then
  echo "Installing Miniforge into ${MINIFORGE_DIR}"
  INSTALLER="${CORROMICS_HOME}/Miniforge3-Linux-x86_64.sh"
  curl -L "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" -o "${INSTALLER}"
  bash "${INSTALLER}" -b -p "${MINIFORGE_DIR}"
  rm -f "${INSTALLER}"
fi

cat > "${ENV_FILE}" <<'YAML'
name: gemelli-standalone
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.10
  - numpy<2
  - scikit-bio=0.5.9
  - iow<1.0.8
  - biom-format
  - pandas
  - scipy
  - scikit-learn
  - h5py
  - click
  - pip
  - pip:
      - gemelli==0.0.12
YAML

source "${MINIFORGE_DIR}/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "Updating ${ENV_NAME}"
  conda env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
else
  echo "Creating ${ENV_NAME}"
  conda env create -f "${ENV_FILE}"
fi

conda clean -afy

WORKER_PYTHON="${MINIFORGE_DIR}/envs/${ENV_NAME}/bin/python"
"${WORKER_PYTHON}" - <<'PY'
import biom
import gemelli
from gemelli.rpca import joint_rpca
print("Gemelli worker import check passed")
PY

echo "${WORKER_PYTHON}" > "${CORROMICS_HOME}/gemelli_worker_python.txt"
echo "Worker Python: ${WORKER_PYTHON}"
