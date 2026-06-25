#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="corromics-main"
PYTHON_VERSION="3.10"
MINIFORGE_DIR="$HOME/miniforge3"

echo "Setting up Corromics with joint-RPCA/Gemelli for WSL/Linux..."

if ! command -v wget >/dev/null 2>&1; then
  echo "Installing wget..."
  sudo apt update
  sudo apt install -y wget
fi

if [ ! -x "$MINIFORGE_DIR/bin/conda" ]; then
  echo "Installing Miniforge into $MINIFORGE_DIR..."
  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O "$HOME/miniforge.sh"
  bash "$HOME/miniforge.sh" -b -p "$MINIFORGE_DIR"
  rm "$HOME/miniforge.sh"
fi

source "$MINIFORGE_DIR/etc/profile.d/conda.sh"
conda config --set channel_priority strict

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Conda environment '$ENV_NAME' already exists. Updating it..."
else
  echo "Creating conda environment '$ENV_NAME'..."
  conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    python="$PYTHON_VERSION" \
    pip \
    "numpy<2" \
    scipy \
    pandas \
    scikit-learn \
    scikit-bio=0.5.9 \
    "iow<1.0.8" \
    biom-format \
    h5py \
    click \
    nose
fi

conda activate "$ENV_NAME"

echo "Installing Corromics app dependencies..."
python -m pip install --upgrade pip
python -m pip install \
  "numpy<2" \
  streamlit==1.44.0 \
  streamlit-extras \
  scipy==1.15.2 \
  statsmodels==0.14.4 \
  pandas \
  pandas_flavor \
  plotly==6.0.1 \
  openpyxl==3.1.5 \
  psutil==7.0.0 \
  networkx==3.4.2

echo "Installing Corromics optional packages..."
python -m pip install git+https://github.com/Wang-Bioinformatics-Lab/GNPSDataPackage.git

echo "Installing Gemelli..."
python -m pip install gemelli==0.0.12

conda env config vars set CORROMICS_GEMELLI_BACKEND=local

echo
echo "Setup complete."
echo "Start Corromics from this WSL/Linux shell with:"
echo "  conda activate $ENV_NAME"
echo "  streamlit run Home.py --server.port 5000 --server.address 127.0.0.1"
