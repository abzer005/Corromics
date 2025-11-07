#!/bin/bash

echo "Setting up Corromics on your Mac..."

# Step 1: Check if Python3 is installed
if ! command -v python3 &> /dev/null; then
    echo "Python3 is not installed. Installing via Homebrew..."
    
    # Check if Homebrew is installed
    if ! command -v brew &> /dev/null; then
        echo "Homebrew is not installed. Installing Homebrew first..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    fi
    
    # Install Python
    brew install python git
fi

# Step 2: Create a virtual environment

echo "Downloading Miniconda..."
ARCH=$(uname -m)
if [ "$ARCH" = "arm64" ]; then
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
else
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
fi
wget $MINICONDA_URL -O ~/miniconda.sh

echo "🔧 Installing Miniconda..."
bash ~/miniconda.sh -b -p $HOME/miniconda
rm ~/miniconda.sh

echo "🔁 Initializing Conda..."
source "$HOME/miniconda/etc/profile.d/conda.sh"
conda init

#Create and activate Conda environment
echo "Creating Conda environment 'corromics_env'..."
conda create -y -n corromics_env python=3.11
conda activate corromics_env

# Step 3: Install dependencies
echo "Installing required dependencies..."
pip install --upgrade pip
pip install streamlit pandas numpy scipy streamlit-extras plotly openpyxl psutil statsmodels networkx

echo "All dependencies installed successfully!"

# Step 4: Clone Corromics repository (if not already downloaded)
if [ ! -d "$HOME/Corromics" ]; then
    echo "Cloning Corromics from GitHub..."
    git clone https://github.com/Functional-Metabolomics-Lab/Corromics.git "$HOME/Corromics"
else
    echo "Corromics folder already exists. Skipping clone."
fi

# Step 5: Run Corromics
echo "Launching Corromics..."
cd $HOME/Corromics
streamlit run Home.py --server.port 8502 --server.address localhost
