![banner](assets/corromics_full_logo.png)

Corromics is a web-based application designed to uncover correlation patterns between microbiome and metabolome data. By integrating microbial (e.g., ASV sequences) or other omics information with metabolite information, Corromics constructs a network where each microbial entity and metabolite is represented as a node. Edges are drawn between nodes to visualize their interactions, helping users discern meaningful relationships between microbial communities and metabolic profiles. This tool facilitates a deeper understanding of microbial-metabolite associations, enabling insights into functional dynamics in diverse biological systems.

# Corromics Installation Guide

## Why Install Corromics Locally?

Running Corromics in the cloud has a built-in safeguard: a **limit of 1 million correlations** to prevent crashes. For larger datasets, we **strongly recommend installing and running Corromics locally**, allowing faster performance and access to more computing power directly on your machine.

---

## 🖥️ 1. Windows Users

You can directly download the **Windows executable (.exe)** from our lab's website:

🔗 [www.functional-metabolomics.com/resources](https://www.functional-metabolomics.com/resources)

- Click the **Download** button next to *Corromics*.
- Run the installer and follow the on-screen instructions.

---

## 🍎 2. macOS Users

Due to macOS security restrictions, we do **not** provide a pre-built app.  
However, Mac users can easily **build and run Corromics** locally using a simple shell script.

### What the Shell Script Does:

- Checks for **Python** and **Git** (installs via Homebrew if missing)
- Installs **Homebrew** if needed
- Creates a **Conda virtual environment** (`corromics_env`) and activates it
- Installs all required Python dependencies
- Clones the **Corromics GitHub repository**
- Launches the app at `http://localhost:8502`

---

## Installation Steps for Mac

### ✅ Step 1: Download the Installation Script

🔗 [Download `install_corromics.sh`](https://www.functional-metabolomics.com/resources)

This file will typically appear in your **Downloads** folder.

---

### ✅ Step 2: Open Terminal and Run the Script

```bash
# Navigate to your Downloads folder
cd ~/Downloads

# Make the script executable
chmod +x install_corromics.sh

# Run the installer script
./install_corromics.sh

```

### ✅  Step 3: Launching the App
Once installation is complete, you will see the following message in your terminal:

```
"browser.browserName" is not a valid config option. If you previously had this config option set, it may have been removed.
  
You can now view your Streamlit app in your browser.
  
URL: http://localhost:8502

```

The app should open automatically in your default web browser. If it doesn’t, manually paste http://localhost:8502 into your browser.


### Running Corromics After Installation
Once installed, you don’t need to run the script again. To start Corromics in the future:

Open Terminal and enter the following commands:

```bash
cd ~/Corromics 
source corromics_env/bin/activate 
streamlit run Home.py --server.port 8502 --server.address localhost
```

- Navigate to Corromics directory
- Activate the virtual environment
- Launch the app


### Closing the App
To stop Corromics, press **Ctrl + C** in the terminal.
Even if you close the browser, the app keeps running in the background until you stop it from the terminal.


