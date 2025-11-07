![banner](assets/corromics_full_logo.png)

Corromics is a web-based application designed to uncover correlation patterns between microbiome and metabolome data. By integrating microbial (e.g., ASV sequences) or other omics information with metabolite information, Corromics constructs a network where each microbial entity and metabolite is represented as a node. Edges are drawn between nodes to visualize their interactions, helping users discern meaningful relationships between microbial communities and metabolic profiles. This tool facilitates a deeper understanding of microbial-metabolite associations, enabling insights into functional dynamics in diverse biological systems.

# Corromics Installation Guide

## Why Install Corromics Locally?

Running Corromics in the [cloud](https://corromics.gnps2.org/) has a built-in safeguard: a **limit of 1 million correlations** to prevent crashes. For larger datasets, we **strongly recommend installing and running Corromics locally**, allowing faster performance and access to more computing power directly on your machine.


## 1. Windows Users

You can directly download the **Windows executable (.exe)** from our lab's website:

🔗 [www.functional-metabolomics.com/resources](https://www.functional-metabolomics.com/resources)

- Click the **Download** button next to *Corromics*.
- Run the installer and follow the on-screen instructions.

---

## 2. macOS Users

We don't have a standalone desktop app for macOS. If your analysis involves fewer than 1 million correlations, you can safely use the hosted web version of the app.

Run Locally (Recommended for Larger Datasets)
For heavy analyses or larger datasets, it’s best to run Corromics locally:

- Clone the repository
````
git clone https://github.com/Functional-Metabolomics-Lab/Corromics.git
cd Corromics
````
 - Install dependencies and launch the app. Make sure you have Python 3.11 installed (same version used in the Windows .exe build).
```
pip install -r requirements.txt
streamlit run Home.py
```
Running locally avoids browser memory limits and gives full control over computation.

---

## Closing the App
To stop Corromics, press **Ctrl + C** in the terminal.
Even if you close the browser, the app keeps running in the background until you stop it from the terminal.


