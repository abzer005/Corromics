import streamlit as st
from src.common import *
import pandas as pd
from streamlit.components.v1 import html

page_setup()

initialize_app()
st.title('CORROMICS')

c1, c2, c3 = st.columns([2, 2, 2])  # Adjust column ratios for centering
with c2:
    try:
        st.image("assets/corromics_icon.png", caption="Corromics Logo", use_container_width=True)
    except TypeError:
        st.image("assets/corromics_icon.png", caption="Corromics Logo", use_column_width=True)

# Introduction
st.markdown("""
<div style="color: #05577C; font-size: 25px;">
<b>Welcome to Corromics!</b> Before proceeding with your analysis, please take a moment to read this homepage.
</div>
""", unsafe_allow_html=True)

st.write(' ')
st.subheader('CorrOmics: A Multi-Omics Correlation Analysis Tool')

st.write("""
         CorrOmics is a web-based application designed to uncover correlation patterns between microbiome and metabolome data. 
         By integrating microbial (e.g., ASV sequences) and metabolite information, CorrOmics constructs a network where each microbial entity and metabolite is represented as a node. 
         Edges are drawn between nodes to visualize their interactions, helping researchers discern meaningful relationships between microbial communities and metabolic profiles. 
         This tool facilitates a deeper understanding of microbial-metabolite associations, enabling insights into functional dynamics in diverse biological systems.""")

# Input Files
st.markdown("""
### Data Preparation Essentials
            
Corromics supports **multiple ways** of providing your data. The app will do as much
automatic restructuring as possible depending on the chosen input method.
            """)
            
with st.expander("🧾 1. Metadata — Linking Omics Datasets", expanded=False):
     st.markdown(
           """
The metadata table is **critical** for linking metabolomics and other omics.
            
**Required columns:**
- `filename`  
    - Filenames of the metabolomics feature table, should have the extensions as well (e.g. `sample1.mzML`).  
    - These must match the sample names used in the metabolomics table  
        (or their basenames, e.g. `spiked_1`, `spiked_1.mzML`).
- `ATTRIBUTE_Corromics_filename`  
    - Filenames used for the other omics quantification (e.g. FASTA/ASV/proteomics IDs).  
    - Used to map between metabolomics samples and corresponding other omics sample

The app internally uses:
- `filename` as the **metadata index** (for matching to metabolomics), and  
- `ATTRIBUTE_Corromics_filename` to map to the **Other Omics Feature Table** column

**Optional (but recommended) metadata columns:**
- Experimental grouping, e.g.:`Sample_type` (control, treatment, blank, QC, etc.)
- Any other relevant attributes (timepoint, tissue, condition, batch, etc)
These optional columns are used in the **Data Filter** to subset data for correlation analysis

Example metadata table:
            
|filename|ATTRIBUTE_Corromics_filename|Sample_Type|Time_Point|
|---|---|---|---|
|sample1.mzML|sample1.fasta|Sample|1h|
|sample2.mzML|sample2.fasta|Sample|2h|
|sample3.mzML|sample3.fasta|Sample|3h| 
|blank.mzML|---|Blank|N/A|
          """
     )

with st.expander("🧬 2. Other Omics Feature Table (Proteomics / Genomics / ASV)", expanded=False):
     st.markdown(
           """
- Upload your **other omics quantification table** (ASV sequencing, proteomics, etc.).  
- The **first column** is treated as the feature identifier, automatically renamed to `feature_ID` and used as the index.
- All remaining columns are treated as **sample columns** or **metadata/taxonomy columns**. The app will:
    - align sample columns with the metadata, and  
    - detect non-sample columns (taxonomy / annotation) automatically for later filtering and ordering.
""")

with st.expander("⚛️ 3. Metabolomics Feature Table", expanded=False):
    st.markdown(
           """
You can provide the metabolomics table via:

- **Use Example Dataset**  
  No preparation needed – preformatted demo data for exploration.

- **Manual Input (Custom Data)**  
  - Expects features in **rows** and samples in **columns**.  
  - Your table **must** contain a column named `feature_ID`.  
  - `feature_ID` will be used as the row index in the app.  
  - Make sure the filename extensions (eg. mzML) match with those in the metadata.  

  Example (manual input):  

  |feature_ID|sample1.mzML|sample2.mzML|blank.mzML|
  |---|---|---|---|
  |1|1000|1100|100|
  |2|2000|2200|200|

- **FBMN-Stats Output (Metabolomics Table)**  
  - Upload the feature table exported from **FBMN-Stats** as it is.  
  - Expected format: **samples as rows**, **features as columns**.  
  - The first (unnamed) column is treated as sample IDs.  
  - The app automatically detects feature columns (e.g., `ID_mz_RT` or `ID_mz_RT&name`),  
    transposes the table, and matches filenames with metadata (e.g., `spiked_1` ↔ `spiked_1.mzML`).

- **MZmine Export (Metabolomics Table)**  
  - Upload the MZmine feature table containing the columns `row ID`, `row m/z`, and `row retention time`.  
  - The app automatically combines these into a unique `feature_ID` (e.g., `45_233.12_5.02`).  
  - If the metadata filenames do not include the text *Peak area*,  
    the app automatically removes that suffix from column names for matching.
        """
    )

# Output Files
st.subheader('Output File Information')
st.write("""
- Upon processing your data, CorrOmics generates an output **edge file** in CSV format. 
- You can download the result as a **GraphML** file as well for Cytoscape visualization.  
                  """)

with st.expander("💡 Tips for Using the GraphML file in Cytoscape"):
    st.markdown("""
    - Simply **drag and drop** the `.graphml` file into Cytoscape.
    - Once loaded, you will see both the **edge table** and the **node table**
    ##### Node Table
    - **Node names** are derived from the `shared_name` column by default.
    - The `Node_Info` column contains detailed annotations:
    - For the example dataset:
        - **Omics_1** refers to features from the metabolomics quantification table (the first table you uploaded).
        - **Omics_2** refers to the binned ASVs from the ASV table.
    - You can use `Node_Info` to **color nodes** or **assign different shapes** based on node type in Cytoscape.
    - The `Original_index` column contains the original identifiers:
    - For metabolomics features: e.g., `ID_mz_RT`.
    - For ASVs: longer binned names like `Bacteria_Proteobacteria_Gammaproteobacteria_Thiotrichales_Thiotrichaceae`
    
    ##### Edge Table
    - Use the `Absolute_Correlation_Score` column to **assign weights** to edges in your network.
    - Use the `Sign_Score` column (`-1` or `+1`) to **color positive and negative correlations** for better visual distinction.
                 """)

# Subheader and Interactive Features
st.subheader('About the App Elements')
st.markdown("""
🔍 **How to Know If the App Is Running?**  
If you're performing a calculation or generating a figure and don't see any updates on the main page, 
check the **top-right corner** of the page. If you see the message **'RUNNING'**, the app is active and processing your request.  
            
💡 **All plots are interactive!**  
- Use your mouse to select areas and zoom in on specific regions.  
- Double-click on the plot to reset the zoom and return to the default view.  
- Save plots using the **camera icon** in the top-right corner of the plot. You can specify the image format (e.g., PNG) in the settings panel before saving.
""")

# Citation and Resources
st.subheader('Citation and Further Resources')

st.warning("""
           Running **CorrOmics in the cloud** at [https://corromics.gnps2.org/](https://corromics.gnps2.org/) comes with a **restriction for over 1 million correlations**, to protect server performance. For larger datasets, we recommend running CorrOmics locally. 
           """)

with st.expander("For Windows Users", expanded=False):
    st.markdown(
        """
Corromics provides a **standalone installer** for Windows systems.

**How to use:**
1. Download the latest `.exe` installer from our lab's website [www.functional-metabolomics.com/resources](https://www.functional-metabolomics.com/resources)**.
2. Click the **Download** button next to CorrOmics.  
3. Run the installer and follow the on-screen instructions.  
4. Once installed, the app can be accessed directly like any normal desktop application.  

**Tip:**  
If you prefer full control or encounter compatibility issues,  
you can also run the app locally by cloning the [CorrOmics repository](https://github.com/Functional-Metabolomics-Lab/Corromics) from GitHub and using:
```bash
pip install -r requirements.txt
streamlit run Home.py
```

    """)

with st.expander("For macOS Users", expanded=False):
    st.markdown(
    """
There is no standalone desktop app for macOS at the moment.

You can run Online (Recommended for smaller datasets). If your analysis involves fewer than 1 million correlations,
you can safely use the hosted web version of the app.

### Run Locally (Recommended for Larger Datasets)
For heavy analyses or larger datasets, it’s best to run Corromics locally:

1. **Clone the repository**
   ```bash
   git clone https://github.com/Functional-Metabolomics-Lab/Corromics.git
   cd Corromics
   ````

2. **Install dependencies and launch the app**.
Make sure you have Python 3.11 installed (same version used in the Windows .exe build).
```bash
pip install -r requirements.txt
streamlit run Home.py
````
Running locally avoids browser memory limits and gives full control over computation.
""")

#st.write('**If you use CorrOmics in your research, please cite:**')
#st.markdown("""
#            * [FBMN-STATS](https://fbmn-statsguide.gnps2.org/) - A statistical pipeline for downstream processing of FBMN results.
#            * Pakkir Shah, A.K., Walter, A., Ottosson, F. et al. Statistical analysis of feature-based molecular networking results from non-targeted metabolomics data. Nat Protoc (2024). https://doi.org/10.1038/s41596-024-01046-3
#            """
#            )
            

# Add more links as needed

# Feedback Section
st.subheader("We Value Your Feedback")
st.markdown("""
            We welcome your feedback and suggestions to improve CorrOmics. Please feel free to [Create an Issue on GitHub](https://github.com/Functional-Metabolomics-Lab/Corromics/issues/new) repository to share your thoughts or report any issues you encounter. 
            Your input is invaluable in making the tool better for everyone.
""")

# Contribution and Follow Us
st.subheader("Contribute and Follow Us")
st.markdown("""
- Interested in contributing? Check out the [GitHub page](https://github.com/Functional-Metabolomics-Lab/Corromics).
- For more about our work, visit our [lab's GitHub page](https://github.com/Functional-Metabolomics-Lab).
""")

# Optional: Footer
st.markdown("---")
st.text("CorrOmics © 2025")
