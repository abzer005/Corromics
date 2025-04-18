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
st.subheader('Input File Requirements')
st.write(""" 
         The accepted formats for the manually uploading the input files are **CSV**, **TXT**, or **XLSX**. The necessary files include:
         1. Feature (Metabolite) Quantification Table – Contains metabolite abundances. Input files must include the `.mzML` extension in their filenames.
         2. Proteomics/Metagenomics Quantification Table – Includes proteomics (protein abundances) or metagenomics (microbial composition, e.g., ASV sequences, taxa abundances).
         3. Metadata
    
         """)

st.write("""
### Data Preparation Essentials
#### 1. Feature Table
- Both your Metabolomics and Proteomics/Genomics feature tables must include a column called: `feature_ID`
- If your table contains a column like `row ID`, `index`, or other identifiers, rename it to `feature_ID`.
This applies to both Omics 1 and Omics 2 tables.
         
Example feature table:  
 
|feature_ID|sample1.mzML|sample2.mzML|blank.mzML|
|---|---|---|---|
|1|1000|1100|100|
|2|2000|2200|200|
""")

st.markdown("""          

#### 2. Metadata
##### Metadata **must include** the following columns: 
- `filename` → Filenames of the metabolite quantification table (**must match .mzML files**).  
- `ATTRIBUTE_Corromics_filename` → Filenames of proteomics/metagenomics data (**e.g., fasta or other omics files**). 

These two columns are mandatory for linking and filtering. 

##### Metadata Optional columns:
- `Sample_type` → Defines sample conditions (**e.g., control, treatment**).  
- Other relevant experimental attributes. 
  
These additional columns can be used for filtering the samples from the metabolomics feature table in the 'Data Filter' step.""")

with st.expander("About ATTRIBUTE_taxa (Optional Hierarchical Info)"):
    st.info(
        "If your **second omics table** contains **hierarchical or taxonomic columns** (e.g., Phylum, Genus, Species), "
        " you can specify them using the `ATTRIBUTE_taxa` column in your metadata table.\n\n"
        "- These should be **column names**, not filenames or sample names.\n"
        "- This field is **independent** of the `filename` column.\n"
        "- You do **not** need to fill it for all rows — just list the names once in any row.\n"
        "- For example, if your omics table has columns `Phylum`, `Genus`, `Species`, you can write:\n"
        "`Phylum, Genus, Species`\n\n"
        "- This helps the app detect and organize hierarchical groupings. Check out the **example dataset** in the app for reference.")

    st.markdown("""
                Example metadata table:
                |filename|ATTRIBUTE_Corromics_filename|Sample_Type|Time_Point|ATTRIBUTE_taxa|
                |---|---|---|---|---|
                |sample1.mzML|sample1.fasta|Sample|1h|Phylum|
                |sample2.mzML|sample2.fasta|Sample|2h|Genus|
                |sample3.mzML|sample3.fasta|Sample|3h| |
                |blank.mzML|---|Blank|N/A|---|
                """)
    
st.write(' ')

st.markdown("""        
Example metadata table:
            
|filename|ATTRIBUTE_Corromics_filename|Sample_Type|Time_Point|
|---|---|---|---|
|sample1.mzML|sample1.fasta|Sample|1h|
|sample2.mzML|sample2.fasta|Sample|2h|
|sample3.mzML|sample3.fasta|Sample|3h| 
|blank.mzML|---|Blank|N/A|
""")

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

with st.expander("Download the CorrOmics app"):
    st.write("""            
    **For Windows users:**
    - **You can directly download the Windows executable (.exe) from our lab's website [www.functional-metabolomics.com/resources](https://www.functional-metabolomics.com/resources)**
    - Click the **Download** button next to CorrOmics.
    - Run the installer and follow the on-screen instructions.

    **For macOS users:**
    Please follow the installation instructions in the [CorrOmics GitHub](https://github.com/abzer005/Corromics) repository

    """)

st.write('**If you use CorrOmics in your research, please cite:**')
st.markdown("""
            * [FBMN-STATS](https://fbmn-statsguide.gnps2.org/) - A statistical pipeline for downstream processing of FBMN results.
            * Pakkir Shah, A.K., Walter, A., Ottosson, F. et al. Statistical analysis of feature-based molecular networking results from non-targeted metabolomics data. Nat Protoc (2024). https://doi.org/10.1038/s41596-024-01046-3
            """
            )
            

# Add more links as needed

# Feedback Section
st.subheader("We Value Your Feedback")
st.markdown("""
            We welcome your feedback and suggestions to improve CorrOmics. Please feel free to [Create an Issue on GitHub](https://github.com/abzer005/Corromics/issues/new) repository to share your thoughts or report any issues you encounter. 
            Your input is invaluable in making the tool better for everyone.
""")

# Contribution and Follow Us
st.subheader("Contribute and Follow Us")
st.markdown("""
- Interested in contributing? Check out the [GitHub page](https://github.com/abzer005/Corromics).
- For more about our work, visit our [lab's GitHub page](https://github.com/Functional-Metabolomics-Lab).
""")

# Optional: Footer
st.markdown("---")
st.text("CorrOmics © 2025")