import streamlit as st
from src.common import *
import pandas as pd
from streamlit.components.v1 import html

page_setup()

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
st.subheader('Corromics: A Multi-Omics Correlation Analysis Tool')

st.write("""
         Corromics is a web-based application designed to uncover correlation patterns between microbiome and metabolome data. 
         By integrating microbial (e.g., ASV sequences) and metabolite information, Corromics constructs a network where each microbial entity and metabolite is represented as a node. 
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
         
**Metadata must include the following columns:**  
- `filename` → Filenames of the metabolite quantification table (**must match .mzML files**).  
- `ATTRIBUTE_Corromics_filename` → Filenames of proteomics/metagenomics data (**e.g., fasta or other omics files**).  

**Metadata can include additional columns such as:**  
- `Replicates` → Identifies biological/technical replicates.  
- `Sample_type` → Defines sample conditions (**e.g., control, treatment**).  
- Other relevant experimental attributes.  
""")

st.markdown("""          
Example feature table:  
 
|feature_ID|sample1.mzML|sample2.mzML|blank.mzML|
|---|---|---|---|
|1|1000|1100|100|
|2|2000|2200|200|
""")

st.write(' ')

st.markdown("""        
Example metadata table:
            
|filename|ATTRIBUTE_Corromics_filename|Sample_Type|Time_Point|
|---|---|---|---|
|sample1.mzML|sample1.fasta|Sample|1h|
|sample2.mzML|sample2.fasta|Sample|2h|
|blank.mzML|---|---|Blank|N/A| 
""")

# Output Files
st.subheader('Output File Information')
st.write("""
- Upon processing your data, Corromics generates an edge file as an output in CSV format. 
- You can download the result as a **GraphML** file as well for Cytoscape visualization.  
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
st.write('If you use Corromics in your research, please cite:')
st.markdown("""
            * [FBMN-STATS](https://fbmn-statsguide.gnps2.org/) - A statistical pipeline for downstream processing of FBMN results.
            * Pakkir Shah, A.K., Walter, A., Ottosson, F. et al. Statistical analysis of feature-based molecular networking results from non-targeted metabolomics data. Nat Protoc (2024). https://doi.org/10.1038/s41596-024-01046-3
            """
            )
            

# Add more links as needed

# Feedback Section
st.subheader("We Value Your Feedback")
st.markdown("""
            We welcome your feedback and suggestions to improve Corromics. Please feel free to create an issue on our GitHub repository to share your thoughts or report any issues you encounter. 
            Your input is invaluable in making the tool better for everyone.

            [Create an Issue on GitHub](https://github.com/abzer005/Corromics/issues/new)
""")

# Contribution and Follow Us
st.subheader("Contribute and Follow Us")
st.markdown("""
- Interested in contributing? Check out the [GitHub page](https://github.com/abzer005/Corromics).
- For more about our work, visit our [lab's GitHub page](https://github.com/Functional-Metabolomics-Lab).
""")

# Optional: Footer
st.markdown("---")
st.text("Corromics © 2025")