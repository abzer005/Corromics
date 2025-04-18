# Import necessary libraries
import streamlit as st
import pandas as pd
from src.common import *        # Importing common functionalities
from src.fileselection import * # Importing file selection functionalities

initialize_app()

# Introduction Section
st.markdown("### Please select your method for data input below.")

# Input Selection Section
input_method = st.selectbox("Select Input Method", 
                            ["Use Example Dataset",
                             "Manual Input",
                             ],
                             index=1,
                             key="input_method",
                             )

# Clearing the session state 
if 'last_input_method' not in st.session_state:
    st.session_state['last_input_method'] = None

# Example Dataset Section
if input_method == "Use Example Dataset":
    
    # Check if input method has changed
    if st.session_state['last_input_method'] != st.session_state["input_method"]:
        st.write()
        # Clear the data
        for key in ['ft', 'omics_ft', 'md', ]:
            st.session_state[key] = None

        # Update the last input method
        st.session_state['last_input_method'] = st.session_state["input_method"]
        
        load_example()  # Load data into session state

    for file_name, key in zip(["Metabolomics Feature Table", 
                               "Proteomics/Genomics Feature Table",
                               "MetaData"],
                              ['ft', 'omics_ft', 'md']):
        if st.session_state.get(key) is not None:
            display_dataframe_with_toggle(key, file_name)# Import necessary libraries


# Manual Input Section
if input_method == "Manual Input":

    if st.session_state['last_input_method'] != input_method:
        # Clear the data
        for key in ['ft', 'omics_ft', 'md']:
            st.session_state[key] = None

        # Update the last input method
        st.session_state['last_input_method'] = input_method

    st.info("💡 Upload tables in txt (tab separated), tsv, csv or xlsx (Excel) format.")

    # Create 3 columns for the ft, md file uploaders
    col1, col2, col3 = st.columns(3)
    with col1:
        ft_file = st.file_uploader("Upload Metabolomics Feature Table", 
                                   type=["csv", "xlsx", "txt", "tsv"],
                                   key="ft_upload",
                                   help = "This table is a key output of LC-MS/MS metabolomics studies. The table presents a list of mass spectral features along with their relative intensities (represented by its integrated peak area) observed across various samples.")
        if ft_file:
            st.session_state['ft'] = load_ft(ft_file)

    with col2:
        omics_ft_file = st.file_uploader("Upload Proteomics/Genomics Feature Table", 
                                        type=["csv", "xlsx", "txt", "tsv"],
                                        key="omics_upload",
                                        help = ("This table represents the key output of proteomics or genomics studies, "
                                                "providing a list of proteins, genes, or other molecular entities along with their "
                                                "quantification (e.g., expression levels or abundances) across various samples. "
                                                "It is essential for integrating and correlating this dataset with metabolomics data."
                                                )
                                                )
        if omics_ft_file:
            st.session_state['omics_ft'] = load_omics_ft(omics_ft_file)
    
    with col3:
        md_file = st.file_uploader("Upload Metadata", 
                                   type=["csv", "xlsx", "txt", "tsv"],
                                   key="md_upload",
                                   help = "The metadata table is created by the user, providing additional context for the measured samples, such as sample type, species, and tissue type, etc.")
        if md_file:
            st.session_state['md'] = load_md(md_file)

    # Display headers and 'View all' buttons for each file
    for file_name, key in zip(["Metabolomics Feature Table", 
                               "Proteomics/Genomics Feature Table", 
                               "Metabolomics MetaData"],
                              ['ft', 'omics_ft', 'md']):
        if st.session_state.get(key) is not None:
            display_dataframe_with_toggle(key, file_name)

else:
    # If data is not available, display a message
    st.warning("Input data not loaded yet. Please load the data first.")


st.markdown("## Data Filter")
# Check if the data is available in the session state
if (
    'ft' in st.session_state and 
    'md' in st.session_state and 
    st.session_state['ft'] is not None and 
    not st.session_state['ft'].empty and 
    st.session_state['md'] is not None and 
    not st.session_state['md'].empty
):

    ft = st.session_state['ft'].copy()
    md = st.session_state['md'].copy()
    
    # Get taxa information from metadata if available
    # Check if the metadata contains the 'ATTRIBUTE_taxa' column
    if 'ATTRIBUTE_taxa' in md.columns:
        st.session_state['taxa_series'] = md['ATTRIBUTE_taxa'].dropna()
    
    # If data is available, proceed with cleanup and checks
    cleaned_ft = clean_up_ft(ft)
    cleaned_md = clean_up_md(md)

    # Check if ft column names and md row names are the same
    cleaned_md, cleaned_ft = check_columns(cleaned_md, cleaned_ft)
    
    st.markdown("### Metadata overview")
    df = inside_levels(cleaned_md)
    st.dataframe(df)

    st.session_state['ft_for_analysis'] = cleaned_ft
    st.session_state['md_for_analysis'] = cleaned_md

    ################################################################################
    st.markdown("### Filter the Metabolomics Data")
    st.markdown("#### Filter by metadata")

    with st.container():

        c1, c2 = st.columns(2)
        # Allow the user to select any column for further filtering
        filter_column = c1.selectbox(
            "**Select the metadata column for filtering the metabolomics data**",
            options=st.session_state['md_for_analysis'].columns,
            key="filter_column"
        )

        # Multi-select for categories in the selected column
        filter_group = c2.multiselect(
            "**Select categories for filtering**",
            options=sorted(st.session_state['md_for_analysis'][filter_column].dropna().unique()),
            key="filter_group"
        )

        # Apply the filter if categories are selected
        if filter_group:
 
            filter_group = list(map(str, filter_group))  # Convert the group to strings if needed for matching
            filter_indices = cleaned_md[cleaned_md[filter_column].astype(str).isin(filter_group)].index
                        
            # Update the feature table and metadata based on additional filtering
            final_ft = cleaned_ft.loc[:, filter_indices]
            final_md = cleaned_md.loc[filter_indices]
        
        else:
            st.info("ℹ️ No groups selected — continuing with the full dataset. You can filter by metadata using the options above.")
            final_ft = cleaned_ft.copy()
            final_md = cleaned_md.copy()

    # Drop fully empty rows and columns
    final_ft = final_ft.dropna(how='all').dropna(axis=1, how='all')

    # Keep only numeric columns
    final_ft = final_ft.select_dtypes(include="number").copy()

    with st.expander(f"Metabolomics feature table {final_ft.shape}"):
        st.dataframe(final_ft)

    with st.expander(f"Metadata {final_md.shape}"):
        st.dataframe(final_md)
       
    # Update session state with the final filtered tables
    st.session_state['metabolome_ft'] = final_ft
    st.session_state['metabolome_md'] = final_md
    
st.markdown("---")
st.markdown("#### Filter the Other Omics Data")

#------------------------------------------------------------------------------        
#### Filter for the Genomics Data
if (
    'omics_ft' in st.session_state and 
    'metabolome_md' in st.session_state and 
    st.session_state['omics_ft'] is not None and 
    not st.session_state['omics_ft'].empty and 
    st.session_state['metabolome_md'] is not None and 
    not st.session_state['metabolome_md'].empty
):

    omics_ft = st.session_state['omics_ft'].copy()
    omics_md = st.session_state['metabolome_md'].copy()
    metabolome_ft = st.session_state['metabolome_ft'].copy()

    taxonomy_columns = []

    if 'taxa_series' in st.session_state:
        taxa_series = st.session_state['taxa_series']
        
        # Ensure column names are strings
        omics_ft.columns = omics_ft.columns.astype(str)

        # Extract valid taxonomy columns
        taxonomy_columns = sorted(set(
            col.strip() for col in taxa_series if col.strip() in omics_ft.columns
        ))

    # Create mapping from 'Corromics_filename' to 'filename'
    metadata = omics_md[['ATTRIBUTE_Corromics_filename']].dropna()
    corromics_to_filename = dict(zip(metadata['ATTRIBUTE_Corromics_filename'], metadata.index))

    # Rename columns in omics_df using that mapping
    omics_ft_renamed = omics_ft.rename(columns=corromics_to_filename)

    # Find shared columns (i.e., samples with both data)
    common_samples = set(corromics_to_filename.values()) & set(metabolome_ft.columns) & set(omics_ft_renamed.columns)

    # Subset both DataFrames to only shared samples
    st.session_state['metabolome_ft'] = metabolome_ft[sorted(common_samples)]
    omics_ft = omics_ft_renamed.copy()

    valid_sample_columns = set(common_samples)
    st.session_state["valid_sample_columns"] = valid_sample_columns
    
    # uncomment the following if you wanna see all the other columns in addition to taxa columns
    # relevant_columns = [
    #     col for col in omics_ft.columns 
    #     if col in valid_sample_columns or col not in corromics_to_filename.keys()
    # ]

    # Subset omics_ft to only relevant columns (sample + metadata columns)
    relevant_columns = taxonomy_columns.copy()
    omics_ft_subset = omics_ft[relevant_columns]
    ordered_columns = order_taxonomic_columns(relevant_columns)

    ###########-------------------------------------------------
    # Ask if the table contains taxonomic information

    has_taxonomic_info = st.radio("Does this quantification table contain any hierarchical or structured metadata (e.g., groupings like Class, Type)?", 
                                  ["Yes", "No"],
                                  index=1
                                  )

    if has_taxonomic_info == "Yes":

        st.session_state['has_taxonomic_info'] = 'Yes'      
        # Allow user to arrange columns
        taxonomic_order = st.multiselect(
            "Reorganize Columns by Hierarchical Levels and Exclude Unnecessary Columns",
            options=ordered_columns,
            default=ordered_columns,
            help=(
                "These are non-sample columns (e.g., taxonomy or metadata). "
                "You can choose which ones to keep and arrange their order for display. "
                "Unselected columns will be excluded from the table."
                )
        )

        st.session_state['taxonomic_order'] = taxonomic_order

        # Get all the sample columns other than the taxonomic columns
        valid_sample_columns = st.session_state.get("valid_sample_columns", set())
        sample_columns = [col for col in valid_sample_columns if col not in taxonomic_order]

        #### Rearranged Omics Quant Table:
        rearranged_table = omics_ft[taxonomic_order + sample_columns]
        
        # Drop rows where all values are None/NaN
        rearranged_table = rearranged_table.dropna(how='all')
        rearranged_table = rearranged_table.dropna(axis=1, how='all')
        
        st.session_state['rearranged_omics_table'] = rearranged_table

        st.markdown("#### Rearranged Omics Quant Table")
        st.write(
            f"The table is rearranged with the first **{len(taxonomic_order)}** columns based on the hierarchy selected by the user, "
            f"followed by **{len(sample_columns)}** sample columns. A total of **{rearranged_table.shape[0]}** unique features (rows) are included."  
            )
        st.dataframe(rearranged_table)

        #here last rows are full of None
        st.info(
            "Once you're satisfied with how the (proteomics/genomics) table looks, you can proceed to the next step: **Correlation Analysis**.\n\n"
            "This rearranged table will be used for correlation with the metabolomics data. "
            "The columns you selected and reordered based on the taxonomic or hierarchical structure will be used to **bin and group the data accordingly** during the correlation analysis.\n\n"
            "To continue, click on **'Correlation Analysis'** from the menu on the left."
        )
                  
    elif has_taxonomic_info == "No":

        st.session_state['has_taxonomic_info'] = 'No'   

        # Step 1: Get sample columns from metabolomics data
        valid_sample_columns = st.session_state['metabolome_ft'].columns.astype(str)

        # Step 2: Ensure omics_ft column names are strings
        omics_ft.columns = omics_ft.columns.astype(str)

        # Step 3: Identify numeric columns that match metabolome sample names
        numeric_columns = omics_ft.select_dtypes(include='number').columns
        matched_numeric_columns = [col for col in numeric_columns if col in valid_sample_columns]

        # Step 4: Identify dropped columns (either non-numeric or not in metabolomics table)
        dropped_columns = [col for col in omics_ft.columns if col not in matched_numeric_columns]

        # Step 5: Notify and filter
        if dropped_columns:
            st.warning(
                f"The following **{len(dropped_columns)} columns**  were removed because they either contain non-numeric data "
                f"or are not present in the metabolomics data:\n\n{', '.join(dropped_columns)}."
            )

        # Step 6: Filter table to matched numeric sample columns only
        omics_ft = omics_ft[matched_numeric_columns]
        omics_ft = omics_ft.dropna(how='all')
        omics_ft = omics_ft.dropna(axis=1, how='all')
        
        st.session_state['rearranged_omics_table'] = omics_ft
        st.write(f"Dimensions of the omics table: {omics_ft.shape[0]} rows × {omics_ft.shape[1]} columns")
        st.dataframe(omics_ft)

        st.info(
            "Once you're satisfied with how the (proteomics/genomics) table looks, you can proceed to the next step: **Correlation Analysis**.\n\n"
            "Since no hierarchical information was provided, a **one-to-one correlation** will be performed directly between the metabolomics data and the above omics data.\n\n"
            "Use the menu on the left to navigate to the next page."
        )
    else:
        # If the user hasn't selected any option yet, show a message
        st.warning("Please select an option to proceed.")

else:
    # If data is not available, display a message
    st.warning("Data for filtering is not available. Please load the data first.")

