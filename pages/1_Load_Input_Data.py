# Import necessary libraries
import streamlit as st
import pandas as pd
from src.common import *        # Importing common functionalities
from src.fileselection import * # Importing file selection functionalities

page_setup()
initialize_app()

# Introduction Section
st.markdown("## Please select your method for data input below.")

# -------------------------------
# Input Selection Section
# -------------------------------
input_method = st.selectbox(
    "Select Input Method", 
    [
        "Use Example Dataset",
        "Manual Input (Custom Data)",
        "FBMN-Stats Output (Metabolomics Table)",
        "MZmine Export (Metabolomics Table)",
    ],
    index=1,
    key="input_method",
)

# --- Initialize session_state keys if not present ---
for key in ['md', 'ft', 'omics_ft', 'last_input_method']:
    if key not in st.session_state:
        st.session_state[key] = None

# --- Detect change in input method ---
if st.session_state.get('last_input_method') != input_method:
    # Clear all previous dataframes
    for key in ['md', 'ft', 'omics_ft']:
        st.session_state[key] = None

    # Update the tracker
    st.session_state['last_input_method'] = input_method
    st.rerun() # Rerun to reset the UI

# ============================================================
# 1. Example Dataset (special: loads all three)
# ============================================================
if input_method == "Use Example Dataset":
    if (
        st.session_state['md'] is None
        and st.session_state['ft'] is None
        and st.session_state['omics_ft'] is None
    ):
        load_example()  # should fill md, ft, omics_ft

    show_input_tables_in_tabs()
    #st.stop()  # don't show uploaders below in this mode

# ==========================
# All other input methods
else:
    # ---- 1. Metadata + Other Omics (global for all non-example modes) ----
    st.markdown("#### Upload Metadata and Other Omics Tables")

    col1, col2 = st.columns(2)

    with col1:
        md_file = st.file_uploader(
            "Upload Metadata",
            type=["csv", "xlsx", "txt", "tsv"],
            key="md_upload",
            help=(
                "Metadata table should contain `filename` and "
                "`ATTRIBUTE_Corromics_filename` columns."
            ),
        )
        if md_file is not None:
            st.session_state['md'] = load_md(md_file)

    with col2:
        omics_ft_file = st.file_uploader(
            "Upload Other Omics Feature Table",
            type=["csv", "xlsx", "txt", "tsv"],
            key="omics_upload",
            help=(
                "Quantification table from ASV sequencing / proteomics / other omics. "
                "The **first column** will be treated as the feature identifier, "
                "renamed to `feature_ID`, and used as the index."
            ),
        )
        if omics_ft_file is not None:
            st.session_state['omics_ft'] = load_omics_ft(omics_ft_file)

    # ---- 2. Metabolomics feature table (depends on input_method) ----
    st.markdown("#### Metabolomics Feature Table")

    # Manual Input
    if input_method == "Manual Input (Custom Data)":
        ft_file = st.file_uploader(
            "Upload Metabolomics Feature Table",
            type=["csv", "xlsx", "txt", "tsv"],
            key="ft_upload",
            help=(
                "The app will assume features as rows "
                "and samples as columns after loading/processing."
            ),
        )
        if ft_file is not None:
            st.session_state['ft'] = load_ft(ft_file)

    # FBMN-Stats Output
    elif input_method == "FBMN-Stats Output (Metabolomics Table)":
        ft_file = st.file_uploader(
            "Upload FBMN-Stats Feature Table",
            type=["csv", "xlsx", "txt", "tsv"],
            key="ft_stats_upload",
            help=(
                "Feature table exported from the FBMN-Stats app. "
                "Expected format: samples as rows, features as columns; "
                "feature names like `ID_mz_RT` or `ID_mz_RT&name`."
            ),
        )
        if ft_file is not None:
            raw_ft = open_stats_ft(ft_file)
            if not raw_ft.empty and st.session_state.get('md') is not None:
                st.session_state['ft'] = load_ft_from_statsapp(
                    raw_ft,
                    st.session_state['md'],
                )
            elif not raw_ft.empty and st.session_state.get('md') is None:
                st.warning("Please upload metadata first so sample names can be matched.")

    # MZmine Export
    elif input_method == "MZmine Export (Metabolomics Table)":
        ft_file = st.file_uploader(
            "Upload MZmine Feature Table",
            type=["csv", "xlsx", "txt", "tsv"],
            key="ft_mzmine_upload",
            help=(
                "Feature table exported from MZmine. Expected format: features as rows, "
                "samples as columns, with 'row ID', 'row m/z', and 'row retention time' "
                "columns to construct a unique `feature_ID`."
            ),
        )
        if ft_file is not None:
            raw_ft = open_df(ft_file)
            if not raw_ft.empty and st.session_state.get('md') is not None:
                st.session_state['ft'] = load_ft_from_mzmine(
                    raw_ft,
                    st.session_state['md'],
                )
            elif not raw_ft.empty and st.session_state.get('md') is None:
                st.warning("Please upload metadata first so sample names can be matched.")

    show_input_tables_in_tabs()

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
#### Filter for the other omics Data
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

    taxonomy_columns = [
        col for col in omics_ft.columns if col not in valid_sample_columns
        ]
    
    taxonomy_columns = sorted(taxonomy_columns)

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

