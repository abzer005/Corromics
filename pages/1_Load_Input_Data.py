# Import necessary libraries
import streamlit as st
import pandas as pd
from src.common import *        # Importing common functionalities
from src.fileselection import * # Importing file selection functionalities

# Introduction Section
st.markdown("### Please select your method for data input below.")

# Input Selection Section
input_method = st.selectbox("Select Input Method", 
                            ["Use Example Dataset",
                             "Manual Input",
                             ],
                             index=0  # This sets "Use Example Dataset" as the default option
                             )

# Clearing the session state 
if 'last_input_method' not in st.session_state:
    st.session_state['last_input_method'] = None

# Example Dataset Section
elif input_method == "Use Example Dataset":
    # Check if input method has changed
    if st.session_state['last_input_method'] != input_method:
        # Clear the data
<<<<<<< HEAD
        for key in ['ft', 'md', 'omics_ft', 'omics_md']:
=======
        for key in ['ft', 'omics_ft', 'md', ]:
>>>>>>> b669908 (FDR added)
            st.session_state[key] = None

        # Update the last input method
        st.session_state['last_input_method'] = input_method
    
    load_example()  # Load data into session state

    for file_name, key in zip(["Metabolomics Feature Table", 
<<<<<<< HEAD
                               "Metabolomics MetaData", 
                               "Proteomics/Genomics Feature Table", 
                               "Proteomics/Genomics MetaData"],
                              ['ft', 'md', 'omics_ft', 'omics_md']):
        display_dataframe_with_toggle(key, file_name)
=======
                               "Proteomics/Genomics Feature Table",
                               "MetaData"],
                              ['ft', 'omics_ft', 'md']):
        display_dataframe_with_toggle(key, file_name)# Import necessary libraries

>>>>>>> b669908 (FDR added)

# Manual Input Section
elif input_method == "Manual Input":
    if st.session_state['last_input_method'] != input_method:
        # Clear the data
<<<<<<< HEAD
        for key in ['ft', 'md', 'omics_ft', 'omics_md']:
=======
        for key in ['ft', 'omics_ft', 'md']:
>>>>>>> b669908 (FDR added)
            st.session_state[key] = None
        # Update the last input method
        st.session_state['last_input_method'] = input_method

    st.info("💡 Upload tables in txt (tab separated), tsv, csv or xlsx (Excel) format.")

<<<<<<< HEAD
    # Create 2 columns for the ft, md file uploaders
    col1, col2 = st.columns(2)
=======
    # Create 3 columns for the ft, md file uploaders
    col1, col2, col3 = st.columns(3)
>>>>>>> b669908 (FDR added)
    with col1:
        ft_file = st.file_uploader("Upload Metabolomics Feature Table", 
                                   type=["csv", "xlsx", "txt", "tsv"],
                                   help = "This table is a key output of LC-MS/MS metabolomics studies. The table presents a list of mass spectral features along with their relative intensities (represented by its integrated peak area) observed across various samples.")
        if ft_file:
            st.session_state['ft'] = load_ft(ft_file).set_index("row ID")

    with col2:
<<<<<<< HEAD
        md_file = st.file_uploader("Upload Metabolomics Metadata", 
                                   type=["csv", "xlsx", "txt", "tsv"],
                                   help = "The metadata table is created by the user, providing additional context for the measured samples, such as sample type, species, and tissue type, etc.")
        if md_file:
            st.session_state['md'] = load_md(md_file).set_index("filename")
    
    # Create 2 columns for the nw, annotation file uploaders
    col3, col4 = st.columns(2)
    with col3:
=======
>>>>>>> b669908 (FDR added)
        omics_ft_file = st.file_uploader("Upload Proteomics/Genomics Feature Table", 
                                        type=["csv", "xlsx", "txt", "tsv"],
                                        help = ("This table represents the key output of proteomics or genomics studies, "
                                                "providing a list of proteins, genes, or other molecular entities along with their "
                                                "quantification (e.g., expression levels or abundances) across various samples. "
                                                "It is essential for integrating and correlating this dataset with metabolomics data."
                                                )
                                                )
        if omics_ft_file:
            st.session_state['omics_ft'] = load_omics_ft(omics_ft_file)
    
<<<<<<< HEAD
    with col4:
        omics_md_file = st.file_uploader("Upload Proteomics/Genomics MetaData", 
                                           type=["csv", "xlsx", "txt", "tsv"],
                                           help = "The metadata table is created by the user, providing additional context for the measured samples, such as sample type, species, and tissue type, etc."
                                           )
        
        if omics_md_file:
            st.session_state['omics_md'] = load_omics_md(omics_md_file)

    # Display headers and 'View all' buttons for each file
    for file_name, key in zip(["Metabolomics Feature Table", 
                               "Metabolomics MetaData", 
                               "Proteomics/Genomics Feature Table", 
                               "Proteomics/Genomics MetaData"],
                              ['ft', 'md', 'omics_ft', 'omics_md']):
=======
    with col3:
        md_file = st.file_uploader("Upload Metadata", 
                                   type=["csv", "xlsx", "txt", "tsv"],
                                   help = "The metadata table is created by the user, providing additional context for the measured samples, such as sample type, species, and tissue type, etc.")
        if md_file:
            st.session_state['md'] = load_md(md_file).set_index("filename")

    # Display headers and 'View all' buttons for each file
    for file_name, key in zip(["Metabolomics Feature Table", 
                               "Proteomics/Genomics Feature Table", 
                               "Metabolomics MetaData"],
                              ['ft', 'omics_ft', 'md']):
>>>>>>> b669908 (FDR added)
        display_dataframe_with_toggle(key, file_name)

else:
    # If data is not available, display a message
    st.warning("Data not loaded. Please load the data first.")


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
<<<<<<< HEAD

=======
    #md = md.dropna(subset=['ATTRIBUTE_Corromics_filename'])
    
>>>>>>> b669908 (FDR added)
    # If data is available, proceed with cleanup and checks
    cleaned_ft = clean_up_ft(ft)
    cleaned_md = clean_up_md(md)

    # Check if ft column names and md row names are the same
    cleaned_md, cleaned_ft = check_columns(cleaned_md, cleaned_ft)
    
<<<<<<< HEAD
    st.markdown("#### Metabolomics Metadata overview")
=======
    st.markdown("#### Metadata overview")
>>>>>>> b669908 (FDR added)
    df = inside_levels(cleaned_md)
    st.dataframe(df)

    st.session_state['ft_for_analysis'] = cleaned_ft
    st.session_state['md_for_analysis'] = cleaned_md

    ################################################################################
    st.markdown("#### Filter the Metabolomics Data")
    with st.container():
        c1, c2 = st.columns(2)
        # Allow the user to select any column for further filtering
        filter_column = c1.selectbox(
            "Select the metadata column for filtering",
            options=st.session_state['md_for_analysis'].columns,
            key="filter_column"
        )

        # Multi-select for categories in the selected column
        filter_group = c2.multiselect(
            "Select categories for filtering",
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

            # Display the final filtered data
            with st.expander(f"Filtered Group {final_ft.shape}"):
                st.dataframe(final_ft)
                st.dataframe(final_md)

            # Update session state with the final filtered tables
            st.session_state['metabolome_ft'] = final_ft
            st.session_state['metabolome_md'] = final_md
        
        else:
            st.warning("No groups selected.")     

#### Filter for the Genomics Data
if (
    'omics_ft' in st.session_state and 
<<<<<<< HEAD
    'omics_md' in st.session_state and 
    st.session_state['omics_ft'] is not None and 
    not st.session_state['omics_ft'].empty and 
    st.session_state['omics_md'] is not None and 
    not st.session_state['omics_md'].empty
):

    omics_ft = st.session_state['omics_ft'].copy()
    omics_md = st.session_state['omics_md'].copy()

    # Remove mzML and mzXML columns (case-insensitive)
    relevant_columns = [col for col in omics_ft.columns if not col.lower().endswith((".mzml", ".mzxml"))]
=======
    'metabolome_md' in st.session_state and 
    st.session_state['omics_ft'] is not None and 
    not st.session_state['omics_ft'].empty and 
    st.session_state['metabolome_md'] is not None and 
    not st.session_state['metabolome_md'].empty
):

    omics_ft = st.session_state['omics_ft'].copy()
    omics_md = st.session_state['metabolome_md'].copy()

    # Convert metadata column values to a set for fast lookup
    metadata_filenames = set(omics_md['ATTRIBUTE_Corromics_filename'].dropna().astype(str).str.strip())

    # Filter columns: Remove those ending with .mzML/.mzXML or present in metadata_filenames
    relevant_columns = [
        col for col in omics_ft.columns 
        if not col.lower().endswith((".mzml", ".mzxml")) and col not in metadata_filenames
    ]

    # Select only the relevant columns
>>>>>>> b669908 (FDR added)
    ordered_columns = order_taxonomic_columns(relevant_columns)

    st.markdown("#### Filter the Other Omics Data")
    # Ask if the table contains taxonomic information
    has_taxonomic_info = st.radio("Does the quant table contain taxonomic information?", ["Yes", "No"])

    if has_taxonomic_info == "Yes":
        # Allow user to arrange columns
        taxonomic_order = st.multiselect(
            "Reorganize Columns by Taxonomic Levels and Exclude Unnecessary Columns",
            options=ordered_columns,
            default=ordered_columns,
            help="The columns other than the samples are listed here"
        )

<<<<<<< HEAD
        # Rearrange and display the table
        if taxonomic_order:
            st.session_state['taxonomic_order'] = taxonomic_order

            st.write("#### Rearranged Omics Quant Table:")
            rearranged_table = omics_ft[taxonomic_order + [col for col in omics_ft.columns if col not in relevant_columns]]
            st.dataframe(rearranged_table)

            # Allow the user to select the taxonomic level from the available options
            selected_level = st.selectbox("Select a taxonomic level to bin the data:", taxonomic_order)

            # Perform binning based on the selected level
            binned_by_level = bin_by_taxonomic_level(rearranged_table, selected_level)
            binned_level_filtered = binned_by_level[binned_by_level['Overall_sum'] > 0]
            st.session_state['binned_omics_table'] = binned_level_filtered

            with st.expander(f"Binned Data at Level: {selected_level} , Original Dimension: {binned_by_level.shape}"):
                st.dataframe(binned_by_level)

            # Display the binned table
            with st.expander(f"Binned Data at Level: {selected_level} , Filtered Dimension: {binned_level_filtered.shape}"):
                st.dataframe(binned_level_filtered)

            
    else:
        st.session_state['binned_omics_table'] = omics_ft


    # If data is available, proceed with cleanup and checks
    cleaned_omics_ft = clean_up_omics_ft(st.session_state['binned_omics_table'])
    #st.dataframe(cleaned_omics_ft)
    cleaned_omics_md = clean_up_omics_md(omics_md)

    # Check if ft column names and md row names are the same
    cleaned_omics_md, cleaned_omics_ft = check_columns(cleaned_omics_md, cleaned_omics_ft)
    st.markdown("#### Metagenomics Metadata overview")
    df = inside_levels(cleaned_omics_md)
    st.dataframe(df)

    st.session_state['omics_ft_analysis'] = cleaned_omics_ft
    st.session_state['omics_md_analysis'] = cleaned_omics_md

    ###########################################################
    st.markdown("#### Filter the Proteomics/Metagenomics Data")
    with st.container():
        c1, c2 = st.columns(2)
        # Allow the user to select any column for further filtering
        omics_filter_column = c1.selectbox(
            "Select the metadata column for filtering",
            options=st.session_state['omics_md_analysis'].columns,
            key="omics_filter_column"
        )

        # Multi-select for categories in the selected column
        omics_filter_category = c2.multiselect(
            "Select categories for filtering",
            options=sorted(st.session_state['omics_md_analysis'][omics_filter_column].dropna().unique()),
            key="omics_filter_category"
        )

        # Apply the filter if categories are selected
        if omics_filter_category:
 
            omics_filter_category = list(map(str, omics_filter_category))  # Convert the group to strings if needed for matching
            omics_filter_indices = cleaned_omics_md[cleaned_omics_md[omics_filter_column].astype(str).isin(omics_filter_category)].index
                        
            # Update the feature table and metadata based on additional filtering
            final_omics_ft = cleaned_omics_ft.loc[:, omics_filter_indices]
            final_omics_md = cleaned_omics_md.loc[omics_filter_indices]

            # Display the final filtered data
            with st.expander(f"Filtered Group {final_omics_ft.shape}"):
                st.dataframe(final_omics_ft)
                st.dataframe(final_omics_md)

            # Update session state with the final filtered tables
            st.session_state['genomics_ft'] = final_omics_ft
            st.session_state['genomics_md'] = final_omics_md
        
        else:
            st.warning("No groups selected.")
            st.session_state['genomics_ft'] = cleaned_omics_ft
            st.session_state['genomics_md'] = cleaned_omics_md





    ################################################################################  

else:
    # If data is not available, display a message
    st.warning("Data not loaded. Please load the data first.")
=======
        st.session_state['taxonomic_order'] = taxonomic_order
        st.write("#### Rearranged Omics Quant Table:")
        rearranged_table = omics_ft[taxonomic_order + [col for col in omics_ft.columns if col not in relevant_columns]]
        st.dataframe(rearranged_table)
        st.session_state['rearranged_omics_table'] = rearranged_table
                    
    else:
        st.session_state['binned_omics_table'] = omics_ft
    
else:
    # If data is not available, display a message
    st.warning("Data not loaded. Please load the data first.")
>>>>>>> b669908 (FDR added)
