import streamlit as st
from .common import *  # Importing common functionalities from the 'common' module
import pandas as pd
import re
import numpy as np

patterns = [
    ["m/z", "mz", "mass over charge"],
    ["rt", "retention time", "retention-time", "retention_time"],
]

allowed_formats = "Allowed formats: csv (comma separated), tsv (tab separated), txt (tab separated), xlsx (Excel file)."

def string_overlap(string, options):
    """
    Check if any of the given options are present in the string.
    Exclude any string containing "mzml".

    Parameters:
    string (str): The string to be checked.
    options (list): A list of substrings to search for in the string.

    Returns:
    bool: True if any option is found in the string and "mzml" is not present, False otherwise.
    """
    if "mzml" in string or "mzxml" in string:
        return False
    
    return any(option in string for option in options)

def load_example():
    """
    Load example datasets into Streamlit's session state.
    """
    # Reset session state data
<<<<<<< HEAD
    for key in ['ft', 'md', 'omics_ft', 'omics_md']:
        st.session_state[key] = None
        
    st.session_state['ft'] = open_df("example-data/Normalised_Quant_table.csv").set_index("feature_ID")
    st.session_state['md'] = open_df("example-data/metadata_metabolomics.csv").set_index("filename")
    st.session_state['omics_ft'] = open_df("example-data/asv_16s_table_with_taxonomic_levels.csv").set_index("feature_ID")
    st.session_state['omics_md'] = open_df("example-data/asv_16s_metadata.csv").set_index("filename")

=======
    for key in ['ft', 'md', 'omics_ft']:
        st.session_state[key] = None
        
    st.session_state['ft'] = open_df("example-data/Normalised_Quant_table.csv").set_index("feature_ID")
    st.session_state['md'] = open_df("example-data/metadata_example.csv").set_index("filename")
    st.session_state['omics_ft'] = open_df("example-data/asv_16s_table_with_taxonomic_levels.csv").set_index("feature_ID")
    
>>>>>>> b669908 (FDR added)
def load_ft(ft_file):
    """
    Load and process the feature table.

    Parameters:
    ft_file (file): The feature table file.

    Returns:
    DataFrame: Processed feature table.
    """
    ft = open_df(ft_file)
    ft = ft.dropna(axis=1)  # Drop columns with missing values
    return ft

def load_md(md_file):
    """
    Load and process metadata. Set 'filename' as the index if present.

    Parameters:
    md_file (file): The metadata file.

    Returns:
    DataFrame: Processed metadata.
    """
    md = open_df(md_file)
    return md

<<<<<<< HEAD
def load_nw(omics_ft_file):
=======
def load_omics_ft(omics_ft_file):
>>>>>>> b669908 (FDR added)
    """
    Load and process the quantification table from proteomics/genomics study. 
    """
    omics_ft = open_df(omics_ft_file)
    return omics_ft

<<<<<<< HEAD
def load_annotation(omics_md_file):
    """
    Load and process the metadata from proteomics/genomics study.
    """
    omics_md = open_df(omics_md_file)
    return omics_md

=======
>>>>>>> b669908 (FDR added)
def display_dataframe_with_toggle(df_key, display_name):
    if df_key in st.session_state and isinstance(st.session_state[df_key], pd.DataFrame):
        st.write(f"### {display_name}")

        col1, col2 = st.columns([0.8, 0.2])

        # Show dimensions
        num_rows, num_cols = st.session_state[df_key].shape
        col1.write(f"Dimension: {num_rows} rows × {num_cols} columns")

        view_all = col2.checkbox("View all", key=f"{df_key}_toggle")

        if view_all:
            st.dataframe(st.session_state[df_key])  # Show full dataframe
        else:
            st.dataframe(st.session_state[df_key].head())  # Show header


##### Cleanup functions

@st.cache_data
def clean_up_md(md):
    md = md.copy()
<<<<<<< HEAD
    md = md.dropna(how="all")
    md.index = [name.strip() for name in md.index]
=======

     # Remove rows that have NA values in ATTRIBUTE_Corromics_filenam
    md = md.dropna(subset=['ATTRIBUTE_Corromics_filename'])

     # Remove rows where all values are NA
    md = md.dropna(how="all")
    md.index = [name.strip() for name in md.index]

>>>>>>> b669908 (FDR added)
    # for each col in md
    # 1) removing the spaces (if any)
    # 2) replace the spaces (in the middle) to underscore
    # 3) converting them all to UPPERCASE
<<<<<<< HEAD
=======

>>>>>>> b669908 (FDR added)
    for col in md.columns:
        if md[col].dtype == str:
            md[col] = [item.strip().replace(" ", "_").upper() for item in md[col]]

    md.index = [
        re.sub(r"\.mzxml|\.mzml", "", i, flags=re.IGNORECASE).replace(" Peak area", "")
        for i in md.index
    ]
<<<<<<< HEAD
    #md.index = [i.replace(".mzXML", "").replace(".mzML", "").replace(" Peak area", "") for i in md.index]
=======
    
>>>>>>> b669908 (FDR added)
    return md


@st.cache_data
def clean_up_ft(ft):
    ft = ft.copy()
    ft = ft.dropna(how="all")

    # drop all columns that are not mzML or mzXML file names
    ft.drop(
        columns=[col for col in ft.columns if not re.search(r"\.mzml|\.mzxml", 
                                                            col, 
                                                            flags=re.IGNORECASE)],
        inplace=True,
    )

    # remove " Peak area" from column names, contained after mzmine pre-processing
    ft.rename(
        columns={
            col: re.sub(r"\.mzxml|\.mzml", "", col, flags=re.IGNORECASE)
                 .replace(" Peak area", "")
                 .strip()
            for col in ft.columns
        },
        inplace=True,
    )

    return ft

@st.cache_data
def clean_up_omics_md(md):
    md = md.copy()
    md = md.dropna(how="all")
    md.index = [name.strip() for name in md.index]
    # for each col in md
    # 1) removing the spaces (if any)
    # 2) replace the spaces (in the middle) to underscore
    # 3) converting them all to UPPERCASE
    for col in md.columns:
        if md[col].dtype == str:
            md[col] = [item.strip().replace(" ", "_").upper() for item in md[col]]

    md.index = [
        re.sub(r"\.mzxml|\.mzml", "", i, flags=re.IGNORECASE).replace(" Peak area", "")
        for i in md.index
    ]

<<<<<<< HEAD
    # Keep taxonomic columns if they exist in session state
    #taxonomic_columns = st.session_state.get("taxonomic_order", [])
    #if taxonomic_columns:
        # Ensure we keep only relevant taxonomic columns
    #    md = md[taxonomic_columns + [col for col in md.columns if col not in taxonomic_columns]]

=======
>>>>>>> b669908 (FDR added)
    return md


@st.cache_data
def clean_up_omics_ft(ft):
    """
    Cleans up an omics quantification table by:
    - Keeping mzML or mzXML file names and taxonomic columns if specified in session state.
    - Removing unwanted substrings from column names (case-insensitive).
    """
    ft = ft.copy()  # Preserve the original file
    ft = ft.dropna(how="all")

    # Check if `taxonomic_order` exists in session state
    taxonomic_columns = st.session_state.get("taxonomic_order", [])

    # Drop all columns that are not mzML/mzXML file names or taxonomic columns
    ft.drop(
        columns=[
            col for col in ft.columns 
            if not (
                re.search(r"\.mzml|\.mzxml", col, flags=re.IGNORECASE) or col in taxonomic_columns
            )
        ],
        inplace=True,
    )

    # Remove " Peak area", ".mzXML", ".mzML" from column names (case-insensitive)
    ft.rename(
        columns={
            col: re.sub(r"\.mzxml|\.mzml", "", col, flags=re.IGNORECASE)
                 .replace(" Peak area", "")
                 .strip()
            for col in ft.columns
        },
        inplace=True,
    )

    return ft

@st.cache_data
def check_columns(md, ft):
    if sorted(ft.columns) != sorted(md.index):
<<<<<<< HEAD
        st.warning("Not all files are present in both meta data & feature table.")
=======
        st.warning("Not all files are present in both metadata & feature table.")
>>>>>>> b669908 (FDR added)
        
        # Find and remove columns in 'ft' that are not in the index of 'md'
        ft_cols_not_in_md = [col for col in ft.columns if col not in md.index]
        if ft_cols_not_in_md:
            st.warning(
<<<<<<< HEAD
                f"These {len(ft_cols_not_in_md)} columns of feature table are not present in metadata table and will be removed:\n{', '.join(ft_cols_not_in_md)}"
            )
=======
                f"These {len(ft_cols_not_in_md)} columns of feature table are not present in metadata column **'ATTRIBUTE_Corromics_filename'** and will be removed:\n\n"
                + f"\n\n{', '.join(ft_cols_not_in_md)}"
                )

>>>>>>> b669908 (FDR added)
            ft = ft.drop(columns=ft_cols_not_in_md)
        
        # Find and remove rows in 'md' that are not in the columns of 'ft'
        md_rows_not_in_ft = [row for row in md.index if row not in ft.columns]
        if md_rows_not_in_ft:
            st.warning(
                f"These {len(md_rows_not_in_ft)} rows of metadata table are not present in feature table and will be removed:\n{', '.join(md_rows_not_in_ft)}"
            )
            md = md.drop(md_rows_not_in_ft)
    return md, ft

<<<<<<< HEAD

# @st.cache_data
# def inside_levels(df):
#     df = pd.DataFrame(
#         {
#             "ATTRIBUTES": df.columns,
#             "LEVELS": [sorted(set(df[col].dropna().astype(str).to_list())) for col in df],
#             "COUNTS": [df[col].value_counts().to_list() for col in df],
#         }
#     )
#     return df

=======
>>>>>>> b669908 (FDR added)
@st.cache_data
def inside_levels(df):

    result = []

    for col in df.columns:
        # Convert all values to string (including NaN as 'NaN') for uniformity
        values = df[col].astype(str)
        
        # Create a dictionary of levels and counts (including NaN)
        level_count_dict = values.value_counts(dropna=False).to_dict()

        # Sort levels alphabetically/numerically, keeping 'NaN' visible
        sorted_levels = sorted(level_count_dict.keys(), key=lambda x: (x != "nan", x))

        # Extract counts corresponding to the sorted levels
        sorted_counts = [level_count_dict[level] for level in sorted_levels]

        # Append to result
        result.append({
            "ATTRIBUTES": col,
            "LEVELS": sorted_levels,
            "COUNTS": sorted_counts,
        })

    # Convert to DataFrame
    df = pd.DataFrame(result)
    return df


@st.cache_data
def order_taxonomic_columns(relevant_columns):
    """
    Orders the relevant columns based on a predefined taxonomic hierarchy.

    Parameters:
    relevant_columns (list): List of column names to be ordered.

    Returns:
    list: Ordered list of column names.
    """
    # Predefined taxonomic hierarchy
    hierarchy = [
        "domain", "domains",
        "kingdom", "kingdoms",
        "phylum", "phyla",
        "class", "classes",
        "order", "orders",
        "family", "families",
        "genus", "genera",
        "species"
    ]

    # Function to determine the hierarchy position of a column
    def get_order(column):
        column_lower = column.lower()  # Convert to lowercase for comparison
        for i, term in enumerate(hierarchy):
            if term in column_lower:
                return i
        return len(hierarchy)  # Place unmatched columns at the end

    # Sort relevant columns based on the hierarchy
    ordered_columns = sorted(relevant_columns, key=get_order)
    return ordered_columns

<<<<<<< HEAD


=======
>>>>>>> b669908 (FDR added)
@st.cache_data
def bin_by_taxonomic_level(df, taxonomic_level):
    """
    Bins the DataFrame based on a selected taxonomic level and sums the remaining columns.

    Parameters:
    df (pd.DataFrame): The input DataFrame with taxonomic columns and numeric columns to sum.
    taxonomic_level (str): The taxonomic level to bin by (e.g., "Domain", "Phylum", "Genus").

    Returns:
    pd.DataFrame: The binned DataFrame.
    """
    taxonomic_columns = st.session_state.get("taxonomic_order", [])

    # Ensure the taxonomic level is valid
    if taxonomic_level not in taxonomic_columns:
        raise ValueError(f"Invalid taxonomic level: {taxonomic_level}. Must be one of {taxonomic_columns}.")

    # Find the index of the taxonomic level
    level_index = taxonomic_columns.index(taxonomic_level)

    # Filter rows where any of the taxonomic columns up to the selected level are NaN
    taxonomic_columns_to_check = taxonomic_columns[:level_index + 1]
    df = df.dropna(subset=taxonomic_columns_to_check)

    # Create a combined column for grouping
    group_column = "_".join(taxonomic_columns[:level_index + 1])
    df[group_column] = df[taxonomic_columns[:level_index + 1]].astype(str).agg("_".join, axis=1)

    # Group by the combined column and sum the numeric values
    numeric_columns = df.select_dtypes(include="number").columns
    binned_df = df.groupby(group_column)[numeric_columns].sum().reset_index()

    # Add 'Overall_sum' column
    binned_df["Overall_sum"] = binned_df[numeric_columns].sum(axis=1)
<<<<<<< HEAD
=======
   
>>>>>>> b669908 (FDR added)
    binned_df.reset_index(inplace=True)
    binned_df.set_index(group_column, inplace=True)

    return binned_df


<<<<<<< HEAD
=======
# Function to check for empty rows, rows with a singular same value, or rows full of zeros

def check_rows(df, exclude_cols=['index', 'overall_sum']):
    # Select only numeric columns, excluding specified ones
    numeric_cols = df.select_dtypes(include=['number']).columns.difference(exclude_cols)

    # 1. Empty rows (all NaN or blank in numeric columns)
    empty_rows = df[df[numeric_cols].isnull().all(axis=1) | 
                    (df[numeric_cols].applymap(lambda x: str(x).strip() == "").all(axis=1))]

    # 2. Rows with a single unique value (ignoring NaN)
    singular_value_rows = df[df[numeric_cols].nunique(axis=1) == 1]

    # 3. Rows full of zeros (after filling NaN with 0)
    zero_rows = df[df[numeric_cols].fillna(0).sum(axis=1) == 0]

    return empty_rows, singular_value_rows, zero_rows
>>>>>>> b669908 (FDR added)

