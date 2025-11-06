import streamlit as st
from .common import *  # Importing common functionalities from the 'common' module
import pandas as pd
import numpy as np
import re

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
    for key in ['md','ft', 'omics_ft']:
        st.session_state[key] = None
        
    st.session_state['md'] = open_df("example-data/metadata_example.csv").set_index("filename")
    st.session_state['ft'] = open_df("example-data/Normalised_Quant_table.csv").set_index("feature_ID")
    st.session_state['omics_ft'] = open_df("example-data/asv_16s_table_with_taxonomic_levels.csv").set_index("feature_ID")

@st.cache_data
def load_md(md_file):
    """
    Load and process metadata. Set 'filename' as the index if present.
    Parameters: md_file (file): The metadata file.
    Returns: Processed metadata.
    """
    md = open_df(md_file).set_index("filename")
    return md
    

@st.cache_data
def load_ft(ft_file):
    """
    Load and process the feature table.
    Parameters: ft_file (file): The feature table file.
    Returns: Processed feature table.
    """
    ft = open_df(ft_file).set_index("feature_ID")
    ft = ft.dropna(axis=1)  # Drop columns with missing values
    return ft


@st.cache_data
def load_ft_from_statsapp(
    ft_df: pd.DataFrame,
    md_df: pd.DataFrame,
):
    """
    Process FBMN-Stats feature table:
    - Remove metadata columns (keep only feature columns starting with a digit)
    - Detect sample identifiers from index
    - Map sample IDs to full filenames (e.g. 'spiked_1.mzML') using metadata
    - Transpose so rows = features, cols = samples
    - Set index name to 'feature_ID'
    """

 # -------------------------
     # --- 1. Use first column as sample IDs and drop it from data ---
    sample_col = ft_df.columns[0]
    ft_df.index = ft_df.iloc[:, 0].astype(str)
    ft_df.drop(columns=[sample_col], inplace=True)
    
    # --- 2. Get filenames from metadata (index = filename) ---
    if md_df.index.name is None:
        st.error("Metadata index is not set. Please set it to the 'filename' column.")
        return None
    
    md_filenames = pd.Series(md_df.index.astype(str), index=md_df.index)
    md_basenames = md_filenames.str.rsplit(".", n=1).str[0]
    basename_to_full = dict(zip(md_basenames, md_filenames.astype(str)))

    # --- 3. Keep only feature columns (start with a digit) ---
    def is_feature_col(col):
        col_str = str(col).strip()
        return len(col_str) > 0 and col_str[0].isdigit()

    feature_cols = [c for c in ft_df.columns if is_feature_col(c)]
    if not feature_cols:
        st.error(
            "No feature columns detected in FBMN-Stats table "
            "(no columns starting with a digit)."
        )
        return None

    ft = ft_df[feature_cols].copy()

    # # --- 4. Map sample IDs (index) to full filenames ---
    new_index = []
    missing = []

    for idx in ft.index.astype(str):
        if idx in basename_to_full:
            new_index.append(basename_to_full[idx])  # e.g. file1 -> file1.mzML
        else:
            new_index.append(idx)
            missing.append(idx)

    ft.index = new_index

    if missing:
        st.warning(
            f"{len(missing)} sample(s) could not be matched to metadata filenames. "
            f"Example unmatched: {list(missing)[:5]}"
        )

    ft_t = ft.T
    ft_t.index.name = "feature_ID"

    return ft_t


@st.cache_data
def load_ft_from_mzmine(ft_df: pd.DataFrame,
                        md_df: pd.DataFrame):
    """
    Process MZmine feature table:
    - Features are rows, samples are columns
    - Build feature_ID from 'row ID', 'row m/z', 'row retention time'
    - Set feature_ID as index
    - Match sample columns to metadata filenames (md_df.index)
      * If metadata filenames contain 'Peak area' -> match directly
      * Otherwise, strip ' Peak area' suffix from MZmine column names
        when matching and rename columns to metadata filenames.
    """

    # --- 1. Check required columns ---
    required_cols = ["row ID", "row m/z", "row retention time"]
    missing_req = [c for c in required_cols if c not in ft_df.columns]
    if missing_req:
        st.error(f"MZmine feature table is missing required columns: {missing_req}")
        return None

    # Metadata must use filenames as index
    if md_df.index.name is None:
        st.error(
            "Metadata index is not set.\n"
            "Please ensure the metadata table uses the filename column as its index."
        )
        return None

    ft = ft_df.copy()

    # --- 2. Build feature_ID and set as index ---
    ft["feature_ID"] = (
        ft["row ID"].astype(str)
        + "_"
        + ft["row m/z"].astype(str)
        + "_"
        + ft["row retention time"].astype(str)
    )

    ft.set_index("feature_ID", inplace=True)
    ft.index.name = "feature_ID"

    # Drop the original ID/mz/RT columns
    ft.drop(columns=required_cols, inplace=True)

    # --- 3. Prepare metadata filenames and Peak area logic ---
    md_filenames = set(md_df.index.astype(str))
    md_has_peak_area = any("peak area" in s.lower() for s in md_filenames)

    sample_cols = []
    rename_map = {}  # old_col -> new_col (metadata filename)

    # --- 4. Select and align sample columns ---
    if md_has_peak_area:
        # Case 1: metadata filenames already include 'Peak area' – match directly
        for c in ft.columns:
            c_str = str(c)
            if c_str in md_filenames:
                sample_cols.append(c)
        info_msg = "Matching MZmine columns directly to metadata filenames."
    else:
        # Case 2: metadata filenames do NOT include 'Peak area'
        # Strip ' Peak area' suffix from MZmine column names when matching
        for c in ft.columns:
            c_str = str(c)
            lower = c_str.lower()
            base = c_str

            # remove ' peak area' suffix at the end, case-insensitive
            suffix = " peak area"
            if lower.endswith(suffix):
                base = c_str[:len(c_str) - len(suffix)].rstrip()

            if base in md_filenames:
                sample_cols.append(c)
                rename_map[c] = base  # rename to pure filename

        info_msg = (
            "Matching MZmine columns to metadata filenames after stripping "
            "' Peak area' suffix where present."
        )

    if sample_cols:
        st.info(info_msg)

    dropped_cols = [c for c in ft.columns if c not in sample_cols]
    if dropped_cols:
        st.info(
            f"Dropping {len(dropped_cols)} non-sample columns from MZmine table "
            "(not present in metadata filenames)."
        )

    # Keep only matching columns
    ft = ft[sample_cols].copy()

    # Apply renaming (for the non-peak-area metadata case)
    if rename_map:
        ft.rename(columns=rename_map, inplace=True)

    return ft



def load_omics_ft(omics_ft_file):
    """
    Load and process the quantification table from proteomics/genomics study. 
    """
 # 1) Read file into DataFrame
    df = open_df(omics_ft_file)
    if df is None or df.empty:
        return df  # or raise / log if you prefer

    # 2) Standardize first column as feature_ID
    first_col = df.columns[0]
    if first_col != "feature_ID":
        df.rename(columns={first_col: "feature_ID"}, inplace=True)

    df.set_index("feature_ID", inplace=True)
    df.index.name = "feature_ID"

    return df

def standardize_other_omics_ft(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize the 'other omics' feature table:
    - Treat the first column as feature IDs
    - Rename it to 'feature_ID'
    - Set it as the index
    """
    df = df.copy()
    first_col = df.columns[0]
    if first_col != "feature_ID":
        df.rename(columns={first_col: "feature_ID"}, inplace=True)
    df.set_index("feature_ID", inplace=True)
    df.index.name = "feature_ID"
    return df

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


def show_input_tables_in_tabs():
    """
    Display Metadata, Metabolomics FT, and Other Omics FT
    in three tabs with shape info and full table view.
    """
    tab_defs = [
        ("md", "Metadata"),
        ("ft", "Metabolomics Feature Table"),
        ("omics_ft", "Other Omics Feature Table"),
    ]

    tabs = st.tabs([label for _, label in tab_defs])

    for tab, (key, label) in zip(tabs, tab_defs):
        with tab:
            df = st.session_state.get(key)

            if isinstance(df, pd.DataFrame) and not df.empty:
                # Show shape info
                num_rows, num_cols = df.shape
                st.caption(f"{num_rows} rows × {num_cols} columns")

                # Show full dataframe
                st.dataframe(df, use_container_width=True)
            else:
                st.info(f"{label} not loaded yet.")

##### Cleanup functions

@st.cache_data
def clean_up_md(md):
    md = md.copy()

     # Remove rows that have NA values in ATTRIBUTE_Corromics_filenam
    md = md.dropna(subset=['ATTRIBUTE_Corromics_filename'])

     # Remove rows where all values are NA
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
        st.warning("Not all files are present in both metadata & feature table.")
        
        # Find and remove columns in 'ft' that are not in the index of 'md'
        ft_cols_not_in_md = [col for col in ft.columns if col not in md.index]
        if ft_cols_not_in_md:
            st.warning(
                f"These {len(ft_cols_not_in_md)} columns of feature table are not present in metadata column **'ATTRIBUTE_Corromics_filename'** and will be removed from the feature table:\n\n"
                + f"\n\n{', '.join(ft_cols_not_in_md)}"
                )

            ft = ft.drop(columns=ft_cols_not_in_md)
        
        # Find and remove rows in 'md' that are not in the columns of 'ft'
        md_rows_not_in_ft = [row for row in md.index if row not in ft.columns]
        if md_rows_not_in_ft:
            st.warning(
                f"These {len(md_rows_not_in_ft)} rows of metadata table are not present in feature table and will be removed from the metadata table:\n{', '.join(md_rows_not_in_ft)}"
            )
            md = md.drop(md_rows_not_in_ft)
    return md, ft

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
   
    binned_df.reset_index(inplace=True)
    binned_df.set_index(group_column, inplace=True)

    return binned_df

@st.cache_data
def bin_by_flat_id(df):
    # Step 1: Drop fully empty rows and columns
    df_clean = df.dropna(how='all').dropna(axis=1, how='all')

    # Step 2: Keep only numeric columns
    numeric_df = df_clean.select_dtypes(include="number").copy()

    # Step 3: Add overall sum across numeric values
    numeric_df["Overall_sum"] = numeric_df.sum(axis=1)

    # Step 4: Return DataFrame with original index
    return numeric_df


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

def filter_by_overall_sum(df, exclude_cols=None, label_prefix=""):
    """
    Adds UI widgets to filter a DataFrame by 'Overall_sum' column.
    
    Parameters:
        df (DataFrame): Must contain an 'Overall_sum' column.
        exclude_cols (list): Columns to exclude from later transformations.
        label_prefix (str): Optional label prefix to distinguish keys if re-used in multiple places.

    Returns:
        filtered_df (DataFrame): DataFrame filtered by threshold sliders.
        exclude_cols (list): The columns to exclude in downstream processing.
    """
    if exclude_cols is None:
        exclude_cols = ['index', 'Overall_sum']

    max_overall_sum = df['Overall_sum'].max()
    filter_col1, filter_col2 = st.columns(2)

    with filter_col1:
        st.session_state[f'{label_prefix}filter_threshold_below'] = st.number_input(
            "**Filter out data with reads/intensities below the user input. Press Enter to apply.**", 
            min_value=0.0, 
            value=0.0, 
            step=1.0,
            help="**All variables with a total read/intensity below this threshold will be removed from analysis.**",
            key=f"{label_prefix}filter_below"
        )

    with filter_col2:
        st.session_state[f'{label_prefix}filter_threshold_above'] = st.number_input(
            "**Filter out data with reads/intensities above this value. Press Enter to apply.**", 
            max_value=float(max_overall_sum),
            value=float(max_overall_sum),
            step=1.0,
            help="**All variables with a total read/intensity ABOVE this threshold will be removed from analysis.**",
            key=f"{label_prefix}filter_above"
        )

    filtered_df = df[
        (df['Overall_sum'] > st.session_state[f'{label_prefix}filter_threshold_below']) &
        (df['Overall_sum'] <= st.session_state[f'{label_prefix}filter_threshold_above'])
    ]

    return filtered_df, exclude_cols



