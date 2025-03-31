import streamlit as st
from .common import *  # Importing common functionalities from the 'common' module
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
import psutil
import os

def is_running_locally():
    """Return True if app is running locally (localhost), False if on cloud."""
    return os.environ.get("REMOTE_ADDR", "") in ("127.0.0.1", "localhost", "0.0.0.0", "")


#import cupy as cp

def combine_dataframes(df1, df2):

    # Align columns by names
    common_columns = df1.columns.intersection(df2.columns)

    # Subset DataFrames to common columns
    df1_aligned = df1[common_columns]
    df2_aligned = df2[common_columns]

    # Combine the DataFrames
    combined = pd.concat([df1_aligned, df2_aligned], axis=0)

    # Reset index and add 'ID' column
    combined.reset_index(inplace=True)
    combined.rename(columns={"index": "ID"}, inplace=True)
    combined.set_index("ID", inplace=True)

    return combined

def calculate_single_asv(asv_index, metabolomics, asvs, method="pearson"):
    """
    Compute correlations for a single ASV column.
    """

    if method not in ["pearson", "spearman"]:
        raise ValueError("Method must be 'pearson' or 'spearman'.")

    # Select the appropriate correlation function
    corr_func = pearsonr if method == "pearson" else spearmanr

     # Compute correlations
    correlations = np.array([corr_func(asvs[:, asv_index], metabolomics[:, j]) 
                     for j in range(metabolomics.shape[1])
                     ])
    
    # Extract correlation coefficients and p-values
    correlation_coefficients = correlations[:, 0]
    p_values = correlations[:, 1]

    # Calculate R^2 values
    r_squared = correlation_coefficients ** 2
    
    # Apply FDR correction to p-values
    _, fdr_corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
    
    # Combine results into a single array
    result = np.column_stack((correlation_coefficients, p_values, fdr_corrected_p_values, r_squared))

    return result


def calculate_correlations_parallel(df, metabolome_ft, genome_ft, num_workers=4):
    """
    Faster correlation calculation using parallel processing.
    """
    method = st.session_state.get("method", "pearson") 

    transposed_df = df.T
    length_metabolome = metabolome_ft.shape[0]
    length_genome = genome_ft.shape[0]

    metabolome_names = metabolome_ft.index
    asv_names = genome_ft.index

    # Extract metabolomics and ASV columns
    metabolomics = transposed_df.iloc[:, :length_metabolome].values
    asvs = transposed_df.iloc[:, length_metabolome:length_metabolome + length_genome].values
    
    # Use multiprocessing to calculate correlations in parallel
    with Pool(processes=num_workers) as pool:
        results = pool.starmap(calculate_single_asv, 
                               [(i, metabolomics, asvs, method) for i in range(asvs.shape[1])])
        
    results_with_indices = {
        asv_names[i]: pd.DataFrame(
            results[i],
            index=metabolome_names,
            columns=["Estimate", "P-value", "BH-Corrected P-Value", "R2"]
        )
        for i in range(len(asv_names))
    }

    return results_with_indices

def merge_asv_correlation_results(results, target_df, genome_ft):
    """
    Merge ASV correlation results into a dictionary with dataframes having the same rows as target_df.

    Parameters:
    results (dict): Dictionary of correlation results for each ASV.
    target_df (pd.DataFrame): The combined metabolomics and genomics dataframe.
    gen_ft (pd.DataFrame): Genomics feature table.

    Returns:
    dict: A dictionary of merged dataframes for each ASV.
    """
    merged_list = {}

    empty_asv_df = pd.DataFrame(
            0,
            index=genome_ft.index,
            columns=["Estimate", "P-value", "BH-Corrected P-Value", "R2"]
        )

    # Loop through each ASV and merge with target_df
    for asv_name in genome_ft.index:
        #get each df in the results dictionary
        df = results[asv_name]
        
        #combine the empty_asv_df to that
        merged_pscores_df = empty_asv_df.combine_first(df)
        
        # Merge with target_df to include all sample columns
        merged_final_df = pd.concat([merged_pscores_df, target_df], axis=1)

        # Add the merged dataframe to the dictionary
        merged_list[asv_name] = merged_final_df

    return merged_list


def melt_correlation_results(results):
    """
    Melt a dictionary of correlation dataframes into a single dataframe.

    Parameters:
    results (dict): Dictionary where keys are ASV names, and values are dataframes 
                    with correlation results.

    Returns:
    pd.DataFrame: A single dataframe with all correlation results.
    """
    melted_results = []

    for asv_name, df in results.items():
        # Add Feature and Variable columns
        df = df.copy()
        df["Feature"] = df.index  # Index as Feature
        df["Variable"] = asv_name  # ASV name as Variable

        # Rearrange columns so Feature and Variable are first
        df = df[["Feature", "Variable"] + list(df.columns[:-2])]  # Reorder columns

        # Append to list
        melted_results.append(df)

    # Concatenate all dataframes
    final_df = pd.concat(melted_results, axis=0, ignore_index=True)

    return final_df

# Function to estimate the size of the DataFrame in MB
def estimate_df_size(df):
    return df.memory_usage(deep=True).sum() / (1024 * 1024)  # Convert to MB

# Function to convert DataFrame to CSV
@st.cache_data
def convert_df_to_csv(df):
    output = io.BytesIO()
    df.to_csv(output, index=False)
    return output.getvalue()

def estimate_run_time(metabolome_df, genome_df):
    """
    Estimate the time to complete the correlation run based on the number of correlations.
    """
    correlations = metabolome_df.shape[0] * genome_df.shape[0]

    # Create the base message
    message = f"🔹 Estimated correlations: <b>{correlations:,}</b><br><br>"

    # Styling for the message box (centered content)
    box_style = """
        <div style="
            background-color: #f0f2f6; 
            padding: 20px; 
            border-radius: 10px; 
            text-align: center; 
            border: 1px solid #d9d9d9; 
            box-shadow: 2px 2px 5px rgba(0,0,0,0.1);
        ">
    """

    # Estimate time and adjust box color based on severity
    if correlations <= 50000:
        message += "⏳ It will take around <b>30s</b> to complete the run."
        color = "#e6f7ff"  # Light blue
    elif correlations <= 100000:
        message += "⏳ It will take around <b>1 min</b> to complete the run."
        color = "#e6f7ff"
    elif correlations <= 500000:
        message += "⏳ It will take around <b>2 mins</b> to complete the run."
        color = "#e6f7ff"
    elif correlations <= 1000000:
        message += "⏳ It will take around <b>4 mins</b> to complete the run."
        color = "#fffbe6"  # Light yellow
    elif correlations <= 5000000:
        message += "⏳ It will take around <b>12 mins</b> to complete the run. ☕ Go get a coffee!"
        color = "#fffbe6"
    elif correlations <= 10000000:
        message += "⏳ It will take around <b>20 mins</b> to complete the run. ☕ Go get a coffee!"
        color = "#fffbe6"
    else:
        message += "⏳ The run time might exceed <b>20 mins</b>."
        color = "#fff1f0"  # Light red for very long times

    # Final HTML message
    styled_message = f"""
        <div style="
            background-color: {color}; 
            padding: 20px; 
            border-radius: 10px; 
            text-align: center; 
            border: 1px solid #d9d9d9; 
            box-shadow: 2px 2px 5px rgba(0,0,0,0.1);
        ">
            {message}
        </div>
    """

    # Display the message using st.markdown with HTML support
    st.markdown(styled_message, unsafe_allow_html=True)

    return correlations


def centered_button(label, key="centered_button", disabled=False):
    """
    A simple function to display a centered button.
    
    Parameters:
    - label (str): The text displayed on the button.
    - key (str): A unique key to identify the button.
    - disabled (bool): Whether the button should be disabled.
    
    Returns:
    - bool: True if the button is clicked, False otherwise.
    """

    # Basic CSS to center the button (no extra styling)
    st.markdown("""
        <style>
        .center-container {
            display: flex;
            justify-content: center;
            margin: 20px 0;
        }
        </style>
    """, unsafe_allow_html=True)

    # Wrap the button inside a simple centered div
    st.markdown('<div class="center-container">', unsafe_allow_html=True)
    button_clicked = st.button(label, key=key, disabled=disabled)
    st.markdown('</div>', unsafe_allow_html=True)

    return button_clicked



#Set a Memory Limit for Correlation Computation
def check_memory_limit(correlations):
    available_gb = psutil.virtual_memory().available / (1024 ** 3)  # Convert bytes to GB

     # Define dynamic memory requirements
    if correlations <= 50000:
        required_gb = 2  # Small jobs need 2GB
    elif correlations <= 500000:
        required_gb = 4  # Medium jobs need 4GB
    elif correlations <= 2000000:
        required_gb = 6  # Large jobs need 6GB
    else:
        required_gb = 8  # Very large jobs need 8GB+

    if available_gb < required_gb:
        st.error(f"Not enough memory! Required: {required_gb} GB, Available: {available_gb:.2f} GB.")
        return False
    return True

def shuffle_samples(df):
    shuffled_df = df.copy()
    for col in df.columns:
        shuffled_df[col] = np.random.permutation(df[col].values)
    return shuffled_df

def double_shuffle(df):
    shuffled_df = df.apply(np.random.permutation, axis=0)
    shuffled_df = shuffled_df.sample(frac=1, axis=1)  # Shuffle columns
    return shuffled_df

def generate_decoy(target_df):
    decoy_df = target_df.copy()

    # Shuffle each column independently
    for col in decoy_df.columns:
        decoy_df[col] = np.random.permutation(decoy_df[col].values)

    return decoy_df

def dynamic_winsorize(series, factor=3):
    mean = np.mean(series)
    std = np.std(series)
    lower = mean - factor * std
    upper = mean + factor * std
    return np.clip(series, lower, upper)

def filter_low_variance_features(df, threshold=0.01):
    # Calculate variance across samples (axis=1 for rows)
    variances = df.var(axis=1)
    
    # Filter out rows with variance below the threshold
    filtered_df = df[variances > threshold]
    
    return filtered_df


def display_correlation_message(original_edges, filtered_edges):
    """
    Displays a styled message box showing the number of correlations before and after applying FDR cutoffs.
    """

    # Message Content
    message = f"""
        Number of correlations (edges) in the original edge table: <b>{original_edges:,}</b><br><br>
        Number of correlations (edges) after applying FDR Cutoffs: <b>{filtered_edges:,}</b>
    """

    color = "#f0f2f6"  

    # Final Styled Message
    styled_message = f"""
        <div style="
            background-color: {color}; 
            padding: 20px; 
            border-radius: 10px; 
            text-align: center; 
            border: 1px solid #d9d9d9; 
            box-shadow: 2px 2px 5px rgba(0,0,0,0.1);
        ">
            {message}
        </div>
    """

    # Display the Styled Message in Streamlit
    st.markdown(styled_message, unsafe_allow_html=True)
