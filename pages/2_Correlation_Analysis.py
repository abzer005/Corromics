# Import necessary libraries
import streamlit as st
import pandas as pd
from src.common import *        # Importing common functionalities
from src.fileselection import * # Importing file selection functionalities
from src.correlation import *
from src.fdr import *
<<<<<<< HEAD
import random
from multiprocessing import Pool
from datetime import datetime

st.markdown("### Metabolomics-Metagenomics Data Combined")

#st.session_state['metabolome_md']
#st.session_state['genomics_md']

if 'metabolome_ft' in st.session_state and 'genomics_ft' in st.session_state:
    met_ft = st.session_state['metabolome_ft']
    gen_ft = st.session_state['genomics_ft']
    gen_ft = gen_ft[~(gen_ft.eq(0).all(axis=1))]

    # Combine the DataFrames
    target_df = combine_dataframes(met_ft, gen_ft)

    #Create a Decoy set
    random.seed(42)

    # Permute the values of genomics_ft
    decoy_gen_df = gen_ft.apply(np.random.permutation)
    decoy_gen_df += np.random.normal(0, 0.1, decoy_gen_df.shape)  # Add noise
    decoy_df = combine_dataframes(met_ft, decoy_gen_df)

    # Shuffle rows and columns of genomics_ft
    #shuffled_columns = random.sample(list(gen_ft.columns), len(gen_ft.columns))
    #shuffled_rows = random.sample(list(gen_ft.index), len(gen_ft.index))

    #gen_ft_decoy = gen_ft.loc[shuffled_rows, shuffled_columns]
    #decoy_df = combine_dataframes(met_ft, gen_ft_decoy)

    # Display the combined DataFrame in Streamlit
    with st.expander(f"Target Dataframe {target_df.shape}"):
                st.dataframe(target_df)
    with st.expander(f"Decoy Dataframe {decoy_df.shape}"):
                st.dataframe(decoy_df)

    # Perform Correlation ########################################################
    # Record the start time
    # start_time = datetime.now()
    # st.write(f"Start Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    #Give the time message to the user
    estimate_run_time(met_ft, gen_ft)
    
    with st.spinner("Calculating correlations..."):
        target_results = calculate_correlations_parallel(target_df, met_ft, gen_ft)
        decoy_results = calculate_correlations_parallel(decoy_df, met_ft, gen_ft)
        #merged_list = merge_asv_correlation_results(results, target_df, gen_ft

    melted_target = melt_correlation_results(target_results)
    melted_decoy = melt_correlation_results(decoy_results)

    st.session_state['Target_scores'] = melted_target
    st.session_state['Decoy_scores'] = melted_decoy

    with st.expander(f"Correlation Scores of Target Dataframe {melted_target.shape}"):
                st.dataframe(melted_target)
    with st.expander(f"Correlation Scores of Decoy Dataframe {melted_decoy.shape}"):
                st.dataframe(melted_decoy)
    
    st.markdown('### False Discovery Rate')

    if 'Target_scores' in st.session_state and 'Decoy_scores' in st.session_state:
        
        target_scores = st.session_state['Target_scores']
        decoy_scores = st.session_state['Decoy_scores']

        st.write('Select the positive and negative cutoffs for the correlation scores based on your FDR-curve')
        overall_fdr_table, fig_histogram, fig_fdr = calculate_fdr(target_scores, decoy_scores, score_range=(-1, 1), bin_size=0.001)

        st.plotly_chart(fig_histogram)
        st.plotly_chart(fig_fdr)

        # # Record the start time
        # end_time = datetime.now()
        # st.write(f"Start Time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")




else:
    st.warning("Please input the data in the first page to continue the analysis here")
=======
from src.molnet import *
import random
from multiprocessing import Pool
from datetime import datetime
import threading
import time
from scipy.stats import mstats
from streamlit_extras.stylable_container import stylable_container

st.markdown("## Metabolomics-Metagenomics Data Combined")

# Rearrange and display the table
if 'taxonomic_order' in st.session_state and 'rearranged_omics_table' in st.session_state:
    taxa_order = st.session_state['taxonomic_order']
    omics_rearranged = st.session_state['rearranged_omics_table']

    # Allow the user to select the taxonomic level from the available options
    selected_level = st.selectbox("Select a taxonomic level to bin the data:", taxa_order)

    # Perform binning based on the selected level
    binned_by_level = bin_by_taxonomic_level(omics_rearranged, selected_level)

    with st.expander(f"Binned Data at Level: {selected_level} , Original Dimension: {binned_by_level.shape}"):
        st.dataframe(binned_by_level)

    # Get the maximum value of 'Overall_sum' for dynamic thresholding
    max_overall_sum = binned_by_level['Overall_sum'].max()

    # Create two columns for user input
    filter_col1, filter_col2 = st.columns(2)

    # User inputs for filters
    with filter_col1:
        # Set a threshold to filter out low-read/intensity values
        st.session_state['filter_threshold_below'] = st.number_input(
            "**Filter out data with reads/intensities below the user input. Press Enter to apply.**", 
            min_value=0.0, 
            value=0.0, 
            step=1.0, 
            help="**All variables with a total read/intensity below this threshold will be removed from analysis.**"
        )

    with filter_col2:
        st.session_state['filter_threshold_above'] = st.number_input(
            "**Filter out data with reads/intensities above this value. Press Enter to apply.**", 
            max_value=float(max_overall_sum),         # Dynamically set the maximum possible value
            value= float(max_overall_sum),             # Default value is the current max
            step=1.0, 
            help="**All variables with a total read/intensity ABOVE this threshold will be removed from analysis.**"
        ) 

    # Apply the filter dynamically based on user input
    binned_level_filtered = binned_by_level[
        (binned_by_level['Overall_sum'] > st.session_state['filter_threshold_below']) & 
        (binned_by_level['Overall_sum'] <= st.session_state['filter_threshold_above'])
        ]
    
    # Specify the columns to exclude
    exclude_cols = ['index', 'Overall_sum']  # Replace with actual column names

    # Checkboxes for selecting transformations
    apply_imputation = st.checkbox("Apply Imputation (Replace 0s with 1)")
    apply_log_transform = st.checkbox("Apply Log Transformation (log10)")

    # Step 1: Apply Imputation (if selected)
    if apply_imputation:
        for col in binned_level_filtered.columns:
            if col not in exclude_cols:
                binned_level_filtered[col] = binned_level_filtered[col].replace(0, 1)

    # Step 2: Apply Log Transformation (if selected)
    if apply_log_transform:
        # Ensure no zeros before log transformation
        for col in binned_level_filtered.columns:
            if col not in exclude_cols:
                binned_level_filtered[col] = binned_level_filtered[col].replace(0, 1)
                binned_level_filtered[col] = np.log10(binned_level_filtered[col])

    # Step 3: Recalculate 'Overall_Sum' after transformations
    numeric_cols = [col for col in binned_level_filtered.columns if col not in exclude_cols]
    binned_level_filtered['Overall_sum'] = binned_level_filtered[numeric_cols].sum(axis=1)

    st.session_state['binned_omics_table'] = binned_level_filtered

    # Display the binned table
    with st.expander(f"Binned Data at Level: {selected_level} , Filtered Dimension: {binned_level_filtered.shape}"):
        st.dataframe(binned_level_filtered)

   
#####################################################################################
if 'metabolome_ft' in st.session_state and 'binned_omics_table' in st.session_state:
    met_ft = st.session_state['metabolome_ft']
    gen_ft = st.session_state['binned_omics_table']
    gen_ft = gen_ft[~(gen_ft.eq(0).all(axis=1))] #remove rows with all zeros

    # Combine the DataFrames
    target_df = combine_dataframes(met_ft, gen_ft)
    
    #Getting the node table as well
    node_met = met_ft.copy()
    node_gen = gen_ft.copy()

    node_met['Node_Info'] = 'Omics_1'
    node_met['Original_index'] = node_met.index.copy()

    node_met.index = node_met.index.to_series().str.extract(r'(^\d+)', expand=False)
    
    node_gen['Node_Info'] = 'Omics_2'
    node_gen['Original_index'] = node_gen.index.copy()

    node_table = combine_dataframes(node_met, node_gen)
    st.session_state['node_table'] = node_table

    # Decoy generation
    random.seed(42)

    decoy_gen_df = generate_decoy(gen_ft)
    
    # Permute the values of genomics_ft
    #decoy_gen_df = gen_ft.apply(lambda x: np.random.permutation(x), axis=0)
    decoy_df = combine_dataframes(met_ft, decoy_gen_df)

    st.session_state['target_omics'] = gen_ft
    st.session_state['decoy_omics'] = decoy_gen_df

    # Display the combined DataFrame in Streamlit
    with st.expander(f"Target Dataframe {target_df.shape}"):
        st.dataframe(target_df)

    with st.expander(f"Decoy Dataframe {decoy_df.shape}"):
        st.dataframe(decoy_df)

    st.session_state['target_dataframe'] = target_df
    st.session_state['decoy_dataframe'] = decoy_df

    #Give the time message to the user
    estimate_run_time(met_ft, gen_ft)

else:
    st.warning("Please input the data in the first page to continue the analysis here")

######## Run Correlation Analysis -------------------------------------------------------

# Custom-Styled "Run Correlation Analysis" Button
with stylable_container(
    key="run_correlation_container",
    css_styles="""
    button {
        background-color: #00689E;  /* Blue */
        color: white;
        padding: 14px 20px;
        margin: 10px 0;
        border: none;
        border-radius: 10px;
        cursor: pointer;
        width: 100%;
        font-size: 18px;
        font-weight: bold;
        transition: background-color 0.3s ease, transform 0.2s ease;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.2);
    }
     button:hover {
        background-color: #003043;  /* Darker Blue on Hover */
        transform: scale(1.03);     /* Slight zoom effect */
    }
    button:disabled {
        background-color: #A9A9A9;  /* Grey when disabled */
        cursor: not-allowed;
    }
    """,
):
    run_corr_clicked = st.button("Run Correlation Analysis", key="run_correlation", disabled=st.session_state.get("processing", False))

# Perform Correlation Analysis Logic
if all(key in st.session_state for key in ['target_dataframe', 'decoy_dataframe', 'metabolome_ft', 'binned_omics_table', 'target_omics', 'decoy_omics']):
    if run_corr_clicked:
        st.session_state["processing"] = True  # Set processing flag

        with st.spinner("Calculating correlations..."):
            st.session_state['target_results'] = calculate_correlations_parallel(
                st.session_state['target_dataframe'], 
                st.session_state['metabolome_ft'], 
                st.session_state['target_omics']
            )
            st.session_state['decoy_results'] = calculate_correlations_parallel(
                st.session_state['decoy_dataframe'], 
                st.session_state['metabolome_ft'], 
                st.session_state['decoy_omics']
            )

        st.session_state["processing"] = False  # Reset flag when done
        st.success("✅ Correlation analysis completed!")

else:
    st.warning("⚠️ Please input the data first to continue the analysis.")

# Show a processing message if the analysis is running
if st.session_state.get("processing", False):
     with st.empty():
        while st.session_state.get("processing", False):
             st.warning("Processing... Please wait while correlations are being calculated.")
             time.sleep(1)


# Check if correlation results are available before processing
if "target_results" in st.session_state and "decoy_results" in st.session_state:
    
    # Melt results only after they exist
    melted_target = melt_correlation_results(st.session_state["target_results"])
    melted_decoy = melt_correlation_results(st.session_state["decoy_results"])

    # Store melted results in session state
    st.session_state["Target_scores"] = melted_target
    st.session_state["Decoy_scores"] = melted_decoy

    # Estimate DataFrame size
    df_target_size = estimate_df_size(melted_target)
    df_decoy_size = estimate_df_size(melted_decoy)

    target_csv_data = convert_df_to_csv(melted_target)
    decoy_csv_data = convert_df_to_csv(melted_decoy)

    # Display or provide download option based on size
    if df_target_size <= 200:
        with st.expander(f"Correlation Scores of Target Dataframe {melted_target.shape}"):
            st.dataframe(melted_target)
            st.download_button(
                label="Download Target DataFrame as CSV",
                data=target_csv_data,
                file_name="target_correlation_results.csv",
                mime="text/csv"
            )
    else:
        st.warning(f"⚠️ The target correlation dataframe is too large to display (**{df_target_size:.1f} MB**). Download it instead.")
        
        st.download_button(
            label="Download Target Correlation Results (CSV)",
            data=target_csv_data,
            file_name="target_correlation_results.csv",
            mime="text/csv"
        )

    if df_decoy_size <= 200:
        with st.expander(f"Correlation Scores of Decoy Dataframe {melted_decoy.shape}"):
            st.dataframe(melted_decoy)
            st.download_button(
                label="Download Decoy Correlation Results (CSV)",
                data=decoy_csv_data,
                file_name="decoy_correlation_results.csv",
                mime="text/csv"
            )
    else:
        st.warning(f"⚠️ The decoy correlation dataframe is too large to display (**{df_decoy_size:.1f} MB**). Download it instead.")

        st.download_button(
            label="Download Decoy Correlation Results (CSV)",
            data=decoy_csv_data,
            file_name="decoy_correlation_results.csv",
            mime="text/csv"
        )

# Add space before FDR analysis
st.divider()
st.markdown('## False Discovery Rate')

# Custom-Styled "FDR" Button
with stylable_container(
    key="run_fdr_container",
    css_styles="""
    button {
        background-color: #00689E;  /* Blue */
        color: white;
        padding: 14px 20px;
        margin: 10px 0;
        border: none;
        border-radius: 10px;
        cursor: pointer;
        width: 100%;
        font-size: 18px;
        font-weight: bold;
        transition: background-color 0.3s ease, transform 0.2s ease;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.2);
    }
     button:hover {
        background-color: #003043;  /* Darker Blue on Hover */
        transform: scale(1.03);     /* Slight zoom effect */
    }
    button:disabled {
        background-color: #A9A9A9;  /* Grey when disabled */
        cursor: not-allowed;
    }
    """,
):
    run_fdr_button_clicked = st.button("Apply FDR", key="run_fdr")


# Initialize session state for button click
if "run_fdr_clicked" not in st.session_state:
    st.session_state["run_fdr_clicked"] = False

# Initialize session state for cutoffs if not already set
if "neg_fdr_cutoff" not in st.session_state:
    st.session_state["neg_fdr_cutoff"] = -0.5  # Default value

if "pos_fdr_cutoff" not in st.session_state:
    st.session_state["pos_fdr_cutoff"] = 0.5  # Default value

if "filtered_target_csv" not in st.session_state:
    st.session_state["filtered_target_csv"] = None


# Ensure FDR analysis only runs if target & decoy scores exists
if 'Target_scores' in st.session_state and 'Decoy_scores' in st.session_state:

    if run_fdr_button_clicked:
        st.session_state["run_fdr_clicked"] = True  # Store button state         
        target_scores = st.session_state['Target_scores']
        decoy_scores = st.session_state['Decoy_scores']

        # Apply Winsorization only to the 'Estimate' column
        #target_scores['Estimate'] = dynamic_winsorize(target_scores['Estimate'])
        #decoy_scores['Estimate'] = dynamic_winsorize(decoy_scores['Estimate'])

        # Step 1: Calculate FDR Data
        overall_fdr_table = calculate_fdr(target_scores, 
                                        decoy_scores,
                                        score_range=(-1, 1), 
                                        bin_size=0.001)
        fig_fdr, fig_histogram = plot_fdr_figures(overall_fdr_table, target_scores, decoy_scores)

        st.session_state["fig_histogram"] = fig_histogram
        st.session_state["fig_fdr"] = fig_fdr

else:
    st.warning("Please run FDR to continue here")
    
if st.session_state["run_fdr_clicked"]:

    st.write('Select the positive and negative cutoffs for the correlation scores based on your FDR-curve')
    st.plotly_chart(st.session_state["fig_histogram"])
    st.plotly_chart(st.session_state["fig_fdr"])

    ######USER DEFINED CUTOFFS----------------------
    # Create two columns for user input
    input_col1, input_col2 = st.columns(2)

    # User inputs for FDR cutoff
    with input_col1:
         # User input for filter threshold (default is 0)
        st.session_state["neg_fdr_cutoff"] = st.number_input("Enter Negative FDR Cutoff (e.g., -0.5):",
                                                             value= st.session_state["neg_fdr_cutoff"],
                                                             step=0.05)

    with input_col2:
        st.session_state["pos_fdr_cutoff"] = st.number_input("Enter Positive FDR Cutoff (e.g., 0.5):",
                                                             value = st.session_state["pos_fdr_cutoff"], 
                                                             step = 0.05)

    # Subset the target dataframe based on user-defined cutoffs
    filtered_target = melted_target[(melted_target['Estimate'] <= st.session_state["neg_fdr_cutoff"]) | 
                                        (melted_target['Estimate'] >= st.session_state["pos_fdr_cutoff"])]
        
    st.session_state['filtered_target'] = filtered_target
    st.session_state['filtered_target_csv'] = filtered_target.to_csv(index=False).encode('utf-8')

    # Display the filtered data
    if st.session_state["filtered_target"] is not None:

        st.markdown('#### Filtered Correlation Scores Based on FDR Cutoffs')
    
        display_correlation_message(melted_target.shape[0], filtered_target.shape[0])
        st.write('')
        if st.session_state["filtered_target_csv"] is not None:
            graphml_data = "This is where your GraphML data would go."
                    
            # Create three columns in the UI for buttons
            c1, c2, c3,c4 = st.columns([1,1,1,1]) # to make the buttons lie close to each other
                        
            # Column 2: Download CSV
            with c2:
                # Provide a download button for CSV
                st.download_button(
                    label="Download Results as CSV",
                    data=st.session_state['filtered_target_csv'],  
                    file_name="filtered_correlations_results.csv",
                    mime='text/csv'
                )

            # Column 3: Download GraphML
            with c3:
                # Assuming node_table and edge_table are loaded into Streamlit session state
                output_file = "correlation_network.graphml"

                if st.button("Generate GraphML"):
                    result = generate_graphml_corromics(st.session_state['node_table'], 
                                                        st.session_state['filtered_target'], 
                                                        output_file)
                    if result:
                        with open(output_file, "rb") as f:
                            st.download_button("Download GraphML", 
                                               f, 
                                               file_name=output_file, 
                                               mime="application/graphml+xml")

                # st.download_button(
                #     label="Download network in graphml",
                #     data=graphml_data,  # Replace with actual GraphML data
                #     file_name="network.graphml",
                #     mime='application/graphml+xml',
                #     key="download_graphml"
                # )

                # Button to trigger GraphML generation
                #if __name__ == '__main__':
                #     generate_graphml_zip_corromics()  # Function to generate and download GraphML file
else:
    st.warning("Please run the FDR before further processing")


# Custom-Styled "FDR" Button
with stylable_container(
    key="fdr_reset_container",
    css_styles="""
    button {
        background-color: #FF000066;  /* Orange */
        color: black;
        padding: 14px 20px;
        margin: 10px 0;
        border: none;
        border-radius: 10px;
        cursor: pointer;
        width: 100%;
        font-size: 18px;
        font-weight: bold;
        transition: background-color 0.3s ease, transform 0.2s ease;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.2);
    }
     button:hover {
        background-color: #F0F0F0;  /* Blue on Hover */
        transform: scale(1.05);     /* Slight zoom effect */
        color: black;
    }
    button:disabled {
        background-color: #A9A9A9;  /* Grey when disabled */
        cursor: not-allowed;
    }
    """,
):
    reset_fdr_clicked = st.button("Reset FDR Analysis", key="reset_fdr")








# ✅ Add a reset button to clear the session state after results are shown
if st.session_state["run_fdr_clicked"]:
    if reset_fdr_clicked:
        st.session_state["run_fdr_clicked"] = False  # Reset state
        st.rerun()  # Force Streamlit to refresh
>>>>>>> b669908 (FDR added)
