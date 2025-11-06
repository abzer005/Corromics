import streamlit as st
import numpy as np
from src.common import *     
from src.fileselection import *
from src.correlation import *
from src.fdr import *
from src.molnet import *
from streamlit_extras.stylable_container import stylable_container
import random
import time

page_setup()
initialize_app()
st.session_state['max_corr'] = get_max_correlation_limit()

##########################################################################

st.markdown("## Omics Data Combined")

if ('rearranged_omics_table' in st.session_state and not st.session_state['rearranged_omics_table'].empty):
    omics_rearranged = st.session_state['rearranged_omics_table']

    if (st.session_state.get("has_taxonomic_info") == "Yes"and st.session_state.get("taxonomic_order")):
        taxa_order = st.session_state['taxonomic_order']

        st.warning(
            "⚠️ By default, the data is binned at the **highest available hierarchical level**."
            " This occurs automatically when the page first loads. Please review the selected level below and change it if needed."
            )

        # Allow the user to select the taxonomic level from the available options
        selected_level = st.selectbox("Select a taxonomic level to bin the data:", taxa_order)

        # Perform binning based on the selected level
        binned_by_level = bin_by_taxonomic_level(omics_rearranged, selected_level)
    
    elif (not st.session_state.get("taxonomic_order") or st.session_state.get("has_taxonomic_info") == "No"):
 
        st.warning("No taxonomic order available for binning. The index column was used as the binning level.")
        # Drop rows and cols where all values are NaN
        omics_numeric = omics_rearranged.select_dtypes(include='number')
        selected_level = omics_numeric.index.name or "Index"

        # Allow the user to select the taxonomic level from the available options
        binned_by_level = bin_by_flat_id(omics_numeric)      
    
    else:
        st.warning("❗ Could not determine binning strategy. Please check your inputs.")
        binned_by_level = None
        selected_level = None
   
    if binned_by_level is not None:
        with st.expander(f"Binned Data at Level: {selected_level}, Original Dimension: {binned_by_level.shape}"):
            st.dataframe(binned_by_level)
            st.info(
                "👉 Scroll to the far right of the table to see the **Overall_sum** column. "
                "This represents the total abundance or count for each row and can be used with the filter options below to refine your data."
                )     
    
    st.markdown('#### Filter Binned Data by intensity')
    st.info(
        "Filtering the data by intensity allows users to reduce the number of features included in the correlation analysis, "
        "which can **substantially affect both the number of computed correlations and the overall runtime**. "
        "By default, features with an **Overall Sum of 0** are excluded, as they result in NA values.\n\n"
        "The filtering is optional and intended solely for computational or exploratory purposes. "
        "**No assumptions are made about the biological relevance of features excluded through these filters.**"
    )

    
    binned_level_filtered, exclude_cols = filter_by_overall_sum(binned_by_level, label_prefix="table2_")

    binned_level_filtered, method_used = apply_transformation(
        binned_level_filtered,
        exclude_cols=["feature_ID"],
        key_prefix="binned"
        )

    # Recalculate 'Overall_Sum' after transformations
    numeric_cols = [col for col in binned_level_filtered.columns if col not in exclude_cols]
    binned_level_filtered['Overall_sum'] = binned_level_filtered[numeric_cols].sum(axis=1)

    st.session_state['binned_omics_table'] = binned_level_filtered

    # Display the binned table
    with st.expander(f"Binned Data at Level: {selected_level}, Filtered Dimension: {binned_level_filtered.shape}"):
        st.dataframe(binned_level_filtered)

else:
    st.warning("Please upload the data in the first page to continue the analysis here")

st.markdown('---')

st.write("## Target vs Decoy Dataframe")   
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
    st.session_state['no_correlations'] = estimate_run_time(met_ft, gen_ft)

else:
    st.warning("Please input the data in the first page to continue the analysis here")

######## Run Correlation Analysis -------------------------------------------------------
st.markdown("<br>", unsafe_allow_html=True)
st.session_state["method"] = st.radio(
    "Choose correlation method",
    options=["pearson", "spearman"],
    format_func=lambda x: "Pearson (linear relationship)" if x == "pearson" else "Spearman (monotonic relationship)",
    help="""
    **Pearson**: Measures linear correlation. Assumes normal distribution.\n
    **Spearman**: Non-parametric, measures monotonic relationships. Useful when data is not normally distributed.
    """,
    horizontal=True
)

# Custom-Styled "Run Correlation Analysis" Button
with stylable_container(
    key="run_correlation_container",
    css_styles="""
    button { background-color: #00689E; color: white; font-weight: bold; width: 100%; }
    button:hover { background-color: #003043; transform: scale(1.05); }
    button:disabled { background-color: #A9A9A9;  cursor: not-allowed;}
    """):
    run_corr_clicked = st.button("Run Correlation Analysis", 
                                 key="run_correlation", 
                                 disabled=st.session_state.get('no_correlations', 0) >= st.session_state.max_corr)

# Perform Correlation Analysis Logic
if all(key in st.session_state for key in ['target_dataframe', 'decoy_dataframe', 'metabolome_ft', 'binned_omics_table', 'target_omics', 'decoy_omics']):
    
    # Estimate total number of correlations
    correlations = st.session_state['no_correlations']

    if run_corr_clicked:

        if not st.session_state.is_restricted_mode and correlations >= 10_000_000:
            st.error(
                "❌ Too many correlations to compute in this environment — more than **10 million**.\n\n"
                "This could crash your local memory and force you to quit other apps.💡 Please reduce the number of features or samples to perform a meaningful analysis and avoid memory overload."
            )
        elif st.session_state.is_restricted_mode and correlations >= 1_000_000:
            st.error(f"❌ Too many correlations to compute in this environment (≥ {st.session_state.max_corr}). Please clone the app and run it locally. This helps avoid memory crashes in the cloud environment.")
            st.stop()
        else:
            st.session_state["processing"] = True  # Set processing flag

            # Placeholders for live progress
            status_text = st.empty()
            progress_bar = st.progress(0.0)

            def corr_progress(done, total, est_left):
                fraction = done / total if total else 0.0
                progress_bar.progress(fraction)

                est_left_sec = int(est_left)
                mins, secs = divmod(est_left_sec, 60)
                if mins > 0:
                    time_str = f"~{mins} min {secs} s left"
                else:
                    time_str = f"~{secs} s left"

                status_text.info(
                    f"Calculating correlations: **{done:,}/{total:,}** "
                    f"({fraction:.1%}): {time_str}"
                    )


            with st.spinner("Calculating correlations..."):
                st.session_state['target_results'] = calculate_correlations_parallel(
                    st.session_state['target_dataframe'], 
                    st.session_state['metabolome_ft'], 
                    st.session_state['target_omics'],
                    num_workers=4,
                    progress_callback=corr_progress,
                )
                st.session_state['decoy_results'] = calculate_correlations_parallel(
                    st.session_state['decoy_dataframe'], 
                    st.session_state['metabolome_ft'], 
                    st.session_state['decoy_omics'],
                    num_workers=4,
                    progress_callback=corr_progress,
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
        st.warning(f"⚠️ The correlation dataframes are too large to display (**{df_target_size:.1f} MB**). Download them instead.")
        
        corr_c1, corr_c2 = st.columns(2)
        with corr_c1:
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
        if 'corr_c1' not in locals():  # use same row if available
            corr_c1, corr_c2 = st.columns(2)
        with corr_c2:
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
        button { background-color: #00689E; color: white; font-weight: bold; width: 100%; }
        button:hover { background-color: #003043; transform: scale(1.05); }
        button:disabled { background-color: #A9A9A9;  cursor: not-allowed;}
    """):

    run_fdr_button_clicked = st.button("Apply FDR", 
                                       key="run_fdr",
                                       disabled= st.session_state.get('no_correlations', 0) >= st.session_state.max_corr
                                       )



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
        
        correlations = st.session_state['no_correlations']
        # Estimate total number of correlations
        if correlations >= st.session_state.max_corr and not st.session_state.is_restricted_mode:
            st.error(f"❌ Too many correlations to compute in this environment (≥ {st.session_state.max_corr}). Please clone or download the app and run it locally. This helps avoid memory crashes in the cloud environment.")
            st.stop()
        else:       
            st.session_state["run_fdr_clicked"] = True

            target_scores = st.session_state['Target_scores']
            decoy_scores = st.session_state['Decoy_scores']

            # Step 1: Calculate FDR Data
            overall_fdr_table = calculate_fdr(target_scores, 
                                            decoy_scores,
                                            score_range=(-1, 1), 
                                            bin_size=0.001)
            fig_fdr, fig_histogram = plot_fdr_figures(overall_fdr_table, target_scores, decoy_scores)

            st.session_state["fig_histogram"] = fig_histogram
            st.session_state["fig_fdr"] = fig_fdr
else:
    st.warning("Please run correlation analysis to continue this step.")


# Display FDR results only if button was clicked and correlations are within limits
if (
    st.session_state.get("run_fdr_clicked", False)
    and (
        st.session_state.get("no_correlations", 0) < st.session_state.max_corr
        or not st.session_state.is_restricted_mode
    )
):
    st.plotly_chart(st.session_state["fig_histogram"])
    st.plotly_chart(st.session_state["fig_fdr"])

    melted_target = st.session_state['Target_scores']

    st.write('Select the positive and negative cutoffs for the correlation scores based on your FDR-curve')

    ######USER DEFINED CUTOFFS----------------------
    # Create two columns for user input
    input_col1, input_col2 = st.columns(2)

    # User inputs for FDR cutoff
    with input_col1:
         # User input for filter threshold (default is 0)
        st.session_state["neg_fdr_cutoff"] = st.number_input("Enter Negative FDR Cutoff (e.g., -0.5):",
                                                             value=st.session_state.get("neg_fdr_cutoff", -0.5),
                                                             min_value=-1.0,
                                                             max_value=0.0,
                                                             step=0.05
                                                         )


    with input_col2:
        st.session_state["pos_fdr_cutoff"] = st.number_input("Enter Positive FDR Cutoff (e.g., 0.5):",
                                                             value = st.session_state.get("pos_fdr_cutoff", 0.5),
                                                             min_value=0.0,
                                                             max_value=1.0, 
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
                # Generate GraphML file
                if 'node_table' in st.session_state:
                    graphml_ready = generate_graphml_corromics(st.session_state['node_table'], 
                                                        st.session_state['filtered_target'], 
                                                        output_file)
                    if graphml_ready:
                        with open(output_file, "rb") as f:
                            st.download_button("Download GraphML", 
                                            f, 
                                            file_name=output_file, 
                                            mime="application/graphml+xml")
                    else:
                        st.error("❌ Failed to generate GraphML.")
                else:
                    st.error("❌ Node table not found in session state. Please run the correlation analysis first.")


elif st.session_state.get("run_fdr_clicked", False) and st.session_state.get("no_correlations", 0) >= st.session_state.max_corr and st.session_state.is_restricted_mode:
    st.error(f"❌ As correlations exceed {st.session_state.max_corr}, no correlations were computed for this level. 💡 Please clone or download the app and run it locally. This helps avoid memory crashes in the cloud environment.")
    st.stop()
        
# Custom-Styled "FDR" Button
with stylable_container(
    key="fdr_reset_container",
    css_styles="""
        button {background-color: #FF000066; color: black; font-weight: bold; width: 100%;}
        button:hover { background-color: #003043; transform: scale(1.05); }
        button:disabled { background-color: #A9A9A9;  cursor: not-allowed;}
    """):
    reset_fdr_clicked = st.button("Reset FDR Analysis", 
                                  key="reset_fdr"
                                  )


# ✅ Add a reset button to clear the session state after results are shown
if st.session_state["run_fdr_clicked"]:
    if reset_fdr_clicked:
        st.session_state["run_fdr_clicked"] = False  # Reset state
        st.session_state["neg_fdr_cutoff"] = -0.5  # Reset to default
        st.session_state["pos_fdr_cutoff"] = 0.5  # Reset to default
        st.session_state["filtered_target"] = None  # Reset filtered target
        st.session_state["filtered_target_csv"] = None
        st.rerun()  # Force Streamlit to refresh
