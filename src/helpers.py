import streamlit as st
from src.fileselection import *
from src.correlation import *
from src.fdr import *
from src.molnet import *
from streamlit_extras.stylable_container import stylable_container
import random
import time
# import random





def load_data_section():
    # Assume that file loading and cleaning were done in another page and stored in session_state
    if 'metabolome_ft' not in st.session_state or 'metabolome_md' not in st.session_state:
        st.warning("⚠️ Please upload and process the data in the first page before continuing.")


def binning_section():
    if 'rearranged_omics_table' not in st.session_state or st.session_state['rearranged_omics_table'].empty:
        return None, None

    omics_rearranged = st.session_state['rearranged_omics_table']

    if st.session_state.get("has_taxonomic_info") == "Yes" and st.session_state.get("taxonomic_order"):
        taxa_order = st.session_state['taxonomic_order']

        st.warning(
            "⚠️ By default, the data is binned at the highest hierarchical level. Please review the selected level below and change it if needed."
            )
        
        selected_level = st.selectbox("Select a taxonomic level to bin the data:", taxa_order)
        binned = bin_by_taxonomic_level(omics_rearranged, selected_level)

    elif (not st.session_state.get("taxonomic_order") or st.session_state.get("has_taxonomic_info") == "No"):
        st.warning("No taxonomic structure detected. Binning will be based on feature index.")
        omics_numeric = omics_rearranged.select_dtypes(include='number')
        selected_level = omics_numeric.index.name or "Index"
        binned = bin_by_flat_id(omics_numeric)
    
    else:
        st.warning("❗ Could not determine binning strategy. Please check your inputs.")
        binned = None
        selected_level = None

    if binned is not None:
        with st.expander(f"Binned Data at Level: {selected_level}, Original Shape: {binned.shape}"):
            st.dataframe(binned)
            st.info("👉 Scroll to the far right of the table to see the **Overall_sum** column. This indicates total abundance/count per feature and can be used for filtering")

    return binned, selected_level


def filter_section_second_table(binned_df):
    if binned_df is None:
        return None, []

    st.markdown("#### Filter Binned Data by Intensity")
    st.info(
        "Filtering the data by intensity allows users to reduce the number of features included in the correlation analysis, "
        "which can **substantially affect both the number of computed correlations and the overall runtime**. "
        "By default, features with an **Overall Sum of 0** are excluded to avoid NA correlations.\n\n"
        "All additional filtering is optional and intended solely for computational or exploratory purposes. "
        "**No assumptions are made about the biological relevance of features excluded through these filters.**"
    )

    filtered_df, exclude_cols = filter_by_overall_sum(binned_df, label_prefix="table2_")

    if st.checkbox("Apply Imputation (Replace 0s with 1)", key="imputation_checkbox"):
        for col in filtered_df.columns:
            if col not in exclude_cols:
                filtered_df[col] = filtered_df[col].replace(0, 1)

    if st.checkbox("Apply Log Transformation (log10)", key="log_transform_checkbox"):
        for col in filtered_df.columns:
            if col not in exclude_cols:
                filtered_df[col] = filtered_df[col].replace(0, 1)
                filtered_df[col] = np.log10(filtered_df[col])

    numeric_cols = [col for col in filtered_df.columns if col not in exclude_cols]
    filtered_df['Overall_sum'] = filtered_df[numeric_cols].sum(axis=1)
    st.session_state['binned_omics_table'] = filtered_df

    with st.expander(f"Filtered Omics Table, Final Shape: {filtered_df.shape}"):
        st.dataframe(filtered_df)

    return filtered_df, exclude_cols

def filter_section_first_table(metabolome_df):
    if metabolome_df is None:
        return None, []

    st.markdown("#### Filter Metabolomics Data by intensity")
    num_cols = metabolome_df.select_dtypes(include="number").columns
    
    # Compute row-wise sum and store as new column
    metabolome_df["Overall_sum"] = metabolome_df[num_cols].sum(axis=1)
    final_df, exclude_cols = filter_by_overall_sum(metabolome_df, label_prefix="table1_")
    st.session_state['metabolome_ft'] = final_df

    with st.expander(f"Filtered Metabolomics Table, Final Shape: {final_df.shape}"):
        st.dataframe(final_df)
        st.info(
            "👉 Scroll to the far right of the table to see the **Overall_sum** column. "
            "This represents the total abundance or count for each row and can be used with the filter options above to refine your data."
            )     

    return final_df, exclude_cols


def preview_target_decoy(gen_ft):
    if 'metabolome_ft' not in st.session_state:
        st.warning("Please load the metabolomics data first.")
        return
    
    elif 'metabolome_ft' in st.session_state and 'binned_omics_table' in st.session_state:
        met_ft = st.session_state['metabolome_ft']
        gen_ft = st.session_state['binned_omics_table']
        gen_ft = gen_ft[~(gen_ft.eq(0).all(axis=1))]  # remove empty rows
        st.session_state['no_correlations'] = estimate_run_time(met_ft, gen_ft)
        
        target_df = combine_dataframes(met_ft, gen_ft)
        random.seed(42)
        decoy_gen_df = generate_decoy(gen_ft)
        decoy_df = combine_dataframes(met_ft, decoy_gen_df)

        st.session_state['target_dataframe'] = target_df
        st.session_state['decoy_dataframe'] = decoy_df
        st.session_state['target_omics'] = gen_ft
        st.session_state['decoy_omics'] = decoy_gen_df
        st.write(' ')
        with st.expander(f"Target Dataframe {target_df.shape}"):
            st.dataframe(target_df)
        with st.expander(f"Decoy Dataframe {decoy_df.shape}"):
            st.dataframe(decoy_df)


def correlation_section():

    method = st.radio("Choose correlation method:", ["pearson", "spearman"],
                      format_func=lambda x: "Pearson (linear)" if x == "pearson" else "Spearman (monotonic)",
                      help="""
                      **Pearson**: Measures linear correlation. Assumes normal distribution.\n
                      **Spearman**: Non-parametric, measures monotonic relationships. Useful when data is not normally distributed.
                      """,
                      horizontal=True)
    
    st.session_state['method'] = method

    corr_count = st.session_state.get("no_correlations", 0)
    max_corr = st.session_state['max_corr']
    disable_button = corr_count >= max_corr

    with stylable_container(
        "run_correlation_container", 
        css_styles="""
        button { background-color: #00689E; color: white; font-weight: bold; width: 100%; }
        button:hover { background-color: #003043; transform: scale(1.05); }
        button:disabled { background-color: #A9A9A9;  cursor: not-allowed;}
    """):
        run_corr = st.button("Run Correlation Analysis", key="run_correlation", disabled=disable_button)

    if run_corr:
        if not st.session_state.is_restricted_mode and corr_count >= 10_000_000:
            st.error(
                "❌ Too many correlations to compute in this environment — more than **10 million**.\n\n"
                "This could crash your local memory and force you to quit other apps.💡 Please reduce the number of features or samples to perform a meaningful analysis and avoid memory overload."
            )
            st.stop()
        elif st.session_state.is_restricted_mode and corr_count >= 1_000_000:
            st.error(f"❌ Too many correlations to compute in this environment (≥ {max_corr}). Please clone the app and run it locally. This helps avoid memory crashes in the cloud environment.")
            st.stop()
        
        else:
            st.session_state["processing"] = True  # Set processing flag

            with st.spinner("Calculating correlations..."):
                st.session_state['target_results'] = calculate_correlations_parallel(
                    st.session_state['target_dataframe'], st.session_state['metabolome_ft'], st.session_state['target_omics']
                )
                st.session_state['decoy_results'] = calculate_correlations_parallel(
                    st.session_state['decoy_dataframe'], st.session_state['metabolome_ft'], st.session_state['decoy_omics']
                )    
            st.session_state['processing'] = False
            st.success("✅ Correlation analysis completed!")


def download_raw_correlation_results(target, decoy):
  
    # Melt results only after they exist
    melted_target = melt_correlation_results(target)
    melted_decoy = melt_correlation_results(decoy)

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

####################################################
# Initialize session and set global correlation cap
initialize_app()
st.session_state['max_corr'] = get_max_correlation_limit()

st.markdown("## Metabolomics-Metagenomics Data Combined")

load_data_section()

binned_by_level, selected_level = binning_section()

binned_level_filtered, exclude_cols = filter_section_second_table(binned_by_level)

if ('metabolome_ft' in st.session_state and not st.session_state['rearranged_omics_table'].empty):
    final_ft = st.session_state['metabolome_ft']
    final_ft, met_exclude_cols = filter_section_first_table(final_ft)

preview_target_decoy(binned_level_filtered)

st.markdown("---")
st.write("## Correlation Analysis")

correlation_section()

if st.session_state.get("processing", False):
    with st.empty():
        while st.session_state.get("processing", False):
            st.warning("Processing... Please wait while correlations are being calculated.")
            time.sleep(1)

if "target_results" in st.session_state and "decoy_results" in st.session_state:
    target = st.session_state["target_results"]
    decoy = st.session_state["decoy_results"]
    download_raw_correlation_results(target, decoy)
