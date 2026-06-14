import streamlit as st
import numpy as np
from src.common import *     
from src.fileselection import *
from src.correlation import *
from src.sparcc import *
from src.fdr import *
from src.molnet import *
from src.dist_correlation import *
from src.joint_rpca import *
from src.progress import *
from src.visualization import *
from src.correlation_sections import association_heatmap_section, false_discovery_rate_section
from streamlit_extras.stylable_container import stylable_container
import random
import time


page_setup()
initialize_app()
st.session_state['max_corr'] = get_max_correlation_limit()

##########################################################################

st.markdown("## Hierarchical Binning (Second Omics Dataset)")

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
    st.warning(
        f"Binned Data (at Level: {selected_level}): {binned_by_level.shape[0]} → {binned_level_filtered.shape[0]} after filtering"
    )

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

st.write("## Correlation Estimate & Feature Filtering")   
#####################################################################################

if (
    "metabolome_ft_unfiltered_for_correlation" not in st.session_state
    and "metabolome_ft" in st.session_state
    and st.session_state["metabolome_ft"] is not None
    and not st.session_state["metabolome_ft"].empty
):
    st.session_state["metabolome_ft_unfiltered_for_correlation"] = st.session_state["metabolome_ft"].copy()

if "metabolome_ft_unfiltered_for_correlation" in st.session_state:
    if st.button("Reset Metabolite Feature Filtering", key="reset_metabolite_feature_filtering"):
        st.session_state["metabolome_ft"] = st.session_state["metabolome_ft_unfiltered_for_correlation"].copy()
        st.session_state["filter_metabolome_intensity"] = False
        st.session_state["filter_metabolome_variance"] = False

        for key in [
            "metabolome_ft_filter_below",
            "metabolome_ft_filter_above",
            "metabolome_ft_filter_threshold_below",
            "metabolome_ft_filter_threshold_above",
            "metabolome_ft_var_variance_min",
            "metabolome_ft_var_variance_max",
        ]:
            st.session_state.pop(key, None)

        clear_correlation_outputs()
        st.rerun()

@st.fragment
def correlation_estimate_and_filtering():

    if 'metabolome_ft' in st.session_state and 'binned_omics_table' in st.session_state:
        met_ft = st.session_state['metabolome_ft']
        gen_ft = st.session_state['binned_omics_table']
        gen_ft = gen_ft[~(gen_ft.eq(0).all(axis=1))]

        st.session_state['no_correlations'] = estimate_run_time(met_ft, gen_ft)
        warn_large_correlation_run(
            st.session_state['no_correlations'],
            met_ft.shape[0],
            gen_ft.shape[0],
        )
        st.write("")

        corromics_feature_filtering_info()
        st.write("")

        met_ft_filtered = met_ft.copy()
        met_ft_filtered["Overall_sum"] = met_ft_filtered[numeric_cols].sum(axis=1)

        with st.expander(f"Metabolite Data: {met_ft_filtered.shape[0]} features"):
            st.dataframe(met_ft_filtered)
            st.info(
                "👉 Scroll to the far right of the table to see the **Overall_sum** column. "
                "This represents the total sum for each row and can be used with the filter options."
            )

        # ---------------- INTENSITY FILTER ----------------
        if st.checkbox(
            "Filter features based on overall intensity (sum across samples)",
            value=False,
            key="filter_metabolome_intensity",
            on_change=clear_processing,
        ):
            met_ft_filtered, exclude_cols = filter_by_overall_sum(
                met_ft_filtered,
                label_prefix="metabolome_ft_"
            )

            st.warning(f"Metabolite features after intensity filtering: {met_ft.shape[0]} → {met_ft_filtered.shape[0]}")
            st.session_state['metabolome_ft'] = met_ft_filtered

            with st.expander(f"Filtered Data: {met_ft_filtered.shape[0]} features"):
                st.dataframe(met_ft_filtered)

        # ---------------- VARIANCE FILTER ----------------
        if st.checkbox(
            "Filter features based on variance across samples",
            value=False,
            key="filter_metabolome_variance",
            on_change=clear_processing,
        ):
            numeric_cols_filtered = met_ft_filtered.select_dtypes(include=[np.number]).columns.tolist()

            numeric_cols_filtered = [
                col for col in numeric_cols_filtered
                if col not in ["Overall_sum", "Variance"]
            ]

            n_before_var = met_ft_filtered.shape[0]
            met_ft_filtered["Variance"] = met_ft_filtered[numeric_cols_filtered].var(axis=1)

            with st.expander(f"Metabolite Data with variance column: {met_ft_filtered.shape[0]} features"):
                st.dataframe(
                    met_ft_filtered.style.format({"Variance": "{:.2e}"})
                )

            met_ft_filtered, exclude_cols = filter_by_variance(
                met_ft_filtered,
                label_prefix="metabolome_ft_var_"
            )
            
            st.warning(f"Metabolite features after variance filtering: {n_before_var} → {met_ft_filtered.shape[0]}")

            with st.expander(f"Filtered Metabolite Data: {met_ft_filtered.shape[0]} features"):
                st.dataframe(met_ft_filtered)

        # cleanup
        met_ft_filtered = met_ft_filtered.drop(columns=["Overall_sum", "Variance"], errors="ignore")
        st.session_state['metabolome_ft'] = met_ft_filtered


    else:
        st.warning("Please input the data in the first page to continue the analysis here")

correlation_estimate_and_filtering()


if st.button("Refresh Target/Decoy Dataframes", key="refresh_target_decoy"):
    clear_correlation_outputs()
    st.rerun()


st.write("## Target vs Decoy Dataframe")
if 'metabolome_ft' in st.session_state and 'binned_omics_table' in st.session_state:
    met_ft = st.session_state['metabolome_ft']
    gen_ft = st.session_state['binned_omics_table']

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
    np.random.seed(42)

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

else:
    st.warning("Please input the data in the first page to continue the analysis here")

######## Run Correlation Analysis -------------------------------------------------------
st.markdown("<br>", unsafe_allow_html=True)
st.session_state["method"] = st.radio(
    "Choose correlation method",
    options=[
        "pearson",
        "spearman",
        "sparcc",
        "distance_corr",
        "joint_rpca",
    ],
    format_func=lambda x: {
        "pearson": "Pearson (linear relationship)",
        "spearman": "Spearman (monotonic relationship)",
        "sparcc": "SparCC (compositional correlation)",
        "distance_corr": "Distance correlation (non-linear relationships)",
        "joint_rpca": "Joint-RPCA (loading-space associations)",
    }[x],
    help="""
    Pearson: Measures linear correlation. Assumes normal distribution..

    Spearman: Rank-based correlation for monotonic relationships.
    
    SparCC: Designed for compositional microbiome data.

    Distance correlation: Detects non-linear dependencies between variables.

    Joint-RPCA: Exploratory loading-space associations from joint multi-omics ordination.

    """,
    horizontal=False,
    on_change=clear_processing,
)

allowed_methods = ["pearson", "spearman", "sparcc", "distance_corr", "joint_rpca"]
method = st.session_state["method"]
sparcc_blocked = False
sparcc_size = None

if "last_method" not in st.session_state:
    st.session_state["last_method"] = method

if st.session_state["last_method"] != method:
    st.session_state["target_results"] = None
    st.session_state["decoy_results"] = None
    st.session_state["Target_scores"] = None
    st.session_state["Decoy_scores"] = None
    st.session_state["last_method"] = method

if method == "sparcc":
    st.warning(
        "SparCC is intended for **compositional datasets**, such as microbiome count or relative-abundance tables. "
        "It estimates correlations while accounting for the constant-sum/compositional structure of the data. "
        "It is not recommended for general metabolomics or proteomics intensity tables unless those data have been "
        "processed into an appropriate compositional representation."
    )

    sparcc_iter = st.number_input(
        "SparCC iterations",
        min_value=1,
        max_value=5,
        value=1,
        step=1,
        key="sparcc_iter",
        help="Higher iterations repeat the SparCC estimation and can improve stability, but runtime increases with each iteration."
    )

    if "metabolome_ft" in st.session_state and "target_omics" in st.session_state:
        sparcc_size = assess_sparcc_run_size(
            st.session_state["metabolome_ft"],
            st.session_state["target_omics"],
            st.session_state.is_restricted_mode,
            sparcc_iter,
        )

        if sparcc_size["blocked"]:
            sparcc_blocked = True
            st.error(
                "SparCC is disabled for this dataset in the web app because the combined feature count is too large.\n\n"
                f"- SparCC memory risk: **{sparcc_size['risk']}**\n"
                f"- Combined features: **{sparcc_size['total_features']:,}**\n"
                f"- Iterations selected: **{sparcc_iter}**\n"
                f"- Web app limit: **{sparcc_size['restricted_limit']:,} combined features**\n"
                f"- Available memory in this environment: **{sparcc_size['available_memory_gb']:.2f} GB**\n\n"
                "SparCC can use many times more memory than the base matrix size during processing. "
                "Please reduce the number of features or run the app locally."
            )
        elif sparcc_size["warn"] or sparcc_size["risk"] != "Low":
            st.warning(
                "This SparCC run may use substantial memory and could fail on machines with limited RAM.\n\n"
                f"- SparCC memory risk: **{sparcc_size['risk']}**\n"
                f"- Combined features: **{sparcc_size['total_features']:,}**\n"
                f"- Iterations selected: **{sparcc_iter}**\n"
                f"- Local warning threshold: **{sparcc_size['local_warning_limit']:,} combined features**\n"
                f"- Available memory in this environment: **{sparcc_size['available_memory_gb']:.2f} GB**\n\n"
                "SparCC can use many times more memory than the base matrix size during processing. "
                "Consider reducing iterations or filtering features before running SparCC."
            )

    if sparcc_size is not None:
        st.info(
            f"SparCC iterations selected: **{sparcc_iter}**. "
            "More iterations can make the SparCC estimate more stable, but they also make the run slower. "
            "If the app becomes slow or fails, try reducing the number of iterations or filtering features first."
        )

if method == "joint_rpca":
    st.warning(
        "Joint-RPCA is an exploratory multi-omics ordination method. The scores reflect feature concordance "
        "in the Joint-RPCA loading space, not raw abundance correlations. FDR filtering is not currently available "
        "for Joint-RPCA results, so a decoy result table is not generated."
    )
    joint_rpca_col1, joint_rpca_col2 = st.columns(2)
    with joint_rpca_col1:
        joint_rpca_min_feature_frequency = st.number_input(
            "Minimum feature frequency",
            min_value=1,
            value=2,
            step=1,
            key="joint_rpca_min_feature_frequency",
            help="Features with total frequency below this value are filtered before Joint-RPCA."
        )
    with joint_rpca_col2:
        joint_rpca_iter = st.number_input(
            "Joint-RPCA iterations",
            min_value=1,
            max_value=10,
            value=5,
            step=1,
            key="joint_rpca_iter",
            help="Higher iterations may improve factorization stability but increase runtime."
        )

if method not in allowed_methods:
    st.info(f"⚠️ **{method}** correlation is not implemented yet. ")

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
                                 disabled= st.session_state.get('no_correlations', 0) >= st.session_state.max_corr
                                 or method not in allowed_methods
                                 or sparcc_blocked)

# Perform Correlation Analysis Logic
if all(key in st.session_state for key in ['target_dataframe', 'decoy_dataframe', 'metabolome_ft', 'binned_omics_table', 'target_omics', 'decoy_omics']):
    
    # Estimate total number of correlations
    correlations = st.session_state['no_correlations']

    if run_corr_clicked and method in ["pearson", "spearman"]:
        st.session_state["processing"] = True
        st.session_state["processing_method"] = st.session_state["method"]

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

            method_label = "Pearson" if method == "pearson" else "Spearman"
            target_progress_cb = make_done_total_progress_callback(
                status_text,
                progress_bar,
                0.0,
                0.5,
                f"Running {method_label} correlation on target",
            )
            decoy_progress_cb = make_done_total_progress_callback(
                status_text,
                progress_bar,
                0.5,
                1.0,
                f"Running {method_label} correlation on decoy",
            )

            with st.spinner("Calculating correlations..."):
                st.session_state['target_results'] = calculate_correlations_parallel(
                    st.session_state['target_dataframe'], 
                    st.session_state['metabolome_ft'], 
                    st.session_state['target_omics'],
                    num_workers=4,
                    progress_callback=target_progress_cb,
                )
                progress_bar.progress(0.5)
                status_text.info(f"Target {method_label} correlation completed. Starting decoy...")

                st.session_state['decoy_results'] = calculate_correlations_parallel(
                    st.session_state['decoy_dataframe'], 
                    st.session_state['metabolome_ft'], 
                    st.session_state['decoy_omics'],
                    num_workers=4,
                    progress_callback=decoy_progress_cb,
                )
                update_method_progress(
                    status_text,
                    progress_bar,
                    f"Running {method_label} correlation on decoy",
                    1.0,
                    0,
                    0.5,
                    1.0,
                )

            st.session_state["processing"] = False  # Reset flag when done
            st.success("✅ Correlation analysis completed!")

    elif run_corr_clicked and method == "sparcc":
        st.session_state["processing"] = True
        st.session_state["processing_method"] = st.session_state["method"]

        try:
            metabolome_ft = st.session_state["metabolome_ft"].copy()
            target_omics = st.session_state["target_omics"].copy()
            iter_num = st.session_state.get("sparcc_iter", 1)

            sparcc_size = assess_sparcc_run_size(
                metabolome_ft,
                target_omics,
                st.session_state.is_restricted_mode,
                iter_num,
            )
            if sparcc_size["blocked"]:
                st.session_state["processing"] = False
                st.error(
                    "SparCC is disabled for this dataset in the web app because the combined feature count is too large. "
                    f"The web app limit is {sparcc_size['restricted_limit']:,} combined features. "
                    "Please reduce the number of features or run the app locally."
                )
                st.stop()

            # --- UI placeholders for live progress ---
            status_text = st.empty()
            progress_bar = st.progress(0.0)
    
            status_text.info("Preparing SparCC analysis...")
    
            def estimate_sparcc_seconds(n_features, iter_num):
                return max(10, (n_features / 1000) ** 2 * iter_num * 3)
            
            target_n_features = metabolome_ft.shape[0] + target_omics.shape[0]
            target_estimated_seconds = estimate_sparcc_seconds(target_n_features, iter_num)

            # --- TARGET ---
            target_progress_cb = make_estimated_progress_callback(
                status_text,
                progress_bar,
                0.0,
                1.0,
                "Running SparCC on target",
            )
            corr_target, target_results, common_samples_target = run_sparcc_for_omics_pair(
                metabolome_ft=metabolome_ft,
                omics_ft=target_omics,
                sparcc_script_path=get_sparcc_script_path(),
                iter_num=iter_num,
                progress_callback=target_progress_cb,
                estimated_seconds=target_estimated_seconds,
            )

            update_method_progress(
                status_text,
                progress_bar,
                "Running SparCC on target",
                1.0,
                0,
                0.0,
                1.0,
            )
            
            st.session_state["corr_target"] = corr_target
            st.session_state["target_results"] = target_results
            st.session_state["corr_decoy"] = None
            st.session_state["decoy_results"] = {}
                                    
            st.success("✅ SparCC analysis completed! (P-values/R2 are placeholders)")
            
        except Exception as e:
            st.error(f"❌ SparCC failed: {e}")
            st.session_state["processing"] = False
            st.stop()

        st.session_state["processing"] = False
    
    elif run_corr_clicked and method == "distance_corr":
        st.session_state["processing"] = True
        st.session_state["processing_method"] = st.session_state["method"]

        try:
            with st.spinner("Running distance correlation..."):
                metabolome_ft = st.session_state["metabolome_ft"].copy()
                target_omics = st.session_state["target_omics"].copy()
                decoy_omics = st.session_state["decoy_omics"].copy()

                # --- UI placeholders for live progress ---
                status_text = st.empty()
                progress_bar = st.progress(0.0)

                status_text.info("Preparing distance correlation analysis...")

                # --- TARGET ---
                target_progress_cb = make_timed_done_total_progress_callback(
                    status_text,
                    progress_bar,
                    0.0,
                    0.5,
                    "Running distance correlation on target",
                )
                target_results, common_samples_target = run_distance_correlation(
                    metabolome_ft=metabolome_ft,
                    omics_ft=target_omics,
                    progress_callback=target_progress_cb,
                )

                progress_bar.progress(0.5)
                status_text.info("Target distance correlation completed. Starting decoy...")

                # --- DECOY ---
                decoy_progress_cb = make_timed_done_total_progress_callback(
                    status_text,
                    progress_bar,
                    0.5,
                    1.0,
                    "Running distance correlation on decoy",
                )
                decoy_results, common_samples_decoy = run_distance_correlation(
                    metabolome_ft=metabolome_ft,
                    omics_ft=decoy_omics,
                    progress_callback=decoy_progress_cb,
                )

                progress_bar.progress(1.0)
                update_method_progress(
                    status_text,
                    progress_bar,
                    "Running distance correlation on decoy",
                    1.0,
                    0,
                    0.5,
                    1.0,
                )

                st.success("✅ Distance correlation analysis completed!")
                
                st.expander("ℹ️ Distance correlation: interpretation & notes", expanded=False).markdown(
                    """
                #### 📊 What does distance correlation measure?
                - The **Estimate** column represents the strength of association between features  
                - Values range from **0 (no association)** to **1 (strong association)**  
                - Captures **linear and nonlinear relationships**

                #### 💡 When should you use distance correlation?
                - To detect **non-linear relationships**
                - When Pearson/Spearman may miss dependencies
                - For **general association discovery**, not directional interpretation

                #### ⚠️ No direction (positive / negative)
                - All values are **non-negative**
                - Distance correlation does **not distinguish** whether a relationship is increasing or decreasing  
                - A high score indicates **strong dependence**, but not direction

                #### 🧪 About P-value and R2
                - **P-value** and **R2** are included as placeholders for consistency across methods  
                - Distance correlation does **not inherently produce these statistics**. These values are set to **NaN**
                - Interpretation should be based on the **Estimate** column

                #### 🎯 Interpretation with FDR (target–decoy)
                - FDR filtering retains **high-confidence associations**
                - Only higher values appear after filtering, but they still represent **strength of association**, not direction
                """
                )

                st.session_state["target_results"] = target_results
                st.session_state["decoy_results"] = decoy_results
                st.session_state["common_samples_target"] = common_samples_target
                st.session_state["common_samples_decoy"] = common_samples_decoy
                st.session_state["processing"] = False

        except Exception as e:
            st.session_state["processing"] = False
            st.error(f"Distance correlation failed: {e}")

    elif run_corr_clicked and method == "joint_rpca":
        st.session_state["processing"] = True
        st.session_state["processing_method"] = st.session_state["method"]

        try:
            status_text = st.empty()
            progress_bar = st.progress(0.0)
            update_method_progress(
                status_text,
                progress_bar,
                "Running Joint-RPCA on target",
                0.05,
                0,
            )

            with st.spinner("Running Joint-RPCA..."):
                joint_rpca_result = run_joint_rpca_for_omics_pair(
                    metabolome_ft=st.session_state["metabolome_ft"].copy(),
                    omics_ft=st.session_state["target_omics"].copy(),
                    max_iterations=st.session_state.get("joint_rpca_iter", 5),
                    min_feature_frequency=st.session_state.get("joint_rpca_min_feature_frequency", 2),
                )

            update_method_progress(
                status_text,
                progress_bar,
                "Running Joint-RPCA on target",
                1.0,
                0,
            )

            st.session_state["target_results"] = joint_rpca_result["scores"]
            st.session_state["decoy_results"] = pd.DataFrame()
            st.session_state["joint_rpca_cross_block"] = joint_rpca_result["cross_block"]
            st.session_state["joint_rpca_all_feature_scores"] = joint_rpca_result["all_feature_scores"]
            st.session_state["common_samples_target"] = joint_rpca_result["common_samples"]
            st.session_state["processing"] = False

            st.success("✅ Joint-RPCA completed! Scores are loading-space associations, not raw abundance correlations.")

        except Exception as e:
            st.session_state["processing"] = False
            st.error(f"Joint-RPCA failed: {e}")

else:
    st.warning("⚠️ Please input the data first to continue the analysis.")

# Show a processing message if the analysis is running

if (
    st.session_state.get("processing", False)
    and st.session_state.get("processing_method") == st.session_state.get("method")
):
    st.warning("⏳ Correlation analysis running… Please wait.")


# Check if correlation results are available before processing
melted_target = None
melted_decoy = None

if ("target_results" in st.session_state
    and "decoy_results" in st.session_state 
    and st.session_state["target_results"] is not None
    and st.session_state["decoy_results"] is not None
):
    if method in ["pearson", "spearman", "sparcc"]:
        melted_target = melt_correlation_results(st.session_state["target_results"], method=method)
        melted_decoy = melt_correlation_results(st.session_state["decoy_results"], method=method)

    elif method in ["distance_corr", "joint_rpca"]:
        if isinstance(st.session_state["target_results"], pd.DataFrame) and isinstance(st.session_state["decoy_results"], pd.DataFrame):
            melted_target = st.session_state["target_results"].copy()
            melted_decoy = st.session_state["decoy_results"].copy()
    
    else:
        st.warning(f"No melting/display handler defined yet for method: {method}")
        melted_target = None
        melted_decoy = None
   
    if melted_target is not None and melted_decoy is not None:
        st.session_state["Target_scores"] = melted_target
        st.session_state["Decoy_scores"] = melted_decoy

        df_target_size = estimate_df_size(melted_target)
        df_decoy_size = estimate_df_size(melted_decoy)
        has_decoy_results = isinstance(melted_decoy, pd.DataFrame) and not melted_decoy.empty
        has_decoy_downloads = (
            method not in ["sparcc", "joint_rpca"]
            and has_graphml_edge_columns(melted_decoy)
        )

        target_csv_data = convert_df_to_csv(melted_target)
        decoy_csv_data = convert_df_to_csv(melted_decoy) if has_decoy_downloads else None

    # Display or provide download option based on size
        if df_target_size <= 200:
            with st.expander(f"Correlation Scores of Target Dataframe {melted_target.shape}"):
                st.dataframe(melted_target)
                target_download_col1, target_download_col2 = st.columns(2)
                with target_download_col1:
                    st.download_button(
                        label="Download Target DataFrame as CSV",
                        data=target_csv_data,
                        file_name="target_correlation_results.csv",
                        mime="text/csv"
                    )
                with target_download_col2:
                    render_correlation_graphml_download(
                        st.session_state.get("node_table"),
                        melted_target,
                        method,
                        "target_correlation_network.graphml",
                        "target_graphml_unfiltered",
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
                render_correlation_graphml_download(
                    st.session_state.get("node_table"),
                    melted_target,
                    method,
                    "target_correlation_network.graphml",
                    "target_graphml_unfiltered_large",
                )

        if has_decoy_results and df_decoy_size <= 200:
            with st.expander(f"Correlation Scores of Decoy Dataframe {melted_decoy.shape}"):
                st.dataframe(melted_decoy)
                if has_decoy_downloads:
                    decoy_download_col1, decoy_download_col2 = st.columns(2)
                    with decoy_download_col1:
                        st.download_button(
                            label="Download Decoy Correlation Results (CSV)",
                            data=decoy_csv_data,
                            file_name="decoy_correlation_results.csv",
                            mime="text/csv"
                        )
                    with decoy_download_col2:
                        render_correlation_graphml_download(
                            st.session_state.get("node_table"),
                            melted_decoy,
                            method,
                            "decoy_correlation_network.graphml",
                            "decoy_graphml_unfiltered",
                        )
        elif has_decoy_results:
            if 'corr_c1' not in locals():  # use same row if available
                corr_c1, corr_c2 = st.columns(2)
            with corr_c2:
                if has_decoy_downloads:
                    st.download_button(
                        label="Download Decoy Correlation Results (CSV)",
                        data=decoy_csv_data,
                        file_name="decoy_correlation_results.csv",
                        mime="text/csv"
                    )
                    render_correlation_graphml_download(
                        st.session_state.get("node_table"),
                        melted_decoy,
                        method,
                        "decoy_correlation_network.graphml",
                        "decoy_graphml_unfiltered_large",
                    )

# Add space before FDR analysis
st.divider()
false_discovery_rate_section(method)

if (
    isinstance(st.session_state.get("Target_scores"), pd.DataFrame)
    and not st.session_state["Target_scores"].empty
):
    st.divider()
    association_heatmap_section(st.session_state["Target_scores"], method)
