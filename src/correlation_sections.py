import streamlit as st
from streamlit_extras.stylable_container import stylable_container

from src.correlation import display_correlation_message
from src.fdr import calculate_fdr, plot_fdr_figures
from src.molnet import generate_graphml_download_data
from src.visualization import (
    CORROMICS_EXAMPLE_MICROBES,
    CORROMICS_EXAMPLE_STANDARD_IDS,
    build_association_heatmap,
    get_split_label_options,
    parse_feature_id_text,
)


@st.fragment
def association_heatmap_section(target_scores, method):
    st.subheader("Association Heatmap")
    heatmap_modes = ["Top associations", "Selected feature IDs"]
    if st.session_state.get("input_dataset") == "corromics_manuscript_example":
        heatmap_modes.append("CorrOmics example standards")

    heatmap_mode = st.radio(
        "Heatmap mode",
        heatmap_modes,
        horizontal=True,
        key="heatmap_mode",
    )

    selected_heatmap_features = []
    selected_heatmap_variables = None
    heatmap_top_n = 50
    heatmap_max_rows = 30
    heatmap_label_mode = "Split by separator"
    heatmap_label_separator = "_"
    heatmap_label_part_index = 0

    heatmap_col1, heatmap_col2 = st.columns(2)
    with heatmap_col1:
        if heatmap_mode == "Top associations":
            heatmap_top_n = st.select_slider(
                "Number of top associations",
                options=[50, 100, 250, 500],
                value=100,
                key="heatmap_top_n",
                help="Associations are ranked by absolute Estimate.",
            )
            heatmap_max_rows = st.slider(
                "Maximum feature rows",
                min_value=5,
                max_value=50,
                value=30,
                step=5,
                key="heatmap_max_rows",
            )
        elif heatmap_mode == "Selected feature IDs":
            feature_id_text = st.text_area(
                "Feature IDs",
                key="heatmap_feature_ids",
                placeholder="Paste up to 30 feature IDs, one per line or separated by commas",
                height=120,
            )
            selected_heatmap_features = parse_feature_id_text(feature_id_text, max_ids=30)
            entered_ids = [
                feature_id.strip()
                for feature_id in feature_id_text.replace(",", "\n").splitlines()
                if feature_id.strip()
            ] if feature_id_text else []
            if len(entered_ids) > 30:
                st.warning("Only the first 30 feature IDs will be used for the heatmap.")
        else:
            selected_heatmap_features = CORROMICS_EXAMPLE_STANDARD_IDS
            selected_heatmap_variables = CORROMICS_EXAMPLE_MICROBES
            st.info(
                "Using the CorrOmics example standard metabolites and microbes in a fixed order."
            )

    with heatmap_col2:
        if heatmap_mode != "CorrOmics example standards":
            heatmap_max_cols = st.slider(
                "Maximum variable columns",
                min_value=5,
                max_value=50,
                value=30,
                step=5,
                key="heatmap_max_cols",
                help="Variables are ranked by their strongest absolute Estimate among the plotted features.",
            )
        heatmap_label_mode = st.radio(
            "Feature row labels",
            ["Split by separator", "Original index"],
            horizontal=True,
            key="heatmap_label_mode",
            help="Choose whether the heatmap y-axis shows the original feature index or one part of a separated feature name.",
        )
        if heatmap_label_mode == "Split by separator":
            heatmap_label_separator = st.text_input(
                "Separator in feature index",
                value="_",
                max_chars=10,
                key="heatmap_label_separator",
                help="For an index like 1231_gzefdw_453, use _ to create name1=1231, name2=gzefdw, name3=453.",
            )
            feature_label_source = (
                target_scores["Feature"].dropna().astype(str).drop_duplicates().tolist()
                if "Feature" in target_scores.columns else []
            )
            heatmap_label_options = get_split_label_options(
                feature_label_source,
                heatmap_label_separator,
            )
            if not heatmap_label_options:
                heatmap_label_options = ["name1"]
            heatmap_label_part = st.selectbox(
                "Show split field",
                heatmap_label_options,
                index=0,
                key="heatmap_label_part",
                help="Default is name1, the first value before the separator.",
            )
            heatmap_label_part_index = heatmap_label_options.index(heatmap_label_part)

    heatmap_fig = build_association_heatmap(
        target_scores,
        method,
        mode=heatmap_mode,
        selected_features=selected_heatmap_features,
        selected_variables=selected_heatmap_variables,
        top_n=heatmap_top_n,
        max_rows=heatmap_max_rows,
        max_cols=heatmap_max_cols,
        feature_label_mode=heatmap_label_mode,
        feature_label_separator=heatmap_label_separator,
        feature_label_part_index=heatmap_label_part_index,
    )

    if heatmap_fig is not None:
        st.plotly_chart(heatmap_fig, use_container_width=True)
    elif heatmap_mode != "Top associations" and selected_heatmap_features:
        st.info("No matching Feature IDs were found in the current target scores.")
    elif heatmap_mode == "Selected feature IDs":
        st.info("Enter feature IDs to generate a selected-feature heatmap.")
    else:
        st.info("No numeric Estimate values are available for the heatmap.")


@st.fragment
def false_discovery_rate_section(method):
    st.markdown("## False Discovery Rate")

    if method == "sparcc":
        st.warning(
            "SparCC produces correlation estimates for compositional data, but this app does not currently compute calibrated SparCC **p-values**. "
            "The target-decoy FDR approach used for Pearson/Spearman has not yet been validated for SparCC results. "
            "Because FDR filtering is unavailable for SparCC here, a decoy result table is not generated. "
            "Please interpret SparCC estimates as exploratory association scores."
        )
    elif method == "joint_rpca":
        st.warning(
            "FDR filtering is currently disabled for Joint-RPCA. Joint-RPCA scores are loading-space associations from "
            "a joint ordination, not raw abundance correlations with calibrated p-values. Please interpret them as "
            "exploratory feature-concordance scores. A decoy result table is not generated for Joint-RPCA."
        )

    with stylable_container(
        key="run_fdr_container",
        css_styles="""
            button { background-color: #00689E; color: white; font-weight: bold; width: 100%; }
            button:hover { background-color: #003043; transform: scale(1.05); }
            button:disabled { background-color: #A9A9A9;  cursor: not-allowed;}
        """,
    ):
        run_fdr_button_clicked = st.button(
            "Apply FDR",
            key="run_fdr",
            disabled=st.session_state.get("no_correlations", 0) >= st.session_state.max_corr
            or method in ["sparcc", "joint_rpca"],
        )

    if "run_fdr_clicked" not in st.session_state:
        st.session_state["run_fdr_clicked"] = False
    if "neg_fdr_cutoff" not in st.session_state:
        st.session_state["neg_fdr_cutoff"] = -0.5
    if "pos_fdr_cutoff" not in st.session_state:
        st.session_state["pos_fdr_cutoff"] = 0.5
    if "filtered_target_csv" not in st.session_state:
        st.session_state["filtered_target_csv"] = None

    if (
        "Target_scores" in st.session_state
        and "Decoy_scores" in st.session_state
        and st.session_state["Target_scores"] is not None
        and st.session_state["Decoy_scores"] is not None
    ):
        if run_fdr_button_clicked:
            correlations = st.session_state["no_correlations"]
            if correlations >= st.session_state.max_corr and not st.session_state.is_restricted_mode:
                st.error(f"❌ Too many correlations to compute in this environment (≥ {st.session_state.max_corr}). Please clone or download the app and run it locally. This helps avoid memory crashes in the cloud environment.")
                st.stop()

            st.session_state["run_fdr_clicked"] = True

            target_scores = st.session_state["Target_scores"]
            decoy_scores = st.session_state["Decoy_scores"]

            overall_fdr_table = calculate_fdr(
                target_scores,
                decoy_scores,
                score_range=(-1, 1),
                bin_size=0.001,
            )
            fig_fdr, fig_histogram = plot_fdr_figures(
                overall_fdr_table,
                target_scores,
                decoy_scores,
            )

            st.session_state["fig_histogram"] = fig_histogram
            st.session_state["fig_fdr"] = fig_fdr
    else:
        st.warning("Please run correlation analysis to continue this step.")

    if (
        st.session_state.get("run_fdr_clicked", False)
        and (
            st.session_state.get("no_correlations", 0) < st.session_state.max_corr
            or not st.session_state.is_restricted_mode
        )
    ):
        st.plotly_chart(st.session_state["fig_histogram"])
        st.plotly_chart(st.session_state["fig_fdr"])

        melted_target = st.session_state["Target_scores"]

        if method == "distance_corr":
            st.info(
                "Distance correlation scores are non-negative. Use the positive cutoff to keep stronger associations; "
                "the negative cutoff is not relevant for this method."
            )
        else:
            st.write("Select the positive and negative cutoffs for the correlation scores based on your FDR-curve")

        input_col1, input_col2 = st.columns(2)

        with input_col1:
            st.session_state["neg_fdr_cutoff"] = st.number_input(
                "Enter Negative FDR Cutoff (e.g., -0.5):",
                value=st.session_state.get("neg_fdr_cutoff", -0.5),
                min_value=-1.0,
                max_value=0.0,
                step=0.05,
            )

        with input_col2:
            st.session_state["pos_fdr_cutoff"] = st.number_input(
                "Enter Positive FDR Cutoff (e.g., 0.5):",
                value=st.session_state.get("pos_fdr_cutoff", 0.5),
                min_value=0.0,
                max_value=1.0,
                step=0.05,
            )

        filtered_target = melted_target[
            (melted_target["Estimate"] <= st.session_state["neg_fdr_cutoff"])
            | (melted_target["Estimate"] >= st.session_state["pos_fdr_cutoff"])
        ]

        st.session_state["filtered_target"] = filtered_target
        st.session_state["filtered_target_csv"] = filtered_target.to_csv(index=False).encode("utf-8")

        if st.session_state["filtered_target"] is not None:
            st.markdown("#### Filtered Correlation Scores Based on FDR Cutoffs")

            display_correlation_message(melted_target.shape[0], filtered_target.shape[0])
            st.write("")
            if st.session_state["filtered_target_csv"] is not None:
                _, csv_col, graphml_col, _ = st.columns([1, 1, 1, 1])

                with csv_col:
                    st.download_button(
                        label="Download Results as CSV",
                        data=st.session_state["filtered_target_csv"],
                        file_name="filtered_correlations_results.csv",
                        mime="text/csv",
                    )

                with graphml_col:
                    if "node_table" in st.session_state:
                        graphml_data = generate_graphml_download_data(
                            st.session_state["node_table"],
                            st.session_state["filtered_target"],
                        )
                        if graphml_data:
                            st.download_button(
                                "Download GraphML",
                                graphml_data,
                                file_name="correlation_network.graphml",
                                mime="application/graphml+xml",
                            )
                        else:
                            st.error("❌ Failed to generate GraphML.")
                    else:
                        st.error("❌ Node table not found in session state. Please run the correlation analysis first.")

    elif (
        st.session_state.get("run_fdr_clicked", False)
        and st.session_state.get("no_correlations", 0) >= st.session_state.max_corr
        and st.session_state.is_restricted_mode
    ):
        st.error(f"❌ As correlations exceed {st.session_state.max_corr}, no correlations were computed for this level. 💡 Please clone or download the app and run it locally. This helps avoid memory crashes in the cloud environment.")
        st.stop()

    with stylable_container(
        key="fdr_reset_container",
        css_styles="""
            button {background-color: #FF000066; color: black; font-weight: bold; width: 100%;}
            button:hover { background-color: #003043; transform: scale(1.05); }
            button:disabled { background-color: #A9A9A9;  cursor: not-allowed;}
        """,
    ):
        reset_fdr_clicked = st.button(
            "Reset FDR Analysis",
            key="reset_fdr",
        )

    if st.session_state["run_fdr_clicked"] and reset_fdr_clicked:
        st.session_state["run_fdr_clicked"] = False
        st.session_state["neg_fdr_cutoff"] = -0.5
        st.session_state["pos_fdr_cutoff"] = 0.5
        st.session_state["filtered_target"] = None
        st.session_state["filtered_target_csv"] = None
        st.rerun(scope="fragment")
