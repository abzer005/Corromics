import pandas as pd
import plotly.graph_objects as go


CORROMICS_EXAMPLE_STANDARD_IDS = ["1873", "197", "1086", "256", "984", "31", "730", "3621"]
CORROMICS_EXAMPLE_MICROBES = [
    "Flavobacterium pectinovorum",
    "Frigoribacterium faeni",
    "Pseudomonas koreensis",
    "Microbacterium proteolyticum",
    "Bacillus altitudinis",
    "Sphingomonas faeni",
    "Nocardioides cavernae",
    "Methylobacterium goesingense",
]


def prepare_top_association_heatmap(scores, top_n=50, max_rows=30, max_cols=30):
    required_columns = {"Feature", "Variable", "Estimate"}
    if scores is None or not required_columns.issubset(scores.columns):
        return pd.DataFrame()

    plot_df = scores[["Feature", "Variable", "Estimate"]].copy()
    plot_df["Estimate"] = pd.to_numeric(plot_df["Estimate"], errors="coerce")
    plot_df = plot_df.dropna(subset=["Feature", "Variable", "Estimate"])

    if plot_df.empty:
        return pd.DataFrame()

    plot_df["Abs_Estimate"] = plot_df["Estimate"].abs()
    top_pairs = plot_df.nlargest(top_n, "Abs_Estimate")

    row_order = (
        top_pairs.groupby("Feature")["Abs_Estimate"]
        .max()
        .sort_values(ascending=False)
        .head(max_rows)
        .index
    )
    col_order = (
        top_pairs.groupby("Variable")["Abs_Estimate"]
        .max()
        .sort_values(ascending=False)
        .head(max_cols)
        .index
    )

    heatmap_df = (
        plot_df[
            plot_df["Feature"].isin(row_order)
            & plot_df["Variable"].isin(col_order)
        ]
        .pivot_table(
            index="Feature",
            columns="Variable",
            values="Estimate",
            aggfunc="mean",
        )
        .reindex(index=row_order, columns=col_order)
    )

    return heatmap_df


def parse_feature_id_text(feature_id_text, max_ids=30):
    if not feature_id_text:
        return []

    raw_ids = feature_id_text.replace(",", "\n").splitlines()
    feature_ids = []
    for feature_id in raw_ids:
        feature_id = feature_id.strip()
        if feature_id and feature_id not in feature_ids:
            feature_ids.append(feature_id)

    return feature_ids[:max_ids]


def _feature_id_prefix(feature_value):
    return str(feature_value).split("_", 1)[0]


def get_split_label_options(feature_values, separator):
    if not separator:
        return []

    max_parts = 0
    for feature in feature_values:
        max_parts = max(max_parts, len(str(feature).split(separator)))

    return [f"name{i}" for i in range(1, max_parts + 1)]


def _split_feature_label(feature_value, separator, part_index):
    feature_text = str(feature_value)
    if not separator:
        return feature_text

    parts = feature_text.split(separator)
    if 0 <= part_index < len(parts) and parts[part_index]:
        return parts[part_index]
    return feature_text


def _normalize_variable_name(variable_value):
    return str(variable_value).replace("_", " ")


def _shorten_feature_labels(
    feature_values,
    label_mode="Split by separator",
    separator="_",
    split_part_index=0,
):
    labels = []
    seen = {}
    for feature in feature_values:
        if label_mode == "Original index":
            short_label = str(feature)
        else:
            short_label = _split_feature_label(feature, separator, split_part_index)
        seen[short_label] = seen.get(short_label, 0) + 1
        if seen[short_label] > 1:
            short_label = f"{short_label} ({seen[short_label]})"
        labels.append(short_label)
    return labels


def prepare_selected_feature_heatmap(
    scores,
    selected_features,
    selected_variables=None,
    max_variables=50,
):
    required_columns = {"Feature", "Variable", "Estimate"}
    if scores is None or not required_columns.issubset(scores.columns) or not selected_features:
        return pd.DataFrame()

    plot_df = scores[["Feature", "Variable", "Estimate"]].copy()
    plot_df["Feature"] = plot_df["Feature"].astype(str)
    plot_df["Variable"] = plot_df["Variable"].astype(str)
    plot_df["Variable_Display"] = plot_df["Variable"].map(_normalize_variable_name)
    plot_df["Feature_ID"] = plot_df["Feature"].map(_feature_id_prefix)
    plot_df["Estimate"] = pd.to_numeric(plot_df["Estimate"], errors="coerce")
    plot_df = plot_df.dropna(subset=["Feature", "Variable", "Estimate"])
    selected_features = [str(feature) for feature in selected_features]
    plot_df = plot_df[plot_df["Feature_ID"].isin(selected_features)]

    if selected_variables:
        selected_variables = [str(variable) for variable in selected_variables]
        plot_df = plot_df[plot_df["Variable_Display"].isin(selected_variables)]

    if plot_df.empty:
        return pd.DataFrame()

    plot_df["Abs_Estimate"] = plot_df["Estimate"].abs()
    selected_order = []
    for feature_id in selected_features:
        matching_features = plot_df.loc[plot_df["Feature_ID"] == feature_id, "Feature"].drop_duplicates()
        selected_order.extend(matching_features.tolist())

    if selected_variables:
        variable_order = [
            variable for variable in selected_variables
            if variable in set(plot_df["Variable_Display"])
        ]
    else:
        variable_order = (
            plot_df.groupby("Variable_Display")["Abs_Estimate"]
            .max()
            .sort_values(ascending=False)
            .head(max_variables)
            .index
        )

    return (
        plot_df[plot_df["Variable_Display"].isin(variable_order)]
        .pivot_table(
            index="Feature",
            columns="Variable_Display",
            values="Estimate",
            aggfunc="mean",
        )
        .reindex(index=selected_order, columns=variable_order)
    )


def build_association_heatmap(
    scores,
    method,
    mode="Top associations",
    selected_features=None,
    selected_variables=None,
    top_n=50,
    max_rows=30,
    max_cols=30,
    feature_label_mode="Split by separator",
    feature_label_separator="_",
    feature_label_part_index=0,
):
    if mode in ["Selected feature IDs", "CorrOmics example standards"]:
        heatmap_df = prepare_selected_feature_heatmap(
            scores,
            selected_features=selected_features,
            selected_variables=selected_variables,
            max_variables=max_cols,
        )
    else:
        heatmap_df = prepare_top_association_heatmap(
            scores,
            top_n=top_n,
            max_rows=max_rows,
            max_cols=max_cols,
        )

    if heatmap_df.empty:
        return None

    signed_methods = {"pearson", "spearman", "sparcc", "joint_rpca"}
    z_values = heatmap_df.values
    full_feature_names = heatmap_df.index.astype(str).tolist()
    short_feature_labels = _shorten_feature_labels(
        full_feature_names,
        label_mode=feature_label_mode,
        separator=feature_label_separator,
        split_part_index=feature_label_part_index,
    )
    variable_labels = heatmap_df.columns.astype(str).tolist()
    max_abs = pd.Series(z_values.ravel()).dropna().abs().max()
    max_abs = 1 if pd.isna(max_abs) or max_abs == 0 else max_abs

    heatmap_kwargs = {
        "z": z_values,
        "x": variable_labels,
        "y": short_feature_labels,
        "customdata": [
            [feature_name for _ in variable_labels]
            for feature_name in full_feature_names
        ],
        "colorbar": {"title": "Estimate"},
        "hovertemplate": (
            "Feature: %{customdata}<br>"
            "Variable: %{x}<br>"
            "Estimate: %{z:.3f}<extra></extra>"
        ),
    }

    if method in signed_methods:
        heatmap_kwargs.update(
            colorscale="RdBu",
            zmin=-max_abs,
            zmax=max_abs,
            zmid=0,
        )
    else:
        heatmap_kwargs.update(colorscale="Viridis", zmin=0)

    fig = go.Figure(data=go.Heatmap(**heatmap_kwargs))
    if mode in ["Selected feature IDs", "CorrOmics example standards"]:
        text_values = heatmap_df.apply(
            lambda column: column.map(
                lambda value: "" if pd.isna(value) else f"{value:.2f}"
            )
        ).to_numpy()
        fig.add_trace(
            go.Scatter(
                x=[
                    variable_label
                    for _ in short_feature_labels
                    for variable_label in variable_labels
                ],
                y=[
                    feature_label
                    for feature_label in short_feature_labels
                    for _ in variable_labels
                ],
                mode="text",
                text=text_values.ravel(),
                textfont=dict(color="black", size=12),
                hoverinfo="skip",
                showlegend=False,
            )
        )
    fig.update_layout(
        height=min(max(450, 24 * len(heatmap_df.index)), 900),
        margin=dict(l=10, r=10, t=30, b=10),
        xaxis=dict(title="", tickangle=45, type="category"),
        yaxis=dict(title="", autorange="reversed", type="category"),
    )
    return fig
