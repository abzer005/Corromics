import streamlit as st
from .common import *  # Importing common functionalities from the 'common' module
import pandas as pd
import numpy as np
import networkx as nx
import tempfile
import os


def _graphml_safe_value(value):
    if pd.isna(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def _add_edge_attribute(attributes, output_name, row, source_name):
    if source_name not in row.index:
        return

    value = _graphml_safe_value(row[source_name])
    if value is not None:
        attributes[output_name] = value


def has_graphml_edge_columns(edge_table):
    required_columns = {"Feature", "Variable", "Estimate"}
    return (
        isinstance(edge_table, pd.DataFrame)
        and not edge_table.empty
        and required_columns.issubset(edge_table.columns)
    )


def generate_graphml_corromics(node_table, edge_table, output_file):
    """
    Generate a GraphML file from the node and edge tables.

    Parameters:
    node_table (pd.DataFrame): DataFrame containing node information.
    edge_table (pd.DataFrame): DataFrame containing edge (correlation) information.
    output_file (str): The name of the GraphML file to generate.
    """

    # Create a directed graph that allows multiple edges between nodes
    G = nx.MultiDiGraph()

    # Fetch node and edge tables from Streamlit session state if not provided
    if node_table is None:
        node_table = st.session_state.get('node_table')
    if edge_table is None:
        edge_table = st.session_state.get("filtered_target")

    if node_table is None or edge_table is None:
        st.error("Missing node or edge table data.")
        return None

    if not has_graphml_edge_columns(edge_table):
        return None
    
    # Add nodes with attributes
    for node_id, node_data in node_table.iterrows():
        
        # Add node with shape attribute
        G.add_node(node_id, **node_data.to_dict())


    edge_table = edge_table.copy()
    edge_table['Feature'] = edge_table['Feature'].astype(str)
    extracted_feature_ids = edge_table['Feature'].str.extract(r'(^\d+)', expand=False)
    edge_table['Feature'] = extracted_feature_ids.fillna(edge_table['Feature'])

    # Add the sign of the correlation estimate
    edge_table['Sign_Estimate'] = edge_table['Estimate'].apply(lambda x: 1 if x > 0 else (-1 if x < 0 else 0))

    # Add nodes with attributes from node_table
    for node_id, node_data in node_table.iterrows():
        G.add_node(node_id, **node_data.to_dict())

    # Add edges with attributes from edge_table
    for _, row in edge_table.iterrows():
        clusterid1 = row['Feature']
        clusterid2 = row['Variable']

        # Add nodes if missing
        if clusterid1 not in G.nodes:
            G.add_node(clusterid1, node_names=clusterid1)
        if clusterid2 not in G.nodes:
            G.add_node(clusterid2, node_names=clusterid2)

        edge_attributes = {
            'Correlation_Score': row['Estimate'],
            'Sign_Score': row['Sign_Estimate'],
            'Absolute_Correlation_Score': abs(row['Estimate']),
        }
        _add_edge_attribute(edge_attributes, 'P_Value', row, 'P-value')
        _add_edge_attribute(edge_attributes, 'BH_Corrected_P_Value', row, 'BH-Corrected P-Value')
        _add_edge_attribute(edge_attributes, 'R2', row, 'R2')
        _add_edge_attribute(edge_attributes, 'Method', row, 'Method')

        # Add the edge
        G.add_edge(clusterid1, clusterid2, key="primary", **edge_attributes)

    # ✅ Write the graph to GraphML
    nx.write_graphml(G, output_file)

    return output_file

def generate_graphml_download_data(node_table, edge_table):
    """
    Generate GraphML bytes without writing a persistent file into the project.
    """

    output_file = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".graphml", delete=False) as tmp:
            output_file = tmp.name

        graphml_path = generate_graphml_corromics(node_table, edge_table, output_file)
        if not graphml_path:
            return None

        with open(graphml_path, "rb") as graphml_file:
            return graphml_file.read()

    finally:
        if output_file and os.path.exists(output_file):
            os.remove(output_file)


def render_correlation_graphml_download(
    node_table,
    edge_table,
    method,
    file_name,
    key,
    allowed_methods=("pearson", "spearman", "sparcc", "distance_corr", "joint_rpca"),
):
    if method not in allowed_methods:
        return

    if node_table is None:
        st.info("GraphML download will be available after the node table is created.")
        return

    if not has_graphml_edge_columns(edge_table):
        return

    graphml_data = generate_graphml_download_data(node_table, edge_table)
    if graphml_data:
        st.download_button(
            "Download GraphML",
            graphml_data,
            file_name=file_name,
            mime="application/graphml+xml",
            key=key,
        )
    else:
        st.error("Failed to generate GraphML.")
