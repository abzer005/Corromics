import streamlit as st
from .common import *  # Importing common functionalities from the 'common' module
import pandas as pd
import numpy as np
import networkx as nx

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
    
    # Add nodes with attributes
    for node_id, node_data in node_table.iterrows():
        
        # Add node with shape attribute
        G.add_node(node_id, **node_data.to_dict())


    # Clean the Feature column to keep only numeric parts
    edge_table['Feature'] = edge_table['Feature'].astype(str).str.extract(r'(^\d+)', expand=False)

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

        # Edge attributes
        edge_attributes = {
            'Correlation_Score': row['Estimate'],
            'Sign_Score': row['Sign_Estimate'],
            'P_Value': row['P-value'],
            'BH_Corrected_P_Value': row['BH-Corrected P-Value'],
            'R2': row['R2'],
            'Absolute_Correlation_Score': abs(row['Estimate']),
        }

        # Add the edge
        G.add_edge(clusterid1, clusterid2, key="primary", **edge_attributes)

    # ✅ Write the graph to GraphML
    nx.write_graphml(G, output_file)

    return output_file


