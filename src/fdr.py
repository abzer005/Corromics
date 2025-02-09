import streamlit as st
from .common import *  # Importing common functionalities from the 'common' module
from .correlation import *
import pandas as pd
import numpy as np
import plotly.graph_objects as go


def calculate_fdr(target, decoy, score_range=(-1, 1), bin_size=0.001):
    
    # Define bins based on specified range and bin size
    bins = np.arange(score_range[0], score_range[1] + bin_size, bin_size)

    # Calculate target and decoy counts per bin
    target_counts, _ = np.histogram(target['Estimate'], bins=bins)
    decoy_counts, _ = np.histogram(decoy['Estimate'], bins=bins)

    # Calculate FDR for each bin
    fdr_df = pd.DataFrame({'Range_min': bins[:-1], 
                           'Range_max': bins[1:], 
                           'Target_counts': target_counts,
                           'Decoy_counts': decoy_counts
                               })
    
    epsilon = 1e-10

    # Separate positive and negative bins
    positive_bins = fdr_df[fdr_df['Range_max'] >= 0].reset_index(drop=True)
    negative_bins = fdr_df[fdr_df['Range_max'] < 0].reset_index(drop=True)

    # Calculate cumulative FDR for positive bins (0 to 1, forward order)
    positive_bins = positive_bins[::-1].reset_index(drop=True)  # Reverse order for cumulative sum
    positive_bins['FDR'] = positive_bins['Decoy_counts'].cumsum() / (
        positive_bins['Target_counts'].cumsum()  + epsilon + positive_bins['Decoy_counts'].cumsum()
    )

    positive_bins = positive_bins[::-1].reset_index(drop=True)  # Restore original order

    # # Calculate cumulative FDR for negative bins (-1 to 0)
    negative_bins['FDR'] = negative_bins['Decoy_counts'].cumsum() / (
        negative_bins['Target_counts'].cumsum() +  epsilon + negative_bins['Decoy_counts'].cumsum()
    )
   
    # Combine positive and negative bins
    combined_fdr_df = pd.concat([negative_bins, positive_bins]).reset_index(drop=True)

    return combined_fdr_df


def plot_fdr_figures(combined_fdr_df, target, decoy):
    # Plot FDR Curve
    fig_fdr = go.Figure()
    fig_fdr.add_trace(go.Scatter(
        x=combined_fdr_df['Range_min'],
        y=combined_fdr_df['FDR'] * 100,
        mode='markers',
        name='FDR Score'
    ))

    # FDR Thresholds
    fdr_thresholds = [1, 5, 10]
    colors = ['green', 'orange', 'red']

    # Add threshold lines
    for fdr_value, color in zip(fdr_thresholds, colors):
        fig_fdr.add_shape(
            type="line",
            x0=combined_fdr_df['Range_min'].min(), 
            x1=combined_fdr_df['Range_min'].max(),
            y0=fdr_value, y1=fdr_value,
            line=dict(color=color, width=1, dash="dash"),
            name=f"{fdr_value}% FDR"
        )

        tolerance = 0.8  # Tolerance for threshold markers
        threshold_points = combined_fdr_df[
            (combined_fdr_df['FDR'] * 100 >= fdr_value - tolerance) &
            (combined_fdr_df['FDR'] * 100 <= fdr_value + tolerance)
        ]

        fig_fdr.add_trace(go.Scatter(
            x=threshold_points['Range_min'],
            y=threshold_points['FDR'] * 100,
            mode='markers',
            marker=dict(color=color, size=10, symbol='circle-open'),
            name=f"{fdr_value}% Threshold"
        ))

    fig_fdr.update_layout(
        title="Overlay of FDR Score for Target-Decoy Sets",
        xaxis_title="Correlation Score Range",
        yaxis_title="FDR (%)",
        width=900,
        height=500
    )

    # Plot Histograms
    fig_histogram = go.Figure()
    fig_histogram.add_trace(go.Histogram(
        x=target['Estimate'],
        name='Target',
        xbins=dict(start=-1, end=1, size=0.1),
        opacity=0.7,
        marker_color='blue'
    ))

    fig_histogram.add_trace(go.Histogram(
        x=decoy['Estimate'],
        name='Decoy',
        xbins=dict(start=-1, end=1, size=0.1),
        opacity=0.3,
        marker_color='red'
    ))

    fig_histogram.update_layout(
        title="Histogram of Target vs. Decoy Scores",
        xaxis_title="Correlation Score Range",
        yaxis_title="Frequency",
        xaxis_range=[-1, 1],
        barmode='overlay',
        width=900,
        height=500
    )

    return fig_fdr, fig_histogram

