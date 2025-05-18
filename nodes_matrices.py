import math
import pickle
from statistics import mean
from typing import List, Dict, Literal

import networkx as nx
import numpy as np
import pandas as pd
from infomap import Infomap
from joypy import joyplot
from matplotlib import pyplot as plt
import seaborn as sns
from mantel import test
from scipy.stats import pearsonr
from scipy.stats import norm

from funcs import load_data, FragmentationResult, assign_node_numbers


def compute_centrality_network(graph: nx.Graph) -> pd.DataFrame:
    """
    Compute degree and betweenness centrality for all nodes in a network.

    :param graph: NetworkX graph instance.
    :return: DataFrame with columns ['node_number', 'degree_centrality', 'betweenness_centrality'].
    """
    # degree_centrality = nx.degree_centrality(graph)
    betweenness_centrality = nx.betweenness_centrality(graph)
    degree_centrality = dict(nx.degree(graph))

    #         'degree': lambda net: dict(nx.degree(net))

    df = pd.DataFrame({
        'node_number': list(degree_centrality.keys()),
        'degree': list(degree_centrality.values()),
        'betweenness': list(betweenness_centrality.values())
    })

    return df


def compute_centrality_replicates(
    networks: List[List[nx.Graph]]
) -> pd.DataFrame:
    """
    Compute node centralities for all replicate-step networks.

    :param networks: Nested list of graphs [replicate][step].
    :return: DataFrame with columns ['replica', 'step', 'node_number', 'degree', 'betweenness'].
    """
    records = []

    for replica_idx, replicate_graphs in enumerate(networks):
        for step_idx, graph in enumerate(replicate_graphs):
            centralities_df = compute_centrality_network(graph)
            centralities_df['replica'] = replica_idx
            centralities_df['step'] = step_idx
            records.append(centralities_df)

    return pd.concat(records, ignore_index=True)


def compute_centrality_types(
    data: Dict[str, FragmentationResult],
    frag_types: list[str]
) -> pd.DataFrame:
    """
    Compute degree and betweenness centralities for all graphs across multiple fragmentation types.

    :param data: Mapping from frag_type to FragmentationResult.
    :param frag_types: List of fragmentation type keys to process.
    :return: DataFrame with columns ['fragmentation_type', 'replica', 'step', 'node_number', 'degree_centrality', 'betweenness_centrality'].
    """
    all_dfs = []
    for frag_type in frag_types:
        frag_res = data[frag_type]
        df = compute_centrality_replicates(frag_res.networks)
        df['frag_type'] = frag_type
        all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)
    cols = ['frag_type', 'replica', 'step', 'node_number', 'degree', 'betweenness']
    combined_df.to_csv(f'./csv_new/centrality.csv', index=False)
    return combined_df[cols]


def merge_centrality_het(
        centrality_df: pd.DataFrame,
        data: Dict[str, FragmentationResult],
        frag_types: list[str]
) -> pd.DataFrame:
    """
    Preprocess and merge the centrality data with heterozygosity data for each fragmentation type.

    :param centrality_df: DataFrame containing 'frag_type', 'replica', 'step', 'node_number', 'degree', 'betweenness'
    :param data: Dictionary mapping frag_type → FragmentationResult
    :param frag_types: List of fragmentation types to process
    :return: Merged DataFrame with centrality and heterozygosity for each node.
    """
    all_data = []

    # Iterate over each fragmentation type
    for frag_type in frag_types:
        # Get the heterozygosity data from FragmentationResult
        frag_res = data[frag_type]
        assign_node_numbers(frag_res.het_dist)
        het_df = frag_res.het_dist

        # Merge the centrality and heterozygosity data on ['replica', 'step', 'node_number']
        merged_df = pd.merge(
            centrality_df[centrality_df['frag_type'] == frag_type],
            het_df[['replica', 'step', 'node_number', 'het']],
            on=['replica', 'step', 'node_number'],
            how='left'  # 'left' join keeps all centrality data and adds 'het' where possible
        )

        all_data.append(merged_df)

    # Concatenate all fragmentation types into a single DataFrame
    final_df = pd.concat(all_data, ignore_index=True)
    final_df.to_csv(f'./csv_new/centrality_het.csv', index=False)

    return final_df




##### plot heterozygisuty vs. node centrality

# fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
# data = load_data(fragmentation_types)

### compute correlation between heterozygosity and centrality for all processes
### change centrality name
# make_het_central_all(data, 'bet')
# for frag in fragmentation_types:
#     df = pd.read_csv(f'./csv/bet_het_{frag}.csv')
#     results_df = compute_correlation(df)
#     results_df.to_csv(f'./csv/bet_het_{frag}_cor.csv', index=False)



#### plot correlation het-degree for all replicates during fragmetiaon
# plt.figure(figsize=(10, 6))
# for frag in fragmentation_types:
#     df = pd.read_csv(f'./csv/bet_het_{frag}_cor.csv')
#     df['step'] = df['step'] / df['step'].max() * 100  # Normalize steps to percentage
#
#     df = df[df['pval'] < 0.05]  # Filter rows with pval < 0.05
#     # Filter steps with at least 5 unique replicas
#     valid_steps = df.groupby('step')['replica'].nunique()
#     valid_steps = valid_steps[valid_steps >= 5].index
#     df = df[df['step'].isin(valid_steps)]
#     sns.lineplot(x='step', y='cor', data=df, errorbar='sd')
#
# plt.xlabel('% fragmentation', fontsize=28)
# plt.ylabel('Correlation (r)', fontsize=28)
# plt.tick_params(axis='both', labelsize=25)
# plt.ylim(-1.05, 1.1)
# # plt.savefig(f'./figs/bet_het_cor_pval.svg', format="svg")
# plt.show()
#





def add_annotation(ax, r: float, p: float) -> None:
    """
    Add annotation with correlation coefficient and p-value to the plot.

    :param ax: Matplotlib Axes object.
    :param r: Pearson correlation coefficient.
    :param p: P-value of the correlation.
    """
    annotation_text = f'r = {r:.2f}\np = {p:.3f}' if p >= 0.001 else f'r = {r:.2f}\np < 0.001'
    ax.annotate(
        annotation_text,
        xy=(0.5, 0.1),  # Position (x, y) as relative plot coordinates
        xycoords='axes fraction',  # Use axes fraction for relative positioning
        fontsize=16,
        style='italic',
        fontname='serif'
    )


def plot_correlation(
    df: pd.DataFrame,
    measure: Literal['degree', 'betweenness'],
    output_path: str
) -> None:
    """
    Plot the correlation between a centrality measure and heterozygosity for
     fragmentation type-step-replica.

    :param df: DataFrame containing the centrality and heterozygosity data.
    :param measure: The centrality measure to correlate ('degree_centrality' or 'betweenness_centrality').
    :param output_path: Path to save the plot.
    """
    # 1. Compute the correlation coefficient (Pearson)
    r, p = pearsonr(df[measure], df['het'])

    plt.figure(figsize=(6, 4))
    sns.regplot(data=df, x=measure, y='het', fit_reg=True)

    # 3. Annotate the plot with r and p-value
    add_annotation(plt.gca(), r, p)

    # 3. Add the correlation coefficient to the plot
    plt.xlabel('Degree', fontsize=18)
    plt.ylabel('Heterozygosity', fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.ylim(-0.05, 1.2)

    plt.savefig(output_path, format='svg')
    plt.show()

def preprocess_centrality_data(df: pd.DataFrame, replica: int, step: int, frag_type: str) -> pd.DataFrame:
    """
    Preprocess the centrality DataFrame by filtering for a specific replica, step, and fragmentation type.

    :param df: DataFrame containing centrality data.
    :param replica: Replica index to filter.
    :param step: Step index to filter.
    :param frag_type: Fragmentation type to filter.
    :return: Filtered DataFrame.
    """
    filtered_df = df[(df['replica'] == replica) & (df['step'] == step)]
    filtered_df = filtered_df[filtered_df['frag_type'] == frag_type]
    return filtered_df


def plot_correlation_steps(
    df: pd.DataFrame,
    frag_type: str,
    replica: int,
    steps: List[int],
    measure: Literal['degree', 'betweenness'],
    output_path: str
) -> None:
    """
    Produce a row of three scatter+regression plots of centrality vs. het,
    for a single frag_type and replica, at the specified steps.

    :param df: DataFrame with columns ['frag_type','replica','step','node_number',
               'degree','betweenness','het'].
    :param frag_type: Fragmentation type to filter on.
    :param replica: Replica index to filter on.
    :param steps: step indices to plot.
    :param measure: Which centrality to plot ('degree' or 'betweenness').
    :param output_path: Where to save the combined figure.
    """
    # set up 1×3 axes
    fig, axes = plt.subplots(1, 3, figsize=(10, 2), sharey=True)
    for ax, step in zip(axes, steps):
        # filter for frag_type, replica, and step
        sub = preprocess_centrality_data(df, replica, step, frag_type)

        sns.regplot(
            data=sub,
            x=measure,
            y='het',
            ax=ax,
            scatter_kws={'alpha':0.7},
        )

        # compute and annotate r & p
        r, p = pearsonr(sub[measure], sub['het'])
        add_annotation(ax, r, p)

        # styling
        ax.set_xlabel(measure.capitalize(), fontsize=14)
        if ax is axes[0]:
            ax.set_ylabel("Heterozygosity", fontsize=14)
        else:
            ax.set_ylabel("")

        ax.tick_params(labelsize=12)
        ax.set_ylim(-0.1, 1.4)

    plt.savefig(output_path, format='svg')
    plt.show()


def compute_het_central_correlation(
    df: pd.DataFrame,
    centrality: Literal['degree', 'betweenness'],
) -> pd.DataFrame:
    """
    Compute Pearson correlation (r) and p-value between centrality and heterozygosity
    for each (fragmentation_type, replica, step) group.

    :param df: DataFrame containing columns:
               ['frag_type', 'replica', 'step', centrality_col, heterozygosity_col]
    :param centrality_col: Name of centrality measure column ('degree' or 'betweenness')
    :param heterozygosity_col: Name of heterozygosity column (default 'het')
    :return: DataFrame with columns:
             ['frag_type', 'replica', 'step', 'r', 'p']
    """
    results = []

    grouped = df.groupby(['frag_type', 'replica', 'step'])

    for (frag_type, replica, step), group in grouped:
        r, p = pearsonr(group[centrality], group['het'])
        results.append({
            'frag_type': frag_type,
            'replica': replica,
            'step': step,
            'r': r,
            'p': p
        })

    corr_df = pd.DataFrame(results)
    corr_df.to_csv(f'./csv_new/het_bet_correlation.csv', index=False)
    return pd.DataFrame(results)


###scripts
###### compute centrality for all fragmentation types
# fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
# data = load_data(fragmentation_types)
# centrality_df = compute_centralities_types(data, fragmentation_types)


##### merge centrality with heterozygosity data
# frag_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
# data = load_data(frag_types)
# centrality_df = pd.read_csv('./csv_new/centrality.csv')
# merged_df = merge_centrality_het(centrality_df, data, frag_types)


##### plot centrality vs. heterozygosity
# centrality_df = pd.read_csv('./csv_new/centrality_het.csv')
# steps = [0, 75, 150]
# plot_correlation_steps(
#     df=centrality_df,
#     frag_type='dist',
#     replica=10,
#     steps=steps,
#     measure='betweenness',
#     output_path='./figs/het_bet_steps.svg'
# )

centrality_df = pd.read_csv('./csv_new/centrality_het.csv')

corr_df = compute_het_central_correlation(
    df=centrality_df,
    centrality='betweenness',
)

print(corr_df)
