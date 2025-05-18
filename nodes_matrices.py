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


#
#
#
# def plot_node_centrality(df, step: int, centrality: str, frag: str, log=bool):
#
#     sns.regplot(x='central', y='het', data=df, fit_reg=True, order=1,logistic=True,
#                 scatter_kws={'s': 50, 'alpha': 0.7, 'color': 'blue'})
#     plt.ylabel("Heterozygosity", fontsize=18)
#
#     plt.tick_params(axis='both', labelsize=16)
#     plt.savefig(f'./figs/node_{centrality}_{step}_{frag}.svg', format="svg")
#     plt.show()
#
#
# def compute_correlation(df):
#     """
#     Calculate the Pearson correlation coefficient between 'het' and 'central'
#     for each combination of 'step' and 'replica' in the DataFrame.
#
#     Parameters:
#         df (pd.DataFrame): Input DataFrame containing 'step', 'het', 'replica', and 'central' columns.
#
#     Returns:
#         pd.DataFrame: DataFrame with columns 'step', 'replica', and 'cor', containing the correlation values.
#     """
#     results = []
#     grouped = df.groupby(['step', 'replica'])
#
#     for (step, replica), group in grouped:
#         if len(group) < 2: # need at least 2 data points to calculate correlation
#             continue
#         cor, pval = pearsonr(group['het'], group['central'])
#         results.append({'step': step, 'replica': replica, 'cor': cor, 'pval': pval})
#
#     results_df = pd.DataFrame(results)
#     return results_df


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


##### plot single correlation, degree-heterozygosity or betweenness-heterozygosity
# df = pd.read_csv(f'./csv/bet_het_dist.csv')
#
# steps = [0, 75, 150]
# fig, axes = plt.subplots(1, 3, figsize=(18, 7), sharey=True)
#
# for i, step in enumerate(steps):
#     df_step = df[(df['step'] == step) & (df['replica'] == 11)]
#     sns.regplot(x='central', y='het', data=df_step, fit_reg=True, order=1, ax=axes[i])
#     cor, pval = pearsonr(df_step['central'], df_step['het'])
#     axes[i].set_xlabel('Population betweenness', fontsize=30)
#     axes[i].set_ylabel('Heterozygosity' if i == 0 else '', fontsize=30)
#     axes[i].tick_params(axis='both', labelsize=25)
#     axes[i].set_ylim(0, 1.4)
#     axes[i].text(0.05, 1.2, f'r={cor:.2f}\np={pval:.2e}', fontsize=20, transform=axes[i].transAxes)
# plt.tight_layout()
# plt.savefig(f'./figs/degree_bet_steps.svg', format="svg")
# plt.show()



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



def plot_correlation(
    df: pd.DataFrame,
    measure: Literal['degree', 'betweenness'],
    frag_type_col: str = 'fragmentation_type',
    output_path: str = './figs/correlation_plot.svg'
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
    plt.annotate(
        f'r = {r:.2f}\np = {p:.2e}',
        xy=(0.05, 0.95),  # Position (x, y) as relative plot coordinates
        xycoords='axes fraction',  # Use axes fraction for relative positioning
        fontsize=14,
        color='black',
        ha='left',
        va='top',  # Align to the top-left corner
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.5')
    )

    # 3. Add the correlation coefficient to the plot
     plt.xlabel('Degree', fontsize=18)
    plt.ylabel('Heterozygosity', fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=14)

    # plt.savefig(output_path, format='svg')
    plt.show()




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
centrality_df = pd.read_csv('./csv_new/centrality_het.csv')
#filter data for random replica and step and frag type
centrality_df = centrality_df[(centrality_df['replica'] == 0) & (centrality_df['step'] == 0)]
centrality_df = centrality_df[centrality_df['frag_type'] == 'rand']

plot_correlation(centrality_df, measure='degree', frag_type_col='frag_type')
