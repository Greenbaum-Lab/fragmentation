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
from funcs import load_data, FragmentationResult, assign_node_numbers, percent_step




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

    # Ensure frag_type maintains its order
    frag_type_order = df['frag_type'].unique()
    df['frag_type'] = pd.Categorical(df['frag_type'], categories=frag_type_order, ordered=True)

    grouped = df.groupby(['frag_type', 'replica', 'step'])

    for (frag_type, replica, step), group in grouped:
        group = group[group[centrality] != 0]  # Exclude rows where centrality is 0
        if group[centrality].nunique() < 2:
            continue
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


def plot_correlation(
    corr_df,
    output_path
):
    """
    Plot correlation coefficient r over steps using Seaborn to compute mean ± SD.

    :param corr_df: DataFrame with columns ['frag_type', 'replica', 'step', 'r', 'p'].
    :param frag_type_col: Column name for fragmentation type.
    :param step_col: Column name for step.
    :param r_col: Column name for correlation coefficient.
    :param output_path: Path to save plot.
    """
    # Convert step to percentage using func percent_step
    corr_df = percent_step(corr_df, step_col='step', pct_col='step_pct')

    plt.figure(figsize=(6, 4))
    sns.lineplot(
        data=corr_df,
        x='step_pct',
        y='r',
        hue='frag_type',
        estimator='mean',
        errorbar='sd',
    )
    plt.xlabel('% fragmentation', fontsize=16)
    plt.ylabel('Correlation (r)', fontsize=16)
    plt.tick_params(axis='both', labelsize=14)
    plt.ylim(-1, 1.1)
    plt.legend().set_visible(False)
    plt.savefig(output_path, format='svg')
    plt.show()



def filter_correlations(
    corr_df: pd.DataFrame,
    min_replicates: int
) -> pd.DataFrame:
    """
    Filter correlation DataFrame to include only significant results (p < threshold)
    and groups with more than min_replicates.

    :param corr_df: DataFrame with correlation results, including p-values.
    :param min_replicates: Minimum number of replicates required per (frag_type, step).
    :return: Filtered DataFrame.
    """
    df_filtered = corr_df[(corr_df['p'] < 0.05) & (corr_df['p'] > 0)]
    # Identify valid (frag_type, step) groups with enough replicates
    valid_groups = (
        df_filtered
        .groupby(['frag_type', 'step'])['replica']
        .nunique()
        .reset_index()
        .query(f"replica >= {min_replicates}")
        [['frag_type', 'step']]
    )

    return df_filtered.merge(valid_groups, on=['frag_type', 'step'], how='inner')

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

##### compute correlation between centrality and heterozygosity
# centrality_df = pd.read_csv('./csv_new/centrality_het.csv')
# corr_df = compute_het_central_correlation(
#     df=centrality_df,
#     centrality='betweenness',
# )

##### plot correlation between centrality and heterozygosity
# corr_df = pd.read_csv('./csv_new/het_bet_correlation.csv')
# filtered_corr_df = filter_correlations(corr_df, min_replicates=5)
# print(filtered_corr_df)
# plot_correlation(
#     corr_df=filtered_corr_df,
#     output_path='./figs/corr_bet.svg'
# )



############ correlation of sigle steps
# def add_annotation(ax, r: float, p: float) -> None:
#     """
#     Add annotation with correlation coefficient and p-value to the plot.
#
#     :param ax: Matplotlib Axes object.
#     :param r: Pearson correlation coefficient.
#     :param p: P-value of the correlation.
#     """
#     annotation_text = f'r = {r:.2f}\np = {p:.3f}' if p >= 0.001 else f'r = {r:.2f}\np < 0.001'
#     ax.annotate(
#         annotation_text,
#         xy=(0.5, 0.1),  # Position (x, y) as relative plot coordinates
#         xycoords='axes fraction',  # Use axes fraction for relative positioning
#         fontsize=16,
#         style='italic',
#         fontname='serif'
#     )
# def plot_correlation(
#     df: pd.DataFrame,
#     measure: Literal['degree', 'betweenness'],
#     output_path: str
# ) -> None:
#     """
#     Plot the correlation between a centrality measure and heterozygosity for
#      fragmentation type-step-replica.
#
#     :param df: DataFrame containing the centrality and heterozygosity data.
#     :param measure: The centrality measure to correlate ('degree_centrality' or 'betweenness_centrality').
#     :param output_path: Path to save the plot.
#     """
#     # 1. Compute the correlation coefficient (Pearson)
#     r, p = pearsonr(df[measure], df['het'])
#
#     plt.figure(figsize=(6, 4))
#     sns.regplot(data=df, x=measure, y='het', fit_reg=True)
#
#     # 3. Annotate the plot with r and p-value
#     add_annotation(plt.gca(), r, p)
#
#     # 3. Add the correlation coefficient to the plot
#     plt.xlabel('Degree', fontsize=18)
#     plt.ylabel('Heterozygosity', fontsize=18)
#     plt.tick_params(axis='both', which='major', labelsize=14)
#     plt.ylim(-0.05, 1.2)
#
#     plt.savefig(output_path, format='svg')
#     plt.show()
# def preprocess_centrality_data(df: pd.DataFrame, replica: int, step: int, frag_type: str) -> pd.DataFrame:
#     """
#     Preprocess the centrality DataFrame by filtering for a specific replica, step, and fragmentation type.
#
#     :param df: DataFrame containing centrality data.
#     :param replica: Replica index to filter.
#     :param step: Step index to filter.
#     :param frag_type: Fragmentation type to filter.
#     :return: Filtered DataFrame.
#     """
#     filtered_df = df[(df['replica'] == replica) & (df['step'] == step)]
#     filtered_df = filtered_df[filtered_df['frag_type'] == frag_type]
#     return filtered_df
# def plot_correlation_steps(
#     df: pd.DataFrame,
#     frag_type: str,
#     replica: int,
#     steps: List[int],
#     measure: Literal['degree', 'betweenness'],
#     output_path: str
# ) -> None:
#     """
#     Produce a row of three scatter+regression plots of centrality vs. het,
#     for a single frag_type and replica, at the specified steps.
#
#     :param df: DataFrame with columns ['frag_type','replica','step','node_number',
#                'degree','betweenness','het'].
#     :param frag_type: Fragmentation type to filter on.
#     :param replica: Replica index to filter on.
#     :param steps: step indices to plot.
#     :param measure: Which centrality to plot ('degree' or 'betweenness').
#     :param output_path: Where to save the combined figure.
#     """
#     # set up 1×3 axes
#     fig, axes = plt.subplots(1, 3, figsize=(10, 2), sharey=True)
#     for ax, step in zip(axes, steps):
#         # filter for frag_type, replica, and step
#         sub = preprocess_centrality_data(df, replica, step, frag_type)
#
#         sns.regplot(
#             data=sub,
#             x=measure,
#             y='het',
#             ax=ax,
#             scatter_kws={'alpha':0.7},
#         )
#
#         # compute and annotate r & p
#         r, p = pearsonr(sub[measure], sub['het'])
#         add_annotation(ax, r, p)
#
#         # styling
#         ax.set_xlabel(measure.capitalize(), fontsize=14)
#         if ax is axes[0]:
#             ax.set_ylabel("Heterozygosity", fontsize=14)
#         else:
#             ax.set_ylabel("")
#
#         ax.tick_params(labelsize=12)
#         ax.set_ylim(-0.1, 1.4)
#
#     plt.savefig(output_path, format='svg')
#     plt.show()
#
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
