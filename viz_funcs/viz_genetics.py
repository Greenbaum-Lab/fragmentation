import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict

from data_manipulation.genetics import process_frag_types, compute_histogram
from funcs import FragmentationResult, percent_step
import numpy as np




def plot_genetics(
    data: Dict[str, FragmentationResult],
    measure: str
):
    """
    Plot mean ± SD of the specified measure across all fragmentation types.

    :param data: Mapping from frag_type to FragmentationResult.
    :param measure: 'het' or 'fst'.
    """
    # Process all frag types to get a unified DataFrame
    df = process_frag_types(data, measure)

    # Plot using seaborn's built-in estimator for mean ± SD
    plt.figure(figsize=(10, 6))
    sns.lineplot(
        data=df,
        x='step_pct',
        y='avg',
        hue='frag_type',
        estimator='mean',
        errorbar='sd'
    )
    plt.xlabel('% fragmentation', fontsize=30)
    plt.ylabel(measure.capitalize(), fontsize=30)
    plt.tick_params(axis='both', labelsize=25)
    plt.legend(title='Type')
    plt.tight_layout()
    plt.savefig(f'./figs/genetics_{measure}.svg', format="svg")
    plt.show()


########### distributions ###########



def plot_distribution(
    df: pd.DataFrame,
    measure: str,
    frag_type: str,
) -> None:
    """
    Plot a ridgeline histogram of heterozygosity for one fragmentation type.

    :param df: DataFrame with columns ['step_pct', 'het'].
    :param frag_type: Identifier for the fragmentation type.
    """
    # 1. Compute histogram layers (reversed so lowest step at top)
    steps, bin_edges, hist_counts = compute_histogram(df, measure=measure)
    # 2. Colors reversed for top-down
    n = len(steps)
    if measure == 'het':
        cmap = plt.get_cmap('YlGnBu')(np.linspace(0, 1, n))
    else:
        cmap = plt.get_cmap('YlOrRd')(np.linspace(0, 1, n))
    # 3. Plot bars with offsets
    fig, ax = plt.subplots(figsize=(4, 2 + 0.5 * n))
    bin_width = bin_edges[1] - bin_edges[0]
    for i, (step, counts) in enumerate(zip(steps, hist_counts)):
        base = i * 6
        ax.bar(
            bin_edges[:-1],
            counts,
            width=bin_width,
            bottom=base,
            color=cmap[i],
            edgecolor='black',
            alpha=0.6,
            align='edge'
        )
        ax.hlines(base, bin_edges[0], bin_edges[-1], color='black', linewidth=0.5)

    ax.set_yticks([])
    ax.set_xlabel('Heterozygosity', fontsize=14)
    ax.set_xlim(bin_edges[0], bin_edges[-1])
    ax.set_ylim(0, 6 * n + max(cnt.max() for cnt in hist_counts))
    ax.tick_params(axis='both', labelsize=12)
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)

    plt.title(f"{frag_type}")
    plt.show()


############### individual nodes ########################


def plot_het_nodes(
    df: pd.DataFrame,
    n_nodes: int = 10,
) -> None:
    """
    Plot the heterozygosity for selected nodes across steps.

    :param df: DataFrame with 'step', 'node_replica_id', and 'het' values.
    :param n_nodes: Number of nodes to plot (choose top `n_nodes` nodes based on their node_replica_id).
    :param measure: The column to plot ('het' or 'fst').
    :param title: The plot's title.
    """
    node_ids = df['id'].unique()

    # Sample n_nodes
    selected_nodes = np.random.choice(node_ids, n_nodes, replace=False)

    df = percent_step(df, step_col='step', pct_col='step_pct')
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(6,4))

    for node_id in selected_nodes:
        # Filter data for each node_replica_id
        node_data = df[df['id'] == node_id]

        # Plot the line for the node's data
        ax.plot(node_data['step_pct'], node_data['het'],color='grey', alpha=0.5)

    # Customize plot
    ax.set_xlabel('Time', fontsize=16)
    ax.set_ylabel("Heterozygosity", fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.show()




def plot_variance(df: pd.DataFrame) -> None:
    """
    Plot variance.

    :param df: DataFrame with columns ['fragmentation_type', 'replica', 'step', 'variance'].
    """
    plt.figure(figsize=(10, 6))
    # Normalize step to percentage
    df = percent_step(df, step_col='step', pct_col='step_pct')
    sns.lineplot(
        data=df,
        x='step_pct',
        y='variance',
        hue='frag_type',
        estimator='mean',
        errorbar='sd'
    )
    plt.xlabel('% fragmentation', fontsize=25)
    plt.ylabel('Variance', fontsize=25)
    plt.tick_params(axis='both', labelsize=20)
    plt.savefig('./figs/variance.svg', format='svg', dpi=300)
    plt.show()



