import pickle
from statistics import mean
##############


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict
from funcs import FragmentationResult, percent_step, load_data, assign_node_numbers
from typing import Literal, List, Tuple
import numpy as np

def process_frag_types(
    data: Dict[str, FragmentationResult],
    measure: str
) -> pd.DataFrame:
    """
    Combine and normalize replicate-level mean data for all fragmentation types.

    :param data: Mapping from frag_type to FragmentationResult.
    :param measure: 'het' or 'fst'.
    :return: DataFrame with columns ['step_pct', 'avg', 'replica', 'frag_type'].
    """
    all_types = []
    for frag_type, frag_res in data.items():
        # Select the appropriate summary stats for each frag_type
        if measure == 'het':
            df = frag_res.het_mean.copy()
        elif measure == 'fst':
            df = frag_res.fst_mean.copy()
        else:
            raise ValueError(f"Unknown measure {measure!r}, expected 'het' or 'fst'.")

        # Compute fragmentation percentage (0–100)
        df = percent_step(df, step_col='step', pct_col='step_pct')

        # Tag the fragmentation type
        df['frag_type'] = frag_type
        # Keep only relevant columns
        all_types.append(df[['step_pct', 'avg', 'replica', 'frag_type']])

    # Concatenate all types into one DataFrame
    return pd.concat(all_types, ignore_index=True)


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
def filter_intervals(
    frag_res: FragmentationResult,
    measure: Literal['het', 'fst'],
    interval_pct: int = 25
) -> pd.DataFrame:
    """
    Select node-level measure data at fixed fragmentation-percent intervals
    (e.g. interval_pct=25 → steps at exactly 0, 25, 50, 75, 100).

    :param frag_res: One fragmentation result.
    :param measure: Which column to filter ('het' or 'fst').
    :param interval_pct: Percentage spacing of intervals (must divide 100 evenly).
    :return: DataFrame with columns ['step_pct','replica', measure].
    """
    # 1. Pick the genetic data distribution
    df = frag_res.het_dist if measure == 'het' else frag_res.fst_dist

    # 2. Compute continuous 0–100 step_pct
    df = percent_step(df, step_col='step', pct_col='step_pct')

    # 3. Snap to nearest interval_pct multiple
    df['step_pct'] = (
        (df['step_pct'] / interval_pct)
        .round()              # round to nearest integer multiple
        .astype(int)          # cast to int
        * interval_pct
    )

    # 4. Define the exact allowed intervals
    allowed = set(range(0, 100, interval_pct))

    # 5. Filter to only those snapped intervals
    sel = df[df['step_pct'].isin(allowed)].copy()

    # 6. Return only the clean columns
    return sel[['step_pct', 'replica', measure]]


def compute_histogram(
    df: pd.DataFrame,
    measure: str,
) -> Tuple[List[int], np.ndarray, List[np.ndarray]]:
    """
    Prepare histogram data for each step_pct layer.

    :param df: DataFrame with columns ['step_pct', measure].
    :param measure: Column to histogram ('het' or 'fst').
    :return:
      - steps: sorted unique step_pct values
      - bin_edges: array of length bins+1
      - hist_counts: list of count arrays for each step
    """
    steps = sorted(df['step_pct'].unique(), reverse=True)
    hist_counts = []
    bin_edges = None

    for step in steps:
        values = df.loc[df['step_pct'] == step, measure].values
        counts, edges = np.histogram(values, bins=40, density=True)
        hist_counts.append(counts)
        bin_edges = edges

    return steps, bin_edges, hist_counts


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

def select_random_nodes(
    df: pd.DataFrame,
    per_replica: int = 1,
    nodes_per_step: int = 50
) -> Dict[int, np.ndarray]:
    """
    For each replica, choose `per_replica` random node indices.

    :param df: DataFrame containing the heterozygosity data.
    :param per_replica: Number of random nodes to select per replica.
    :param nodes_per_step: Number of nodes in each step.
    :return: A dictionary with replica ids as keys and arrays of selected node indices as values.
    """
    # Ensure `node_number` is assigned
    df = assign_node_numbers(df, nodes_per_step)

    selections = {}
    for rep, sub in df.groupby("replica"):
        n_nodes = sub["node_number"].nunique()
        picks = np.random.choice(n_nodes, min(per_replica, n_nodes), replace=False)
        selections[int(rep)] = picks
    return selections


def extract_nodes(
    df: pd.DataFrame,
    selections: Dict[int, np.ndarray],
    nodes_per_step: int = 50
) -> pd.DataFrame:
    """
    Extract the heterozygosity data of selected nodes for each replica and step.

    :param df: DataFrame with 'node_number', 'step', 'replica', and 'het' values.
    :param selections: A dictionary with replicas as keys and lists of selected node indices as values.
    :param nodes_per_step: Number of nodes per step in the data (should match original assignment).
    :return: DataFrame containing only selected nodes, including a 'node_replica_id'.
    """
    out = []
    for rep, nodes in selections.items():
        sub = df[df["replica"] == rep]
        for node in nodes:
            node_df = sub[sub["node_number"] == node].copy()
            node_df["id"] = f"n{node}_r{rep}"
            out.append(node_df)

    return pd.concat(out, ignore_index=True).drop(columns=['replica', 'node_number'])


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


############### variance ####################

def variance_per_replica_step(
    data: Dict[str, FragmentationResult],
    frag_type: str
) -> pd.DataFrame:
    """
    Calculate heterozygosity variance per (replica, step) for a given frag_type.

    :param data: Dict of fragmentation results.
    :param frag_type: Fragmentation type key.
    :return: DataFrame with columns ['replica', 'step', 'variance'].
    """
    frag_res = data[frag_type]
    df = frag_res.het_dist
    return (
        df.groupby(['replica', 'step'])['het']
          .var(ddof=1)
          .reset_index(name='variance')
    )


def process_variance(
    data: Dict[str, FragmentationResult],
    fragmentation_types: List[str]
) -> pd.DataFrame:
    """
    Prepare concatenated per-replica variance data for Seaborn.

    :param data: Dict of frag_type → FragmentationResult.
    :param fragmentation_types: List of frag_types to process.
    :return: DataFrame with columns ['fragmentation_type', 'replica', 'step', 'variance'].
    """
    dfs = []
    for frag_type in fragmentation_types:
        var_df = variance_per_replica_step(data, frag_type)
        var_df['frag_type'] = frag_type
        dfs.append(var_df)
    return pd.concat(dfs, ignore_index=True)


def plot_variance(df: pd.DataFrame) -> None:
    """
    Plot variance using Seaborn’s estimator and errorbar.

    :param df: DataFrame with columns ['fragmentation_type', 'replica', 'step', 'variance'].
    :param step_col: Column name for x axis.
    :param var_col: Column name for variance values.
    :param frag_type_col: Column name for hue.
    :param output_path: Path to save the figure.
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



#######################
####################### plot data
# fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt','wrst']
# data = load_data(fragmentation_types)
# plot_genetics(data, measure='het')

##############plot distributions
##### one frag type each time
# fragmentation_types = ['rand']
# data = load_data(fragmentation_types)
# df = filter_intervals(data['rand'], measure='fst', interval_pct=25)
# plot_distribution(df,measure='fst', frag_type='rand')

####################### plot individual nodes
################
# fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt', 'wrst']
# fragmentation_types = ['rand']
# data = load_data(fragmentation_types)
# df = data.get('rand').het_dist
# x=select_random_nodes(df, per_replica=1)
# df_selected = extract_nodes(df, x)
# plot_het_nodes(df_selected, n_nodes=10)

############# calculate and plot variance across nodes in the network
fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt','wrst']
data = load_data(fragmentation_types)
#for single frag type use fragmentation_types[x]
var = process_variance(data, fragmentation_types)
plot_variance(var)