from typing import Dict, List, Tuple, Literal

import numpy as np
import pandas as pd

from calculations.general import variance_per_replica_step
from funcs import FragmentationResult, assign_node_numbers, percent_step


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



def process_variance(
    data: Dict[str, FragmentationResult],
    fragmentation_types: List[str]
) -> pd.DataFrame:
    """
    Prepare concatenated per-replica variance data.

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
