from typing import List, Dict, Literal

import networkx as nx
import pandas as pd
from scipy.stats import pearsonr

from funcs import FragmentationResult


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
