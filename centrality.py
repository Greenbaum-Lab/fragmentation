"""
centrality.py

This module provides functions to compute node centralities (degree and betweenness)
for networks and collections of networks, especially in the context of network fragmentation
analyses. It uses NetworkX for graph processing and pandas for data handling. The main
functions allow centrality calculations for a single network, across replicates (e.g., 
simulation steps), and across multiple fragmentation types.
"""

from typing import List, Dict, Literal

import networkx as nx
import pandas as pd
from scipy.stats import pearsonr

from funcs import FragmentationResult


def compute_centrality_network(graph: nx.Graph) -> pd.DataFrame:
    if graph.is_directed():
        degree = dict(nx.degree(graph, weight='weight'))
        H = graph.copy()
        for u, v, d in H.edges(data=True):
            w = float(d.get('weight', 1.0))
            d['length'] = 1.0 / w
        betweenness = nx.betweenness_centrality(H, weight='length', normalized=True)
    else:
        degree = dict(nx.degree(graph, weight='weight'))
        betweenness = nx.betweenness_centrality(graph, weight='weight', normalized=True)
    return pd.DataFrame({
        'node_number': list(graph.nodes()),
        'degree': [degree[n] for n in graph.nodes()],
        'betweenness': [betweenness[n] for n in graph.nodes()],
    })

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
    frag_types: List[str]
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
    print(combined_df)
    cols = ['frag_type', 'replica', 'step', 'node_number', 'degree', 'betweenness']
    combined_df.to_csv(f'centrality_sig03.csv', index=False)
    return combined_df[cols]

