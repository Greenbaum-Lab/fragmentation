from typing import Dict

import pandas as pd

from funcs import FragmentationResult, assign_node_numbers


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
