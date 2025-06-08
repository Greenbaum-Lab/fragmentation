from typing import Dict, List
import pandas as pd
from funcs import FragmentationResult, giant_component_over_steps,load_data




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


def het_component(
    data: Dict[str, FragmentationResult],
    frag_type: str
) -> pd.DataFrame:
    """
    For one fragmentation type, grab its precomputed mean het per replica-step
    and its giant-component fraction, merge them, and drop any zero-component rows.

    :param data: Mapping frag_type → FragmentationResult
    :param frag_type: Key of the fragmentation type to process
    :return: DataFrame with columns
             ['replica','step','avg_het','component'], filtered to component > 0.
    """
    frag_res = data[frag_type]
    # Mean heterozygosity per replica-step
    het_rep = (
        frag_res.het_mean
    )

    # 2. Giant-component fraction per replica-step
    comp_df = giant_component_over_steps(frag_res.networks)

    # 3. Merge and drop zero-size
    merged = pd.merge(het_rep, comp_df, on=['replica', 'step'], how='inner')
    return merged[merged['component'] > 0].reset_index(drop=True)

def het_component_types(
    data: Dict[str, FragmentationResult],
    frag_types: List[str] = None
) -> Dict[str, pd.DataFrame]:
    """
    For each fragmentation type, compute the merged heterozygosity vs. giant-component
    """
    all_types = {
        ft: het_component(data, ft)
        for ft in frag_types
    }
    return all_types


def bin_het_component(
    df: pd.DataFrame,
    n_bins: int = 20
) -> pd.DataFrame:
    """
    Bin component fractions into n_bins and compute mean±sd of avg_het in each bin.
    """
    binned = pd.cut(df['component'], bins=n_bins)
    stats = (
        df
        .groupby(binned)['avg']
        .agg(mean_het='mean', sd_het='std')
        .reset_index()
    )
    stats['component_mid'] = stats['component'].apply(lambda interval: interval.mid)
    return stats[['component_mid','mean_het','sd_het']]
