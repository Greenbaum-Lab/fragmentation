import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict

from funcs import FragmentationResult, percent_step


####### main plotting of heterozygosity and fst across fragmentation types #######
def mean_het_fst(
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
        if measure == 'het':
            df = frag_res.het_mean.copy()
        elif measure == 'fst':
            df = frag_res.fst_mean.copy()
        else:
            raise ValueError(f"Unknown measure {measure!r}, expected 'het' or 'fst'.")

        df = percent_step(df, step_col='step', pct_col='step_pct')
        df['frag_type'] = frag_type
        all_types.append(df[['step_pct', 'avg', 'replica', 'frag_type']])

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
    df = mean_het_fst(data, measure)

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
    # plt.savefig(f'./figs/genetics_{measure}.svg', format="svg")
    plt.show()

