
from typing import Dict, List, Tuple, Literal

import numpy as np
import pandas as pd

from funcs import FragmentationResult, assign_node_numbers, percent_step


####### main plotting of heterozygosity and fst across fragmentation types #######
def process_het_fst(
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



import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict

from data_manipulation.manp_genetics import process_frag_types, compute_histogram
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
    # plt.savefig(f'./figs/genetics_{measure}.svg', format="svg")
    plt.show()

