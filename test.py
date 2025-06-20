
from typing import Dict, List, Tuple, Literal

import numpy as np
import pandas as pd

from funcs import FragmentationResult, assign_node_numbers, percent_step


####### main plotting of heterozygosity and fst across fragmentation types #######
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
