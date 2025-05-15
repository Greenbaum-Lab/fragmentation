import pickle
from statistics import mean
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict
from funcs import FragmentationResult, percent_step, load_data, assign_node_numbers
from typing import Literal, List, Tuple
import numpy as np




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