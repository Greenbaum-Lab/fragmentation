# -----------------------------------------------------------------------------
# This script provides utility functions for analyzing the distributions of
# genetic diversity measures such as heterozygosity (het) and fixation index (fst)
# across different levels of network fragmentation.

# -----------------------------------------------------------------------------


import pandas as pd
import numpy as np
from typing import Literal, Tuple, List, Optional, TYPE_CHECKING

from matplotlib import pyplot as plt

from funcs import FragmentationResult, percent_step


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


