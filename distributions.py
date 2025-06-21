
####### distributions of heterozygosity and fst for single fragmentation type #######

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

