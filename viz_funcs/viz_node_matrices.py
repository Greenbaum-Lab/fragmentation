from matplotlib import pyplot as plt

from funcs import percent_step


def plot_correlation(
    corr_df,
    output_path
, sns=None):
    """
    Plot correlation coefficient r over steps using Seaborn to compute mean ± SD.

    :param corr_df: DataFrame with columns ['frag_type', 'replica', 'step', 'r', 'p'].
    :param frag_type_col: Column name for fragmentation type.
    :param step_col: Column name for step.
    :param r_col: Column name for correlation coefficient.
    :param output_path: Path to save plot.
    """
    # Convert step to percentage using func percent_step
    corr_df = percent_step(corr_df, step_col='step', pct_col='step_pct')

    plt.figure(figsize=(6, 4))
    sns.lineplot(
        data=corr_df,
        x='step_pct',
        y='r',
        hue='frag_type',
        estimator='mean',
        errorbar='sd',
    )
    plt.xlabel('% fragmentation', fontsize=16)
    plt.ylabel('Correlation (r)', fontsize=16)
    plt.tick_params(axis='both', labelsize=14)
    plt.ylim(-1, 1.1)
    plt.legend().set_visible(False)
    plt.savefig(output_path, format='svg')
    plt.show()