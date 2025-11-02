import pickle
import matplotlib.pyplot as plt
from typing import List, Dict
import logging
import gc

from funcs import FragmentationResult, percent_step
from mean_genetics import mean_het_fst

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def load_single_file(sigma_value: str, frag_types: List[str]) -> Dict[str, FragmentationResult]:
    """Load a single sigma file for specified fragmentation types"""
    results = {}
    for ft in frag_types:
        file_path = f"RGG, {ft}_asymm_sig{sigma_value}.pickle"
        with open(file_path, "rb") as f:
            raw = pickle.load(f)
        results[ft] = FragmentationResult(
            n_steps=raw[0],
            networks=raw[1],
            het_dist=raw[2],
            het_mean=raw[3],
            fst_dist=raw[4],
            fst_mean=raw[5],
            coalescence_list=raw[6],
            fst_matrices=raw[7],
        )
    logger.info(f"Loaded data for sigma {sigma_value}")
    return results


def plot_combined_genetics():
    """Plot het and fst for multiple sigma values in a 4×2 grid"""
    fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
    fragmentation_types = ['rand', 'cor']
    sigma_values = ["00", "01", "03", "05"]
    measures = ['het', 'fst']

    fig, axes = plt.subplots(4, 2, figsize=(16, 20))

    plt.subplots_adjust(
        hspace=0.4,  # Vertical space between subplots
        wspace=0.3   # Horizontal space between subplots
    )
    
    for i, sigma in enumerate(sigma_values):
        # Load single file
        data = load_single_file(sigma, fragmentation_types)

        for j, measure in enumerate(measures):
            # Set current subplot
            plt.sca(axes[i, j])

            # Process data for plotting
            df = mean_het_fst(data, measure)

            # Plot on the current subplot
            plot_genetics_on_axis(df, measure, f"σ = 0.{sigma}")

        # Clear memory after plotting both measures
        del data
        gc.collect()

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('figs/combined_genetics_by_sigma.svg', dpi=300)
    # plt.show()


def plot_genetics_on_axis(df, measure, title):
    """Plot genetics data on the current axis"""
    import seaborn as sns

    sns.lineplot(
        data=df,
        x='step_pct',
        y='avg',
        hue='frag_type',
        estimator='mean',
        errorbar="sd"
    )
    plt.xlabel('% fragmentation', fontsize=30)
    if measure == 'het':
        plt.ylabel('Heterozygosity', fontsize=30)
    else:  # measure == 'fst'
        plt.ylabel('F_ST', fontsize=30)
    plt.title(title, fontsize=20)
    plt.tick_params(axis='both', labelsize=25)

    plt.savefig('SUP_genetics_sigma.svg', dpi=300)
   # plt.show()
    
if __name__ == "__main__":
    plot_combined_genetics()