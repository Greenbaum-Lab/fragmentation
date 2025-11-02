# import pickle
# import matplotlib.pyplot as plt
# from typing import List, Dict
# import logging
# import gc
#
# from funcs import FragmentationResult, percent_step
# from mean_genetics import mean_het_fst
#
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)
#
#
# def load_single_file(sigma_value: str, frag_types: List[str]) -> Dict[str, FragmentationResult]:
#     """Load a single sigma file for specified fragmentation types"""
#     results = {}
#     for ft in frag_types:
#         file_path = f"RGG, {ft}_asymm_sig{sigma_value}.pickle"
#         with open(file_path, "rb") as f:
#             raw = pickle.load(f)
#         results[ft] = FragmentationResult(
#             n_steps=raw[0],
#             networks=raw[1],
#             het_dist=raw[2],
#             het_mean=raw[3],
#             fst_dist=raw[4],
#             fst_mean=raw[5],
#             coalescence_list=raw[6],
#             fst_matrices=raw[7],
#         )
#     logger.info(f"Loaded data for sigma {sigma_value}")
#     return results
#
#
# def plot_combined_genetics():
#     """Plot het and fst for multiple sigma values in a 4×2 grid"""
#     fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
#     sigma_values = ["00", "01", "03", "05"]
#     measures = ['het', 'fst']
#
#     fig, axes = plt.subplots(4, 2, figsize=(12, 18))
#
#     plt.subplots_adjust(
#         hspace=0.4,  # Vertical space between subplots
#         wspace=0.3  # Horizontal space between subplots
#     )
#
#     for i, sigma in enumerate(sigma_values):
#         # Load single file
#         data = load_single_file(sigma, fragmentation_types)
#
#         for j, measure in enumerate(measures):
#             # Set current subplot
#             plt.sca(axes[i, j])
#
#             # Process data for plotting
#             df = mean_het_fst(data, measure)
#
#             # Plot on the current subplot
#             plot_genetics_on_axis(df, measure, f"σ = 0.{sigma}")
#
#         # Clear memory after plotting both measures
#         del data
#         gc.collect()
#
#     # Adjust layout and save
#     plt.tight_layout()
#     plt.savefig('figs/combined_genetics_by_sigma.svg', dpi=300)
#     # plt.show()
#
#
# def plot_genetics_on_axis(df, measure, title):
#     """Plot genetics data on the current axis"""
#     import seaborn as sns
#
#     sns.lineplot(
#         data=df,
#         x='step_pct',
#         y='avg',
#         hue='frag_type',
#         estimator='mean',
#         errorbar="sd",
#         legend=False
#     )
#     plt.xlabel('% fragmentation', fontsize=18)
#     if measure == 'het':
#         plt.ylabel('Heterozygosity', fontsize=18)
#     else:  # measure == 'fst'
#         plt.ylabel('F$_{ST}$', fontsize=18)
#         plt.ylim(-0.05, 1.05)
#     plt.title(title, fontsize=20)
#     plt.tick_params(axis='both', labelsize=15)
#
#     plt.savefig('SUP_genetics_sigma.svg', dpi=300)
#
#
# # plt.show()
#
# if __name__ == "__main__":
#     plot_combined_genetics()


import pickle
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from typing import Dict, List

from funcs import FragmentationResult
from centrality_corr import compute_het_central_correlation, filter_correlations, merge_centrality_het
from centrality import compute_centrality_types


def load_sigma_data(sigma_value: str, frag_types: List[str]) -> Dict[str, FragmentationResult]:
    """Load data for a specific sigma value and fragmentation types"""
    results = {}
    for ft in frag_types:
        file_path = f"../RGG, {ft}_asymm_sig{sigma_value}.pickle"
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
    print(f"Loaded data for sigma {sigma_value}")
    return results


def run_correlation_pipeline(data, frag_types, centrality_type, sigma_value, ax):
    """Run the centrality correlation pipeline for specific data and plot on given axis"""
    # Compute centrality
    centrality_df = compute_centrality_types(data, frag_types, centrality_types=[centrality_type])

    # Merge centrality with heterozygosity data
    merged_df = merge_centrality_het(centrality_df, data, frag_types)

    # Compute correlation
    corr_df = compute_het_central_correlation(df=merged_df, centrality=centrality_type)

    # Filter correlations
    filtered_corr_df = filter_correlations(corr_df, min_replicates=5)

    # Plot on the given axis
    plt.sca(ax)
    # Modify plot_correlation to work with an axis rather than creating a new figure
    import seaborn as sns
    sns.lineplot(
        data=filtered_corr_df,
        x='step_pct',
        y='correlation',
        hue='frag_type',
        errorbar=('ci', 95),
        ax=ax
    )
    ax.set_xlabel('% fragmentation')
    ax.set_ylabel(f'Correlation ({centrality_type})')
    ax.set_title(f'σ = 0.{sigma_value}')

    return filtered_corr_df


def main():
    fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
    sigma_values = ["00", "01", "03", "05"]
    centrality_types = ['degree', 'betweenness']

    # Create a 4×2 figure
    fig, axes = plt.subplots(4, 2, figsize=(16, 20))

    for i, sigma in enumerate(sigma_values):
        # Load data for this sigma value
        print(f"Processing sigma {sigma}...")
        data = load_sigma_data(sigma, fragmentation_types)

        for j, centrality_type in enumerate(centrality_types):
            print(f"  Processing {centrality_type} centrality...")
            # Run correlation pipeline and plot on the appropriate subplot
            corr_df = run_correlation_pipeline(data, fragmentation_types, centrality_type, sigma, axes[i, j])

            # Save the correlation data
            output_file = f'het_{centrality_type}_correlation_sig{sigma}.csv'
            corr_df.to_csv(output_file, index=False)
            print(f"  Saved {output_file}")

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig('centrality_correlations_combined.svg', dpi=300)
    plt.savefig('centrality_correlations_combined.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()