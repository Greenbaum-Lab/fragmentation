from typing import Dict, List

from matplotlib import pyplot as plt

from data_manipulation.manp_net import het_component, bin_het_component
from funcs import FragmentationResult


def plot_het_component(
    data: Dict[str, FragmentationResult],
    frag_types: List[str] = None,
    n_bins: int = 20,
    output: str = './figs/het_component.svg'
):
    """
    For each fragmentation type, prepare het vs. component data, bin it,
    and plot mean ± SD heterozygosity against fraction in the largest component.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    palette = plt.get_cmap('tab10')

    for i, ft in enumerate(frag_types):
        merged = het_component(data, ft)
        # 2. Bin component fractions and compute mean±SD het
        binned = bin_het_component(merged, n_bins=n_bins)

        color = palette(i)
        ax.scatter(
            binned['component_mid'],
            binned['mean_het'],
            label=ft,
            color=color
        )
        ax.errorbar(
            binned['component_mid'],
            binned['mean_het'],
            yerr=binned['sd_het'],
            fmt='o',
            color=color,
            alpha=0.7
        )

    ax.set_xlabel('Fraction of nodes in largest component', fontsize=16)
    ax.set_ylabel('Heterozygosity', fontsize=16)
    ax.tick_params(labelsize=12)
    plt.tight_layout()
    plt.savefig(output, format='svg')
    plt.show()






