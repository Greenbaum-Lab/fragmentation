import math
import pickle
from statistics import mean

import networkx as nx
import numpy as np
import pandas as pd
from infomap import Infomap
from joypy import joyplot
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from scipy.stats import norm

from funcs import calculate_statistics, load_data, calculate_centrality, normalize_steps
from network_analysis import measure_giant_component
from nodes_matrices import calculate_node_centrality


# pd.set_option('display.max_rows', None)


def giant_component_replicates(all_nets: list) -> pd.DataFrame:
    """
    measure the no. of nodes in the giant component for a list of networks
    :param all_nets: list of networks
    :return: dataframe
    """
    data = []
    for i, networks_list in enumerate(all_nets):
        for step, network in enumerate(networks_list):
            size_giant_component = measure_giant_component(network)
            data.append({'replicate': i, 'step': step, 'avg': size_giant_component})

    df = pd.DataFrame(data)
    return df


def plot_component_genetics(data):
    """
    Plot the fraction of nodes in a giant component with heterozygosity
    along fragmentation. Do it for all fragmentation types.
    :param data:
    :return:
    """
    fragmentation_types = list(data.keys())
    # to allow plotting any number of frag types
    num_rows = math.ceil(len(fragmentation_types) / 3)

    fig, axes = plt.subplots(num_rows, 3, figsize=(20, 4 * num_rows))
    axes = axes.flatten()  # Flatten the axes array for easy indexing

    for idx, frag_type in enumerate(fragmentation_types):
        data_frag = data[frag_type]

        giant_component = giant_component_replicates(data_frag[1])
        mean_gc_rand, conf_gc_rand = calculate_statistics(giant_component)
        mean_het_rand, conf_het_rand = calculate_statistics(data_frag[3])

        ax = axes[idx]

        ax.plot(mean_het_rand, label='Heterozygosity')
        ax.plot(mean_gc_rand, label='Giant component')

        ax.fill_between(mean_het_rand.index, mean_het_rand - conf_het_rand,
                        mean_het_rand + conf_het_rand, alpha=0.2)
        ax.fill_between(mean_gc_rand.index, mean_gc_rand - conf_gc_rand,
                        mean_gc_rand + conf_gc_rand, alpha=0.2)

        ax.set_xlabel('Step', fontsize=20)
        ax.set_ylabel('GC/Heterozygosity', fontsize=20)
        ax.set_title(frag_type, fontsize=20, ha='left', loc='left')

        ax.legend()
    plt.savefig('./figs/giant_component.jpg')
    plt.show()


def plot_centrality(data, fragmentation_types, centrality='modularity'):
    """
    Plots centrality measures and their confidence
    intervals against fragmentation steps.

    :param data: Dictionary containing network data for all fragmentation types.
    :param centrality: The centrality measure to plot.
    """

    for frag_type in fragmentation_types:

        data_centrality = calculate_centrality(data[frag_type][1], measure=[centrality])
        data_centrality = calculate_statistics(data_centrality)
        data_centrality = normalize_steps(data_centrality)

        plt.plot(data_centrality.index, data_centrality[f'{centrality}_mean'], label=frag_type)
        lower_bound = data_centrality[f'{centrality}_mean'] - data_centrality[f'{centrality}_ci']
        upper_bound = data_centrality[f'{centrality}_mean'] + data_centrality[f'{centrality}_ci']
        plt.fill_between(data_centrality.index, lower_bound, upper_bound, alpha=0.2)

    plt.xlabel('Fragmentation (%)', fontsize=22)
    plt.ylabel(centrality.capitalize(), fontsize=22)
    if centrality == 'modularity':
        plt.gca().invert_yaxis()
    plt.legend()
    plt.savefig(f'./figs/{centrality}.jpg')
    plt.show()


######## plot centrality along fragmentation

# fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt']
fragmentation_types = ['rand']
net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)

plot_centrality(data, fragmentation_types, centrality='component')


def weighted_algebraic_connectivity(G):
    """Calculate the weighted algebraic connectivity of a graph with disconnected components."""
    components = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    total_weight = sum(len(comp) for comp in components)  # Total weight based on the number of nodes in each component
    weighted_connectivity = sum(len(comp) * nx.algebraic_connectivity(comp) for comp in components) / total_weight
    return weighted_connectivity


def weighted_algebraic_connectivity1(net):
    """Calculate the weighted algebraic connectivity of a graph with disconnected components."""
    if nx.is_connected(net):
        return nx.algebraic_connectivity(net)

    components = list(nx.connected_components(net))
    total_weight = nx.number_of_nodes(net)

    connectivity_sum = 0
    for comp in components:
        comp_size = len(comp)
        comp_subgraph = net.subgraph(comp)
        # Avoid calculation for single-node components as algebraic connectivity would be 0
        if comp_size > 1:
            comp_connectivity = nx.algebraic_connectivity(comp_subgraph)
            connectivity_sum += comp_size * comp_connectivity
        # For a single-node component, you could decide to add or not add to the sum, depending on interpretation
        # In this context, skipping as algebraic connectivity is not defined for single nodes in a meaningful way

    weighted_connectivity = connectivity_sum / total_weight
    return weighted_connectivity


def plot_het_central(data: dict, measure: str, save=bool):
    fragmentation_types = list(data.keys())

    plt.figure()

    for frag_type in fragmentation_types:
        het = calculate_statistics(data[frag_type][3])[0]
        central = calculate_centrality(data[frag_type][1], measure=measure)[0]

        plt.plot(het, central, label=frag_type.capitalize())

    plt.xlabel('Heterozygosity', fontsize=16)
    plt.ylabel('Modularity', fontsize=16)
    plt.legend()
    plt.gca().invert_xaxis()

    if save == True:
        plt.savefig(f'./figs/het_{measure}.jpg', format="jpg")
    plt.show()




# [1-all networks][replica no.][step number]
# [2-all heterozygosity][replica no.][step number]
