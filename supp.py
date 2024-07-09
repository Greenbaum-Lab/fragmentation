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
from mantel import test
from scipy.stats import pearsonr
from scipy.stats import norm

from funcs import calculate_statistics, load_data
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



def plot_centrality(data, centrality='modularity'):
    """
    Plots centrality measures and their confidence
    intervals against fragmentation steps.

    :param data: Dictionary containing network data for all fragmentation types.
    :param centrality: The centrality measure to plot.
    """
    names = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
    labels = ['Random', 'Correlated', 'Intrusive', 'Distance', 'Regressive', 'Divisive', 'Optimal']

    plt.figure(figsize=(10, 6))

    for name, label in zip(names, labels):
        # Calculate centrality and its standard deviation using your function
        mean_centrality, std_centrality = calculate_centrality(data[name][1], measure=[centrality])

        steps = mean_centrality.index

        # Plotting the centrality measure for the current fragmentation type
        plt.plot(steps, mean_centrality[centrality], label=label)

        lower_bound = mean_centrality[centrality] - std_centrality[centrality]
        upper_bound = mean_centrality[centrality] + std_centrality[centrality]

        # Plotting the confidence interval as a shaded area
        plt.fill_between(steps, lower_bound, upper_bound, alpha=0.2)

    plt.xlabel('Step', fontsize=22)
    plt.ylabel(centrality.capitalize(), fontsize=18)
    plt.title(f'{centrality.capitalize()} along Fragmentation', fontsize=22)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./figs/{centrality}.jpg')
    plt.show()


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




def get_distance_matrix(net, default_distance=50):
    nodes = list(net.nodes())
    n = len(nodes)
    distance_matrix = np.full((n, n), default_distance)  # Initialize matrix with default distance
    node_index = {node: idx for idx, node in enumerate(nodes)}  # Map nodes to indices

    # Calculate shortest paths using Floyd-Warshall algorithm
    # This considers all path lengths and sets distances for all connected pairs
    path_lengths = dict(nx.all_pairs_dijkstra_path_length(net))

    for i, distances in path_lengths.items():
        for j, dist in distances.items():
            distance_matrix[node_index[i]][node_index[j]] = dist

    return distance_matrix



def plot_matrix_relationship(distance_matrix, fst_matrix, method='pearson', perms=999):
    #calculate the mantel test correlation
    perform_mantel_test(distance_matrix, fst_matrix, perms, method)
    # Flatten the matrices for plotting, ignoring NaN values
    flat_matrix1 = distance_matrix.flatten()
    flat_matrix2 = fst_matrix.flatten()

    valid_indices = ~np.isnan(flat_matrix1) & ~np.isnan(flat_matrix2)
    flat_matrix1 = flat_matrix1[valid_indices]
    flat_matrix2 = flat_matrix2[valid_indices]

    # Create scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(flat_matrix1, flat_matrix2, color='blue', edgecolor='k', alpha=0.7)

    # Add labels and title
    plt.xlabel('Euclidean distance')
    plt.ylabel('Fst')

    # Add a line of best fit
    m, b = np.polyfit(flat_matrix1, flat_matrix2, 1)
    plt.plot(flat_matrix1, m * flat_matrix1 + b, color='red')
    #
    # coeffs = np.polyfit(flat_matrix1, flat_matrix2, 2)  # Quadratic fit
    # p = np.poly1d(coeffs)  # Create polynomial function
    # t = np.linspace(min(flat_matrix1), max(flat_matrix1), 500)
    # plt.plot(t, p(t), color='red')
    plt.savefig(f'./figs/distance_fst.jpg', format="jpg")

    plt.show()


def get_euclidean_matrix(net):
    # Extract node positions into a numpy array
    nodes = list(net.nodes())
    positions = np.array([net.nodes[node]['pos'] for node in nodes])

    # Calculate the difference matrix for each dimension
    diff = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]

    # Compute the Euclidean distance matrix
    distance_matrix = np.linalg.norm(diff, axis=-1)

    return distance_matrix

def perform_mantel_test(distance_matrix, fst_matrix, perms, method,print=bool):
    # Convert all zeros to NaN in both matrices
    # Convert 50 to NaN. 50 is the default value for isolated nodes
    distance_matrix = np.where((distance_matrix == 0) | (distance_matrix == 50), np.nan, distance_matrix)
    fst_matrix = np.where(fst_matrix == 0, np.nan, fst_matrix)

    # Perform Mantel test, expecting a dictionary as a return value
    result = test(fst_matrix, distance_matrix, perms=perms, method=method, ignore_nans=True)

    if print:
        print(f"Correlation: {result[0]}")
        print(f"P-value: {result[1]}")

    return result[0], result[1]




# [1-all networks][replica no.][step number]
# [2-all heterozygosity][replica no.][step number]


