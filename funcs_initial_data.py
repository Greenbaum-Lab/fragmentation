import pickle
from statistics import mean

import numpy as np
import pandas as pd
import networkx as nx
from multiprocessing import Pool

from Transformation import transform_matrix, conservative_from_normal

from processes import (
    remove_edge_random, remove_edge_correlated, remove_edge_distance,
    remove_edge_intrusive, remove_edge_divisive, remove_edge_regressive,
    remove_edge_optimal, remove_edge_optimal_no_update, remove_edge_worst
)



def normalize(matrix: np.array) -> np.array:

    # # Convert to numpy array fron networkx graph
    matrix = nx.attr_matrix(matrix)[0]

    # Calculate row sums
    row_sums = matrix.sum(axis=1)

    # Find the minimum non-zero row sum
    min_row_sum = np.min(row_sums[row_sums > 0])

    # Initialize normalized matrix as a copy of the original
    normalized_matrix = matrix.copy()

    # Get indices of non-zero rows
    non_zero_rows = row_sums > 0

    # Normalize only non-zero rows
    normalized_matrix[non_zero_rows] = matrix[non_zero_rows] / row_sums[non_zero_rows, None] * min_row_sum

    return normalized_matrix


def normalize_list(migration_list: list):
    new_list = list(map(lambda x: normalize(x), migration_list))
    return new_list


def calculate_genetics(migration_list: list) -> tuple:
    """Calculate genetics from migration list."""
    het_list = []
    fst_list = []

    for M in migration_list:
        # M = normalize(M)
        M = nx.attr_matrix(M)[0]
        # Transform matrix to coalescence and fst function
        T, F = transform_matrix(M)

        # Calculate heterozygosity from diagonal values of T matrix
        het = np.diag(T) / len(M)
        het_list.append(het.tolist())

        fst_list.append(F)

    return het_list, fst_list


def make_fst_dist(f: list, ignore: bool = False) -> pd.DataFrame:
    """
    Takes a list of F metrics and returns a DataFrame without diagonal values (zero).
    If ignore_ones is set to True, it will ignore all values of 1.

    Args:
    f : List of FST metrics
    ignore_ones : If True, ignores all values of 1

    Returns:
    DataFrame with a column of all the pairwise FST values and the corresponding fragmentation step.
    """

    fst_values = []
    steps = []

    for i in range(len(f)):
        upper_triangle = np.triu(f[i])

        # Create a boolean mask for non-zero elements
        non_zero_mask = upper_triangle != 0

        # Apply the mask and flatten the array
        upper_triangle_list = upper_triangle[non_zero_mask].flatten()

        if ignore:
            upper_triangle_list = upper_triangle_list[upper_triangle_list != 1]  # ignore all values of 1

        fst_values.extend(upper_triangle_list)  # extend the list with values
        steps.extend([i]*len(upper_triangle_list))  # extend the list with corresponding step

    # Create a DataFrame
    df = pd.DataFrame({
        'step': steps,
        'fst': fst_values
    })

    # Sort by 'step'
    df = df.sort_values(by='step')
    df = df.reset_index(drop=True)

    return df


def make_het_dist(het_list: list, ignore: bool=False) -> pd.DataFrame:
    """
    Takes a list of heterozygosity values and returns a DataFrame.
    If ignore_ones is set to True, it will ignore all values of 1.

    Args:
    het_list : List of heterozygosity vectors
    ignore_ones : If True, ignores all values of 1

    Returns:
    DataFrame with a column of all the heterozygosity values and the corresponding fragmentation step.
    """
    df = pd.DataFrame(het_list)

    if ignore:
        df = df.replace(1, np.nan)  # replace all 1s with NaN

    df = df.stack().rename_axis(('step', 'delete')).reset_index(name='het')
    df = df.drop(columns=['delete'])

    return df



def make_fst_stat(f: pd.DataFrame) -> pd.DataFrame:
    """
     calculate the mean and median fst of each step
     :param f: dataframe of fst distribution in each step
     :param ignore_isolated: if True, ignores rows where 'fst' value is 1
     :return: dataframe of average and median for each step
    """

    # Get unique steps
    unique_steps = f['step'].unique()

    # Calculate mean and median fst for each step
    avg = f.groupby('step')['fst'].mean().reindex(unique_steps).tolist()
    med = f.groupby('step')['fst'].median().reindex(unique_steps).tolist()

    d = {'step': unique_steps, 'avg': avg, 'median': med}
    df = pd.DataFrame(data=d)
    return df


def make_het_stat(f: pd.DataFrame) -> pd.DataFrame:
    """
     calculate the mean and median heterozygosity of each step
    :param f: dataframe of heterozygosity distribution in each step
    :param ignore_isolated: if True, ignores rows where 'het' value is 1
    :return: dataframe of average and median for each step
    """

    # Get unique steps
    unique_steps = f['step'].unique()

    # Calculate mean and median heterozygosity for each step
    avg = f.groupby('step')['het'].mean().reindex(unique_steps).tolist()
    med = f.groupby('step')['het'].median().reindex(unique_steps).tolist()

    d = {'step': unique_steps, 'avg': avg, 'median': med}
    df = pd.DataFrame(data=d)
    return df

from processes import remove_edge_intrusive
from processes import remove_edge_divisive
from processes import remove_edge_regressive
# create short name to call the desired function
function_mapping = {
    'rand': remove_edge_random,
    'cor': remove_edge_correlated,
    'intr': remove_edge_intrusive,
    'reg': remove_edge_regressive,
    'div': remove_edge_divisive,
    'dist': remove_edge_distance,
    'opt': remove_edge_optimal,
    'opt2': remove_edge_optimal_no_update,
    'wrst': remove_edge_worst
}


def make_fragmentation(net: nx.Graph, frag_type: str, ignore: bool, replica: int) -> tuple:
    """
    run the radom fragmentation pipeline to get genetic data
    :param net: network
    :return: df with heterozygosity statistics
    """
    frag_type = function_mapping[frag_type]

    migration = frag_type(net=net)
    nets_number = len(migration)

    # migration2 = intervals(migration1) # take bins of the process

    genetics_coal, genetics_fst = calculate_genetics(migration_list=migration)

    # calculate heterozygosity
    het_dens = make_het_dist(genetics_coal, ignore=ignore)
    het_dens['replica'] = replica

    het_stat = make_het_stat(het_dens)
    het_stat['replica'] = replica

    # calculate fst
    fst_dens = make_fst_dist(genetics_fst, ignore=ignore)
    fst_dens['replica'] = replica
    fst_stat = make_fst_stat(fst_dens)
    fst_stat['replica'] = replica
    return nets_number, migration, het_dens, het_stat, fst_dens, fst_stat, genetics_coal, genetics_fst


def make_replicates(nets: list, frag_type: str, ignore: bool) -> tuple:
    """
    run multiple iterations of the fragmentation process
    :param ignore_isolated: ignore isolated populations (1)
    :param frag_type: fragmentation type
    :param nets: list of networks
    :return: tuple of all values in dataframes
    """
    nets_number = []
    all_nets = []
    het_dens = []
    het_stat = []
    fst_dens = []
    fst_stat = []
    genetics_coal = []
    genetics_fst = []

    for i in range(len(nets)):
        net = make_fragmentation(net=nets[i], frag_type=frag_type, ignore=ignore)
        nets_number.append(net[0])
        all_nets.append(net[1])
        het_dens.append(net[2])
        het_stat.append(net[3])
        fst_dens.append(net[4])
        fst_stat.append(net[5])
        genetics_coal.append(net[6])
        genetics_fst.append(net[7])

    # Combine the dataframes into a single dataframe
    nets_number = mean(nets_number)
    het_dens = pd.concat(het_dens)
    het_stat = pd.concat(het_stat)
    fst_dens = pd.concat(fst_dens)
    fst_stat = pd.concat(fst_stat)


    return nets_number, all_nets, het_dens, het_stat, fst_dens, fst_stat, genetics_coal, genetics_fst


# Function to apply to each network in the list
def apply_make_fragmentation(args):
    net, frag_type, ignore, replica = args
    return make_fragmentation(net=net, frag_type=frag_type, ignore=ignore, replica=replica)


def make_replicates_new(nets: list, frag_type: str, ignore: bool) -> tuple:
    """
    run multiple iterations of the fragmentation process
    :param ignore_isolated: ignore isolated populations (1)
    :param frag_type: fragmentation type
    :param nets: list of networks
    :return: tuple of all values in dataframes
    """

    # Prepare arguments for the apply_make_fragmentation function
    args = [(net, frag_type, ignore, i) for i, net in enumerate(nets)]

    # Use a pool of workers
    with Pool() as p:
        results = p.map(apply_make_fragmentation, args)

    # Unpack the results
    nets_number, all_nets, het_dens, het_stat, fst_dens, fst_stat, genetics_coal, genetics_fst = zip(*results)

    # Combine the dataframes into a single dataframe
    nets_number = np.mean(nets_number)
    het_dens = pd.concat(het_dens)
    het_stat = pd.concat(het_stat)
    fst_dens = pd.concat(fst_dens)
    fst_stat = pd.concat(fst_stat)

    return nets_number, all_nets, het_dens, het_stat, fst_dens, fst_stat, genetics_coal, genetics_fst


def symm_to_assym_net(net: nx.Graph) -> np.array:
    """
    Convert a networkx graph to a numpy array, apply a conservative transformation, and return the new network.
    """
    matrix = nx.attr_matrix(net)[0]
    assymetric_matrix = conservative_from_normal(matrix,mu=1, sigma=0.4, lower=0.2, upper=4)
    new_net = nx.from_numpy_array(assymetric_matrix)
    return new_net


def make_rgg(n_nets: int, n_nodes: int, target_edges: int,assymetric:bool) -> list:
    """
    Create a list of networks with a target number of edges.
    :param n_nets: number of networks
    :param n_nodes: number of nodes
    :param net_type: type of network: ER, RGG, or SF
    :param target_edges: target number of edges
    :return: list of networks
    """
    nets = []
    for _ in range(n_nets):
        while True:
            net = nx.random_geometric_graph(n=n_nodes, radius=0.3)
            print(net.number_of_edges())

            if net.number_of_edges() == target_edges:
                if assymetric:
                    net_assym = symm_to_assym_net(net)
                    nets.append(net_assym)
                else:
                    nets.append(net)
                break
    return nets



# def spatial_sw(N, k, p):
#     G = nx.Graph()
#     positions = {i: np.random.rand(2) for i in range(N)}
#
#     # Add nodes with their positions
#     for i in range(N):
#         G.add_node(i, pos=positions[i])
#
#     # Ensure each node is connected to exactly k nearest neighbors
#     for i in range(N):
#         distances = np.array([np.linalg.norm(positions[i] - positions[j]) for j in range(N)])
#         nearest_neighbors = np.argsort(distances)[1:k + 1]  # Get the k nearest neighbors
#         for j in nearest_neighbors:
#             G.add_edge(i, j)
#
#     # Rewire edges with probability p
#     for i in range(N):
#         neighbors = list(G.neighbors(i))
#         for neighbor in neighbors:
#             if np.random.rand() < p:
#                 non_neighbors = [node for node in range(N) if node not in neighbors and node != i]
#                 if non_neighbors:
#                     new_neighbor = non_neighbors[np.random.randint(len(non_neighbors))]
#                     G.remove_edge(i, neighbor)
#                     G.add_edge(i, new_neighbor)
#
#     return G, positions
# def make_spatial_sw_nets(n_networks, n_nodes, k, p):
#     networks = []
#     for _ in range(n_networks):
#         G, positions = spatial_sw(n_nodes, k, p)
#         networks.append((G, positions))
#         nx.draw_networkx(G, pos=positions)
#         plt.show()
#     return networks



def make_spatial_ER(n, p=0.1):
    net = nx.Graph()
    positions = {i: np.random.rand(2) for i in range(n)}
    net = nx.erdos_renyi_graph(n, p)
    nx.set_node_attributes(net, positions, 'pos')
    return net

def make_spatial_ER_nets(n_nets, n_nodes, p):
    nets = []
    for _ in range(n_nets):
        net = make_spatial_ER(n_nodes, p)
        nets.append(net)
    return nets


def make_spatial_SW(dim=2, p=0.015):
    net = nx.grid_graph(dim=dim,periodic=False)
    mapping = {node: i for i, node in enumerate(net.nodes())}
    net = nx.relabel_nodes(net, mapping)
    pos = {mapping[node]: node for node in mapping}
    nx.set_node_attributes(net, pos, 'pos')

    nodes = list(net.nodes())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if np.random.rand() < p:
                net.add_edge(nodes[i], nodes[j])

    return net

def make_spatial_SW_nets(n_nets, dim=2, p=0.015):
    nets = []
    for _ in range(n_nets):
        net = make_spatial_SW(dim, p)
        nets.append(net)
    return nets

