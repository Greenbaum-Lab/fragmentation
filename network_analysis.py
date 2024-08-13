from random import random

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from funcs import load_data, normalize_steps, calculate_statistics, compute_modularity, calculate_centrality, \
    measure_giant_component, access_networks, access_fst_matrices
from mantel import test
import random


def plot_het_central(data: dict, measure: str):
    fragmentation_types = list(data.keys())
    plt.figure()
    color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette

    for i, frag_type in enumerate(fragmentation_types):
        het = data[frag_type][3]
        central = calculate_centrality(data[frag_type][1], measure=measure)
        merged = pd.merge(het, central, how='outer')
        merged = merged[merged[measure] != 0]

        # merged['modularity'] = np.log10(merged['modularity'])
        # merged['avg'] = np.log10(merged['avg'])
        if measure == 'component':
            sns.regplot(x='component', y='avg', data=merged, fit_reg=True, order=2,
                        truncate=True,
                        scatter_kws={'rasterized': True, 's': 50, 'alpha': 0.01, 'color': color_palette(i)},
                        line_kws={'lw': 2, 'label': frag_type})

            # add a diagonal line
            plt.plot([0.05, 1], [0.05, 1], linestyle='--', color='black',linewidth=1)
            plt.xlabel('Fraction of nodes in the largest component', fontsize=16)

        if measure == 'modularity':
            sns.regplot(x='modularity', y='avg', data=merged, fit_reg=True, order=2,
                        truncate=True,
                        scatter_kws={'rasterized': True, 's': 50, 'alpha': 0.01, 'color': color_palette(i)},
                        line_kws={'lw': 2, 'label': frag_type})

            plt.gca().invert_xaxis()
            # plt.ylim(-0.05, 1.05)
            plt.xlabel('Modularity', fontsize=16)

    plt.ylabel('Heterozygosity', fontsize=16)
    plt.tick_params(axis='both', labelsize=16)

    plt.legend()

    plt.savefig(f'./figs/het_{measure}.svg', format="svg", dpi=300)
    plt.show()


################################
################################ stack plot
def measure_isolated_nodes(network: nx.Graph) -> int:
    """
    Measure the number of isolated nodes in the network.
    :param network: NetworkX graph
    :return: Number of isolated nodes
    """
    isolated_nodes = list(nx.isolates(network))
    return len(isolated_nodes) / len(network)


def measure_components(network: nx.Graph, min_size: int = 4) -> int:
    """
    Measure the number of components with a size greater than or equal to a given threshold,
    excluding the giant component.
    :param network: NetworkX graph
    :param min_size: Minimum size of components to be counted
    :return: Number of nodes in large components excluding the giant component
    """
    largest_component = max(nx.connected_components(network), key=len)

    components = [
        comp for comp in nx.connected_components(network)
        if (comp != largest_component or len(comp) == min_size) and len(comp) >= min_size
    ]

    return sum(len(comp) for comp in components) / len(network)


def measure_waste(network: nx.Graph, max_size: int = 3, min_size: int = 2) -> int:
    """
    Measure the number of components with a size greater than or equal to a given threshold,
    excluding the giant component.
    :param network: NetworkX graph
    :param min_size: Minimum size of components to be counted
    :return: Number of nodes in large components excluding the giant component
    """
    components = [comp for comp in nx.connected_components(network) if min_size <= len(comp) <= max_size]
    num_nodes_in_medium_components = sum(len(comp) for comp in components)
    return num_nodes_in_medium_components / len(network)


def measure_network_metrics(networks: list) -> pd.DataFrame:
    """
    Measure various metrics of the networks and return them as a DataFrame:
    - Size of the giant component
    - Number of isolated nodes
    - Number of components with 4 or more nodes excluding the giant component
    :param networks: List of NetworkX graphs
    :return: DataFrame with metrics for each network
    """
    metrics = []

    for step, network in enumerate(networks):
        giant_component = measure_giant_component(network)
        isolated_nodes = measure_isolated_nodes(network)
        components = measure_components(network)
        waste = measure_waste(network)

        total = giant_component + isolated_nodes + components + waste

        # Round the first three metrics
        giant = round(giant_component / total, 2)
        isolated = round(isolated_nodes / total, 2)
        components = round(components / total, 2)

        # Adjust the last metric so the total sums up to 1
        waste = 1 - giant - isolated - components

        scaled_metrics = {
            "step": step,
            "giant": giant,
            "isolated": isolated,
            "components": components,
            "waste": waste,
        }

        metrics.append(scaled_metrics)

    return pd.DataFrame(metrics)


def measure_network_metrics_replicas(replicas: list) -> pd.DataFrame:
    """
    Measure metrics for a list of lists of networks (replicas) and return a DataFrame
    including a column for the replica index.
    :param replicas: List of lists of NetworkX graphs
    :return: DataFrame with metrics for each network and replica
    """
    all_metrics = []

    for replica_index, networks in enumerate(replicas):
        replica_metrics = measure_network_metrics(networks)
        replica_metrics['replica'] = replica_index
        all_metrics.append(replica_metrics)

    return pd.concat(all_metrics, ignore_index=True)


def calculate_statistics(df):
    """Calculate mean and 95% confidence interval for all columns in the dataframe."""
    result = []

    # Select all columns except 'step' and 'replica'
    columns_to_analyze = df.columns.difference(['step', 'replica'])

    for column in columns_to_analyze:
        mean_values = round(df.groupby('step')[column].mean(), 3)

        # Create a DataFrame for this column's statistics
        column_stats = pd.DataFrame({
            'step': mean_values.index,
            f'{column}': mean_values.values,
            # f'{column}_ci': confidence_interval.values
        })

        result.append(column_stats)

    # Concatenate all column statistics DataFrames along the 'step' index
    result_df = pd.concat(result, axis=1)

    # Remove duplicate 'step' columns
    result_df = result_df.loc[:, ~result_df.columns.duplicated()]
    # result_df['waste'] = 1 - result_df['giant'] - result_df['isolated'] - result_df['components']

    return result_df


def plot_network_stacked_area(df: pd.DataFrame, frag: str):
    """
    Plot the metrics as stacked area charts.
    :param df: DataFrame containing the metrics to plot
    :param frag: Fragmentation type
    """

    # Ensure the DataFrame is sorted by 'step'
    df = df.sort_values(by='step')
    df = normalize_steps(df)

    # Create a new figure and axes with a specific size
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define the columns to plot and the colors to use
    columns = ['waste', 'isolated', 'components', 'giant']
    colors = plt.cm.Dark2.colors[:len(columns)]

    # Prepare the data for the stackplot
    x_values = df['step'].values
    y_values = [df[col].values for col in columns]

    # Create the stackplot
    ax.stackplot(x_values, y_values, labels=columns, colors=colors, alpha=0.8)

    # Set parameters
    ax.set_xlabel('Fragmentation (%)', fontsize=20)
    ax.set_ylabel('Proportion of the network', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=18)  # Increase the size of the tick labels
    ax.set_title(frag)
    plt.ylim(0, 1)
    # plt.legend(loc='upper left')

    plt.savefig(f'./figs/stack_{frag}.jpg')
    plt.show()



###########################################@@#####################
###############################  analysis  #######################

# fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
# # fragmentation_types = ['div']
# data = load_data(fragmentation_types)
#
# ########################plot centrality vs heterozygosity
# plot_het_central(data, measure='component')


##############################plot stacks
# networks = data[fragmentation_types[0]][1]
# matrices = measure_network_metrics_replicas(networks)
# stats = calculate_statistics(matrices)
# plot_network_stacked_area(stats,frag=fragmentation_types[0])

##############################plot centrality vs fragmnetation
# plot_centrality(data,centrality='connectivity')

# stats= calculate_statistics(x)
# print(stats)
# plot_correlation_with_ci(stats, fragmentation_types)


##fst-distance relationship

def get_shortest_path_matrix(net):
    """
    calculate the shortest path length between all pairs of nodes in the network.
    unconnected nodes are marked with inf.
    :return: distance matrix of edges between nodes
    """
    n = nx.number_of_nodes(net)
    distance_matrix = np.full((n, n), np.inf)
    np.fill_diagonal(distance_matrix, 0)  # Distance to self is 0

    # use dijkstra algorithm to calculate the shortest path length
    for source, paths in nx.shortest_path_length(net):
        for target, length in paths.items():
            distance_matrix[source, target] = length

    return distance_matrix


def get_euclidean_matrix(net):
    """
        Calculate the Euclidean distance between all pairs of nodes in the network.
        unconnected nodes are marked with inf.

        :param net:
        :return:
        """
    n = range(nx.number_of_nodes(net))
    # Extract node positions into a numpy array
    positions = nx.get_node_attributes(net, 'pos')
    pos_array = np.array([positions[node] for node in sorted(net.nodes())])
    # Calculate the Euclidean distance matrix using broadcasting
    diff = pos_array[:, np.newaxis, :] - pos_array[np.newaxis, :, :]
    distance_matrix = np.sqrt(np.sum(diff ** 2, axis=-1))

    # If the network is connected, return the distance matrix
    if not nx.is_connected(net):
        pairs = list(combinations(n, 2))
        # insert inf for each pair of nodes that are not connected
        for u, v in pairs:
            if not nx.has_path(net, source=u, target=v):
                distance_matrix[u, v] = np.inf
                distance_matrix[v, u] = np.inf

    return distance_matrix


def random_walk(net, start, end):
    current_node = start
    steps = 0
    while current_node != end:
        neighbors = list(net.neighbors(current_node))
        current_node = random.choice(neighbors)
        steps += 1
    return steps



# def get_random_walk_matrix(net):
#     """
#     calculate the random walk distance between all pairs of nodes in the network.
#     :return: distance matrix of edges between nodes
#     """
#     n = range(nx.number_of_nodes(net))
#     # create a martix of inf in the size of the number of nodes
#     distance_matrix = np.full(((max(n)+1), (max(n))+1), np.inf)
#     np.fill_diagonal(distance_matrix, 0)
#     #get all node pair combinations
#     pairs = list(combinations(n, 2))
#     # insert the random walk distance for each pair of nodes that are connected
#     steps = []
#     for u, v in pairs:
#         if nx.has_path(net, source=u, target=v):
#             steps = []
#
#             for i in range(50):
#                 itreration = random_walk(net, u, v)
#                 steps.append(itreration)
#             res = np.mean(steps)
#             distance_matrix[u, v] = res
#             distance_matrix[v, u] = distance_matrix[u, v]
#     return distance_matrix
#


import concurrent.futures
from itertools import combinations

import concurrent.futures
from itertools import combinations
import numpy as np
import networkx as nx
import random

def random_walk(net, start, end):
    current_node = start
    steps = 0
    while current_node != end:
        neighbors = list(net.neighbors(current_node))
        current_node = random.choice(neighbors)
        steps += 1
    return steps

def compute_random_walk_distance(net, u, v):
    steps = [random_walk(net, u, v) for _ in range(50)]
    return u, v, np.mean(steps)

def get_random_walk_matrix(net, num_workers=None):
    """
    Calculate the random walk distance between all pairs of nodes in the network.
    :param num_workers: Number of threads to use for parallel processing.
    :return: Distance matrix of edges between nodes
    """
    n = range(nx.number_of_nodes(net))
    distance_matrix = np.full((max(n) + 1, max(n) + 1), np.inf)
    np.fill_diagonal(distance_matrix, 0)
    pairs = list(combinations(n, 2))

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(compute_random_walk_distance, net, u, v) for u, v in pairs if nx.has_path(net, u, v)]
        for future in concurrent.futures.as_completed(futures):
            u, v, res = future.result()
            distance_matrix[u, v] = res
            distance_matrix[v, u] = res

    return distance_matrix

def find_connected_components(net):
    """
    Find all connected components in the network.
    use only compoents bigger than 2 nodes.
    :param matrix: distance matrix.
    :return: adjacency matrices of all components.
    """
    indices = []

    components = list(nx.connected_components(net))
    for comp in components:
        comp_nodes = list(comp)
        if len(comp_nodes) > 3:
            indices.append(comp_nodes)
    return indices


def calculate_mantel(net, fst_matrix, dist_type, perms):
    """
    Perform a Mantel test for two distance matrices.

    :param net: NetworkX graph
    :param fst_matrix: Genetic distance matrix
    :param dist_type: Type of distance ('euclidean' or 'path')
    :param perms: Number of permutations for the Mantel test
    :return: Weighted mean of the correlation and p-value by component size
    """
    # calculate the distance matrix based on network connectivity
    if dist_type == 'euclidean':
        distance_matrix = get_euclidean_matrix(net)
    if dist_type == 'path':
        distance_matrix = get_shortest_path_matrix(net)
    if dist_type == 'random':
        distance_matrix = get_random_walk_matrix(net)
    # if the network is connected, calculate the mantel test directly
    if nx.is_connected(net):
        r, p, _ = test(X=distance_matrix, Y=fst_matrix, perms=perms, method='pearson', ignore_nans=True)
        return r, p

    r_values, p_values, weights = [], [], []
    # in case the network is not connected, calculate the mantel test for each component
    for comp in find_connected_components(net):
        comp_dist_matrix = distance_matrix[np.ix_(comp, comp)]
        comp_fst_matrix = fst_matrix[np.ix_(comp, comp)]
        r, p, _ = test(X=comp_dist_matrix, Y=comp_fst_matrix, perms=perms, method='pearson', ignore_nans=True)
        r_values.append(r)
        p_values.append(p)
        if r is np.nan:
            weights.append(0)
        else:
            weights.append(len(comp))
    # Check if r_values or p_values is empty (happens when the network has its last component >3)
    if not r_values or not p_values:
        return None
    # symmetric 3*3 matrices get na in shortest path distance, so drop them in all distance matrices
    if len(r_values) == 1:
        return r_values[0], p_values[0]
    # Mask NaN values that result from symmetric matrices with no variation
    masked_r_values = np.ma.masked_array(r_values, np.isnan(r_values))
    masked_p_values = np.ma.masked_array(p_values, np.isnan(p_values))

    # calculate the weighted mean of the correlation and p-value
    weighted_r = np.ma.average(masked_r_values, weights=weights)
    weighted_p = np.ma.average(masked_p_values, weights=weights)

    return weighted_r, weighted_p


########### analysis of fst-distance

def calculate_mantel_for_process(data, perms, dist_type, replica):
    """"
    calculate mantel correlation and p value for each step along fragmentation.
    :param data: raw data of fragmentation type
    """
    results = []
    networks = access_networks(data)[replica]
    fst_matrices = access_fst_matrices(data)[replica]

    for step, (net, fst) in enumerate(zip(networks, fst_matrices)):
        result = calculate_mantel(net=net, perms=perms, fst_matrix=fst, dist_type=dist_type)
        if result is None:
            break
        r, p = result
        results.append({'step': step, 'r_val': r, 'p_val': p, 'replica': replica})
    cor_data = pd.DataFrame(results)
    return cor_data


def calculate_mantel_replicas(data,perms, dist_type):
    """"
    calculate mantel correlation and p value across fragmentation for all replicas.
    """
    results = []
    networks = access_networks(data)

    for replica in range(len(networks)):
        print(replica)
        cor_data = calculate_mantel_for_process(data, perms, dist_type=dist_type, replica=replica)
        results.append(cor_data)
    cor_data = pd.concat(results)

    return cor_data


def calculate_mantel_all(data, perms, dist_type='euclidean'):
    """"
    calculate mantel correlation and p value across fragmentation for all fragmnetation types.
    """
    results = []
    for frag_type in data.keys():
        print(frag_type)
        cor_data = calculate_mantel_replicas(data[frag_type], perms, dist_type)
        cor_data['fragmentation_type'] = frag_type
        results.append(cor_data)
        cor_data.to_csv(f'./corl_fst_{frag_type}_{dist_type}.csv', index=False)

    cor_data = pd.concat(results)

    # write data as csv
    cor_data.to_csv(f'./corl_fst_{dist_type}_all.csv', index=False)
    return cor_data


def plot_cor_fst(df):
    """
    plot mantel correlation for a single type or replica.
    mini test.
    """
    # df = normalize_steps(df)
    plt.figure(figsize=(10, 6))
    sns.lineplot(x='step', y='r_val', data=df)
    plt.xlabel('Fragmentation (%)')
    plt.ylabel('Correlation')
    plt.show()


def plot_mantel_all(df):
    """
    plot mantel correlation for all fragmentation types.
    """
    df = normalize_steps(df)
    plt.figure(figsize=(10, 6))
    sns.lineplot(x='step', y='r_val', hue='fragmentation_type', data=df)
    plt.xlabel('Fragmentation (%)')
    plt.ylabel('Correlation')
    plt.tick_params(axis='both', labelsize=16)

    plt.savefig(f'./figs/cor_fst_path.svg', format="svg")

    plt.show()


fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
# fragmentation_types = ['rand', 'cor']
data = load_data(fragmentation_types)
# data = data[fragmentation_types[0]]


#### plot single correlation fst-distance
# net = data[1][0][0]
# fst=data[7][0][0]
# distance_matrix = get_random_walk_matrix(net)
# flat_matrix1 = distance_matrix.flatten()
# flat_matrix2 = fst.flatten()
# flat_matrix1 = flat_matrix1[flat_matrix1 != 0]
# flat_matrix2 = flat_matrix2[flat_matrix2 != 0]
# df = pd.DataFrame({'distance': flat_matrix1, 'fst': flat_matrix2})
# ax = sns.regplot(x='distance', y='fst', data=df,
#             fit_reg=True, order=1)
# plt.ylabel('fst', fontsize=20)
# plt.xlabel('random walk', fontsize=20)
# plt.tick_params(axis='both', labelsize=16)
# plt.savefig(f'./figs/random_fst_single.svg', format="svg")
# plt.show()
# print(calculate_mantel(net=net,fst_matrix=fst, dist_type='random',perms=3000))

# df = calculate_mantel_for_process(data, replica=1, dist_type='random', perms=999)
# print(df)
# df = calculate_mantel_replicas(data, dist_type='random', perms=990)
# plot_cor_fst(df)



### plot mantel correletion for all processes
# df = pd.read_csv('./csv/cor_fst_path.csv')
# df = df.replace('--', np.nan)
# df['r_val'] = df['r_val'].astype(float)
# plt.figure(figsize=(10, 6))
# sns.lineplot(x='step', y='r_val', data=df,hue='fragmentation_type')
# plt.xlabel('Fragmentation (%)')
# plt.ylabel('Correlation')
# plt.tick_params(axis='both', labelsize=16)
# plt.savefig(f'./figs/cor_fst_path.svg', format="svg")
# plt.show()

df = calculate_mantel_all(data=data,perms=999, dist_type='random')

# plot_mantel_all(df)
