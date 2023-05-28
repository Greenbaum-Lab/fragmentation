import networkx as nx
import numpy as np
import pandas as pd
from joypy import joyplot
from matplotlib import pyplot as plt
from Transformation import m_to_f
from Transformation import m_to_t
from statistics import mean
from statistics import median
from fragmentation import remove_edge_random
from fragmentation import remove_edge_correlated
from fragmentation import remove_edge_distance
from fragmentation import remove_edge_random_giant_comp
from fragmentation import remove_edge_correlated_giant_comp
from fragmentation import remove_edge_distance_giant_comp


def make_fst_list(migration_list: list) -> list:
    fst_list = []
    for i in range(len(migration_list)):
        # M = nx.attr_matrix(migration_list[i])[0]  # take the matrix of the net
        M = migration_list[i]
        F = m_to_f(M)  # migration to fst function
        fst_list.append(F.copy())  # add another network step to the list

    return fst_list


def make_het_list(migration_list: list) -> list:
    het_list = []
    for i in range(len(migration_list)):
        # M = nx.attr_matrix(migration_list[i])[0]  # take the matrix of the net
        M = migration_list[i]
        T = m_to_t(M)  # migration to coalescence
        het = np.diag(T)  # take diagonals values (within pop coalesence time=heterozygosity)
        het = np.ndarray.tolist(het)
        het_list.append(het.copy())  # add another network step to the list

    return het_list


def make_fst_dist(f: list) -> pd.DataFrame:
    """
    take a list of F metrics and return a dataframe without diagonal values (zero)
    :param f: list of fst metrics
    :return: dataframe with a column of all the pairwise fst values
     and the corresponding fragmentation step
    """
    fst_dens = []
    for i in range(len(f)):
        F_no_diag = f[i][~np.eye(len(f[i]), dtype=bool)]  # remove diagonals of zero and concatenate array
        F_no_diag = np.ndarray.tolist(F_no_diag)  # transform to list
        fst_dens.append(F_no_diag)  # add another item (fragmentation step) to the list
    df = pd.DataFrame(fst_dens)
    df = df.transpose()
    df = df.stack().rename_axis(('delete', 'step')).reset_index(name='fst')
    df = df.drop(columns=['delete'])
    df = df.sort_values(by='step')
    return df


def make_het_dist(het_list: list) -> pd.DataFrame:
    """
    take a list of heterozygosity values and return a dataframe
    :param het_list: list of heterozygosity vectors
    :return: dataframe with a column of all the heterozygosity values
    and their corresponding fragmentation step
    """
    df = pd.DataFrame(het_list)
    df = df.stack().rename_axis(('step', 'delete')).reset_index(name='het')
    df = df.drop(columns=['delete'])
    return df


def make_fst_stat(f: pd.DataFrame) -> pd.DataFrame:
    """
     calculate the mean and median fst of each step
    :param f: dataframe of fst distribution in each step
    :return: dataframe of average and median for each step
    """
    avg = []
    med = []
    for i in range(max(f['step'])):
        fst_avg = f[f['step'] == i]['fst']
        avg.append(mean(fst_avg))
        fst_med = f[f['step'] == i]['fst']
        med.append(median(fst_med))
        step = range(max(f['step']))
    d = {'step': step, 'avg': avg, 'median': med}
    df = pd.DataFrame(data=d)
    return df


def make_het_stat(f: pd.DataFrame) -> pd.DataFrame:
    """
     calculate the mean and median heterozygosity of each step
    :param f: dataframe of heterozygosity distribution in each step
    :return: dataframe of average and median for each step
    """
    avg = []
    med = []
    for i in range(max(f['step'])):
        het_avg = f[f['step'] == i]['het']
        avg.append(mean(het_avg.copy()))
        het_med = f[f['step'] == i]['het']
        med.append(median(het_med))
        step = range(max(f['step']))
    d = {'step': step, 'avg': avg, 'median': med}
    df = pd.DataFrame(data=d)
    return df


def calculate_centrality(net: list) -> pd.DataFrame:
    """
    Calculate the degree of network
    :param net: list of migration networks
    :return: dataframe of degree and clustering for each step
    """
    m = net.copy()
    clustering = list(map(lambda x: nx.average_clustering(x), m))
    betweenness = list(map(lambda x: sum(nx.betweenness_centrality(x).values()) / len(x), m))
    step = range(len(m))
    d = {'step': step, 'clustering': clustering, 'betweenness': betweenness}
    df = pd.DataFrame(data=d)
    return df


# def frag_random(net, n_frag: int):
#     """
#     run the random fragmentation pipeline to get fst data
#     :param net: network
#     :param n_frag: no. of steps
#     :return: df with fst statistics
#     """
#     migration = remove_edge_random(net=net, n=n_frag)
#     fst = make_fst_list(migration_list=migration)
#     fst_dens = make_fst_dist(fst)
#     fst_stat = make_fst_stat(fst_dens)
#
#     return fst_stat
#
#
# def frag_cor(net, n_frag: int):
#     """
#     run the correlated fragmentation pipeline to get fst data
#     :param net: network
#     :param n_frag: no. of steps
#     :return: df with fst statistics
#     """
#     migration = remove_edge_correlated(net=net, n=n_frag)
#     fst = make_fst_list(migration_list=migration)
#     fst_dens = make_fst_dist(fst)
#     fst_stat = make_fst_stat(fst_dens)
#
#     return fst_stat
#
#
# def frag_dist(net, n_frag: int):
#     """
#     run the distance fragmentation pipeline to get fst data
#     :param net: network
#     :param n_frag: no. of steps
#     :return: df with fst statistics
#     """
#     migration = remove_edge_distance(net=net, n=n_frag)
#     fst = make_fst_list(migration_list=migration)
#     fst_dens = make_fst_dist(fst)
#     fst_stat = make_fst_stat(fst_dens)
#
#     return fst_stat


def frag_random_giant_comp(net):
    """
    run the random fragmentation pipeline to get fst data
    :param net: network
    :param n_frag: no. of steps
    :return: df with fst statistics
    """
    migration1 = remove_edge_random_giant_comp(net=net)
    migration2 = intervals(migration1)
    migration3 = [nx.attr_matrix(net)[0] for net in migration2]
    migration4 = normalize_list(migration3)
    fst = make_fst_list(migration_list=migration4)
    fst_dens = make_fst_dist(fst)
    fst_stat = make_fst_stat(fst_dens)

    return fst_stat, fst_dens, migration2


def frag_cor_giant_comp(net):
    """
    run the correlated fragmentation pipeline to get fst data
    :param net: network
    :return: df with fst statistics
    """
    migration1 = remove_edge_correlated_giant_comp(net=net)
    migration2 = intervals(migration1)
    migration3 = [nx.attr_matrix(net)[0] for net in migration2]
    migration4 = normalize_list(migration3)
    fst = make_fst_list(migration_list=migration4)
    fst_dens = make_fst_dist(fst)
    fst_stat = make_fst_stat(fst_dens)

    return fst_stat, fst_dens, migration2


def frag_dist_giant_comp(net: nx.Graph):
    """
    run the distance-dependent fragmentation pipeline to get fst data
    :param net: network
    :return: df with fst statistics
    """
    migration1 = remove_edge_distance_giant_comp(net=net)
    migration2 = intervals(migration1)
    migration3 = [nx.attr_matrix(net)[0] for net in migration2]
    migration4 = normalize_list(migration3)
    fst = make_fst_list(migration_list=migration4)
    fst_dens = make_fst_dist(fst)
    fst_stat = make_fst_stat(fst_dens)

    return fst_stat, fst_dens, migration2


def het_rand(net: nx.Graph):
    """
    run the radom fragmentation pipeline to get heterozygosity data
    :param net: network
    :return: df with heterozygosity statistics
    """
    migration1 = remove_edge_random_giant_comp(net=net)
    migration2 = intervals(migration1)
    migration3 = [nx.attr_matrix(net)[0] for net in migration2]
    migration4 = normalize_list(migration3)
    het = make_het_list(migration_list=migration4)
    het_dens = make_het_dist(het)
    het_stat = make_het_stat(het_dens)

    return het_stat, het_dens, migration2


def het_cor(net: nx.Graph):
    """
    run the correlated fragmentation pipeline to get fst data
    :param net: network
    :return: df with heterozygosity statistics
    """
    migration1 = remove_edge_correlated_giant_comp(net=net)
    migration2 = intervals(migration1)
    migration3 = [nx.attr_matrix(net)[0] for net in migration2]
    migration4 = normalize_list(migration3)
    het = make_het_list(migration_list=migration4)
    het_dens = make_het_dist(het)
    het_stat = make_het_stat(het_dens)

    return het_stat, het_dens, migration2


def het_dist(net: nx.Graph):
    """
    run the distance-dependent fragmentation pipeline to get fst data
    :param net: network
    :return: df with fst statistics
    """
    migration1 = remove_edge_distance_giant_comp(net=net)
    migration2 = intervals(migration1)
    migration3 = [nx.attr_matrix(net)[0] for net in migration2]
    migration4 = normalize_list(migration3)
    het = make_het_list(migration_list=migration4)
    het_dens = make_het_dist(het)
    het_stat = make_het_stat(het_dens)

    return het_stat, het_dens, migration2


def normalize(array: np.array) -> np.array:
    """
    normalize the migration matrix so that all row sums will be the same
    the sum will be equal to the sum of the row with the lowest sum
    :param array:  migration network with 0,1
    :return: scaled migration network
    """
    # Find the row with the lowest sum
    min_sum_row = np.argmin(np.sum(array, axis=1))

    # Calculate the desired row sum
    desired_sum = np.sum(array[min_sum_row])

    # Calculate the current row sums
    row_sums = np.sum(array, axis=1)

    # Calculate the scaling factors needed for each row
    scaling_factors = desired_sum / row_sums

    # Replace the values in the array with the scaled values
    scaled_arr = array * scaling_factors[:, np.newaxis]

    return scaled_arr


def normalize_list(migration_list: list):
    new_list = list(map(lambda x: normalize(x), migration_list))
    return new_list


def intervals(lst):
    if len(lst) <= 20:
        return lst

    interval = max((len(lst) - 1) // 19, 1)
    return lst[:19 * interval:interval] + [lst[-1]]


def find_breaking_point(lst):
    """
    find the index of the list where the network is no longer connected
    """
    first_element = lst[0]

    for i, element in enumerate(lst):
        if len(element) < len(first_element):
            return i

    return -1  # Return -1 if no element is found


def make_networks(n_nets: int, n_nodes: int, connectivity: float, net_type) -> list:
    nets = []
    for net in range(n_nets):
        if net_type == 'ER':
            net = nx.erdos_renyi_graph(n=n_nodes, p=connectivity)
            nets.append(net)
        if net_type == 'RGG':
            net = nx.random_geometric_graph(n=n_nodes, radius=connectivity)
            nets.append(net)
        if net_type == 'SF':
            net = nx.barabasi_albert_graph(n=n_nodes, m=2)
            nets.append(net)
    return nets


nets=make_networks(n_nets=3, n_nodes=20, connectivity=0.8, net_type='ER')


def make_iterations(nets: list, fragmentation) -> pd.DataFrame:
    """
    run multiple iterations of the fragmentation process
    :param nets: list of networks
    :param fragmentation: fragmentation process
    :return: dataframe of avg fst for each net
    """
    all_stat = []
    all_dens = []

    # Get the corresponding function based on the nickname
    selected_frag = function_mapping[fragmentation]

    # calculate fst for each network in the list
    for i in range(len(nets)):
        net = selected_frag(net=nets[i])
        all_stat.append(net[0])
        all_dens.append(net[1])

    # Combine the dataframes into a single dataframe
    combined_stat = pd.concat(all_stat)
    combined_dens = pd.concat(all_dens)
    return combined_stat, combined_dens


#

# create short name to call the desired function
function_mapping = {
    'rand': frag_random_giant_comp,
    'dist': frag_dist_giant_comp,
    'cor': frag_cor_giant_comp
}

rand = make_iterations(nets, fragmentation='rand')
# cor = make_iterations(nets, fragmentation='cor')
# dist = make_iterations(nets, fragmentation='dist')


# # plot distribution of Fst values - ridge-lines
plt.figure()
joyplot(
    data=rand[1][['fst', 'step']],
    by='step', ylim=0, overlap=0.5,
    colormap=plt.cm.autumn, fade=True,
    figsize=(12, 8)
)

plt.title('pairwise Fst along random fragmentation', fontsize=16)
plt.show()

# # Calculate the mean and median values over the 'step' column
# mean_rand = rand[0].groupby('step')['avg'].mean()
# mean_cor = cor[0].groupby('step')['avg'].mean()
# mean_dist = dist[0].groupby('step')['avg'].mean()
#
# # Calculate the confidence interval
# confidence_rand = rand[0].groupby('step')['avg'].std()
# confidence_cor = cor[0].groupby('step')['avg'].std()
# confidence_dist = dist[0].groupby('step')['avg'].std()
#
# # Plotting the line graph with mean and median values
# plt.plot(mean_rand, label='Rand')
# plt.plot(mean_cor, label='Cor')
# plt.plot(mean_dist, label='Dist')
#
# # Adding the confidence interval
# plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
# plt.fill_between(mean_cor.index, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
# plt.fill_between(mean_dist.index, mean_dist - confidence_dist, mean_dist + confidence_dist, alpha=0.2)
#
# # Add labels and legend
# plt.xlabel('Fragmentation step')
# plt.ylabel('Pairwise Fst')
# plt.legend()
#
# # Display the plot
# plt.show()
