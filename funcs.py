import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Transformation import m_to_f
from Transformation import m_to_t
from statistics import mean
from statistics import median


def calculate_fst(m: np.ndarray, frag_process, n: int) -> list:
    """
    Calculate Fst of M matrix after each random edge removal
    :param frag_process: type of fragmentation (random, correlated)
    :param m: initial migration network M of networkx
    :return: list of lists of Fst matrices for each step of fragmentation
    """
    fst_list = []
    for i in range(n):
        M = nx.attr_matrix(m)[0]  # take the matrix of the net
        F = m_to_f(M)  # migration to fst function
        frag_process(m=m, n=1)  # use the remove edge function
        fst_list.append(F)  # add another item (fragmentation step) to the list

    return fst_list


def calculate_fst_and_plot(m: np.ndarray, frag_process, n: int) -> list:
    """
    Calculate Fst of M matrix after each random edge removal
    :param frag_process: type of fragmentation (random, correlated)
    :param m: initial migration network M of networkx
    :return: list of lists of Fst matrices for each step of fragmentation
    """
    fst_list = []
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    pos = nx.spring_layout(m, seed=55)
    for i in range(n):
        if i == 0:
            nx.draw_networkx(m, pos=pos, ax=axs[i])
            axs[i].set_title(f"Original network")
        if not nx.is_connected(m):
            nx.draw_networkx(m, pos=pos, ax=axs[1])
            axs[1].set_title(f"Network after {i} edges removed")
            break
        M = nx.attr_matrix(m)[0]  # take the matrix of the net
        F = m_to_f(M)  # migration to fst function
        frag_process(m=m, n=1)  # use the remove edge function
        fst_list.append(F)  # add another item (fragmentation step) to the list
        print(f'I removed {i} edges')
    plt.show()
    return fst_list


def calculate_het(m: np.ndarray, frag_process, n: int) -> list:
    """
    Calculate heterozygosity based on coalescence matrix diagonal of M matrix after each random edge removal
    :param frag_process: type of fragmentation (random, correlated)
    :param m: initial migration network M of networkx
    :return: list of lists of coalescence matrices for each step of fragmentation
    """
    het_list = []
    for i in range(n):
        M = nx.attr_matrix(m)[0]  # take the matrix of the net
        T = m_to_t(M)  # migration to coalescence
        h = np.diag(T)
        h = np.ndarray.tolist(h)
        frag_process(m=m, n=1)  # use the remove edge function
        het_list.append(h)  # add another item (fragmentation step) to the list

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
        avg.append(mean(het_avg))
        het_med = f[f['step'] == i]['het']
        med.append(median(het_med))
        step = range(max(f['step']))
    d = {'step': step, 'avg': avg, 'median': med}
    df = pd.DataFrame(data=d)
    return df


import random

n = 15  # no. of nodes
p = 0.8  # probability to connect nodes
seed = 666
net = nx.erdos_renyi_graph(n, p, seed=seed)  # create network
n_frag = 25  # no. of fragmentation steps
pos = nx.spring_layout(net, seed=55)  # set the fixed position for plotting the network
random.seed(5)


def remove_edge_random(m, n: int) -> list:
    """
    Remove a random edge from net m of type networkx
    :param m: initial migration net m
    :param n: no. of fragmentation steps
    :return: net after edge removal
    """
    migration = []
    for i in range(n):
        edges = list(nx.edges(m))
        edges_to_remove = (random.sample(edges, k=1))  # choose a random edge
        m.remove_edge(*(edges_to_remove[0]))
        migration.append(m.copy())

    return migration


nx.draw(net, pos=pos, with_labels=True)
plt.show()

x = remove_edge_random(m=net, n=n_frag)

nx.draw(x[n_frag - 1], pos=pos, with_labels=True)
plt.show()

#
# def calculate_centrality(m: list) -> pd.DataFrame:
#     """
#      calculate the degree of network
#     :param m: list of migration networks
#     :return: dataframe of degree and clustering for each step
#     """
#     betweenes = []
#     clustering = []
#     for i in range(len(m)):
#         temp = nx.average_clustering(m[i])  # calculate avg clustering
#         clustering.append(temp)
#         temp = nx.betweenness_centrality(m[i])  # calculate betweenes of all nodes
#         temp = mean(list(temp.values()))  # calculate the mean of nodes from dict
#         betweenes.append(temp)
#     step = range(len(m))
#     d = {'step': step, 'clustering': clustering, 'betweenes': betweenes}  # create dict of data
#     df = pd.DataFrame(data=d)
#     print(df)
#     return df



def calculate_centrality(m: list) -> pd.DataFrame:
    """
    Calculate the degree of network
    :param m: list of migration networks
    :return: dataframe of degree and clustering for each step
    """
    clustering = list(map(lambda x: nx.average_clustering(x), m))
    betweenes = list(map(lambda x: mean(list(nx.betweenness_centrality(x).values())), m))
    step = range(len(m))
    d = {'step': step, 'clustering': clustering, 'betweenes': betweenes}
    df = pd.DataFrame(data=d)
    print(df)
    return df


calculate_centrality(x)