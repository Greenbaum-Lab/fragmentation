import networkx as nx
import numpy as np
import pandas as pd
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


def make_fst_list(migration_list: list) -> list:
    fst_list = []
    for i in range(len(migration_list)):
        M = nx.attr_matrix(migration_list[i])[0]  # take the matrix of the net
        F = m_to_f(M)  # migration to fst function
        fst_list.append(F.copy())  # add another network step to the list

    return fst_list


def make_het_list(migration_list: list) -> list:
    het_list = []
    for i in range(len(migration_list)):
        M = nx.attr_matrix(migration_list[i])[0]  # take the matrix of the net
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


def frag_random(net, n_frag: int):
    """
    run the random fragmentation pipeline to get fst data
    :param net: network
    :param n_frag: no. of steps
    :return: df with fst statistics
    """
    migration = remove_edge_random(net=net, n=n_frag)
    fst = make_fst_list(migration_list=migration)
    fst_dens = make_fst_dist(fst)
    fst_stat = make_fst_stat(fst_dens)

    return fst_stat


def frag_cor(net, n_frag: int):
    """
    run the correlated fragmentation pipeline to get fst data
    :param net: network
    :param n_frag: no. of steps
    :return: df with fst statistics
    """
    migration = remove_edge_correlated(net=net, n=n_frag)
    fst = make_fst_list(migration_list=migration)
    fst_dens = make_fst_dist(fst)
    fst_stat = make_fst_stat(fst_dens)

    return fst_stat


def frag_dist(net, n_frag: int):
    """
    run the distance fragmentation pipeline to get fst data
    :param net: network
    :param n_frag: no. of steps
    :return: df with fst statistics
    """
    migration = remove_edge_distance(net=net, n=n_frag)
    fst = make_fst_list(migration_list=migration)
    fst_dens = make_fst_dist(fst)
    fst_stat = make_fst_stat(fst_dens)

    return fst_stat


def frag_random_giant_comp(net):
    """
    run the random fragmentation pipeline to get fst data
    :param net: network
    :param n_frag: no. of steps
    :return: df with fst statistics
    """
    migration = remove_edge_random_giant_comp(net=net)
    fst = make_fst_list(migration_list=migration)
    fst_dens = make_fst_dist(fst)
    fst_stat = make_fst_stat(fst_dens)

    return fst_stat, fst_dens


def frag_cor_giant_comp(net):
    """
    run the correlated fragmentation pipeline to get fst data
    :param net: network
    :return: df with fst statistics
    """
    migration = remove_edge_correlated_giant_comp(net=net)
    fst = make_fst_list(migration_list=migration)
    fst_dens = make_fst_dist(fst)
    fst_stat = make_fst_stat(fst_dens)

    return fst_stat, fst_dens



def frag_dist_giant_comp(net):
    """
    run the distance-dependent fragmentation pipeline to get fst data
    :param net: network
    :return: df with fst statistics
    """
    migration = remove_edge_correlated_giant_comp(net=net)
    fst = make_fst_list(migration_list=migration)
    fst_dens = make_fst_dist(fst)
    fst_stat = make_fst_stat(fst_dens)

    return fst_stat, fst_dens


