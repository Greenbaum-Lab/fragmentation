from typing import Tuple, List, Any

from Transformation import transform_matrix
from processes import remove_edge_random, remove_edge_correlated, remove_edge_distance
import numpy as np
import networkx as nx
import pandas as pd
from statistics import mean, median
from multiprocessing import Pool


def normalize(matrix: np.array) -> np.array:
    # Convert to numpy array in case it's a list
    matrix = np.array(matrix)

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
    fst_list = []
    het_list = []

    for i in range(len(migration_list)):
        M = migration_list[i]
        M = nx.attr_matrix(M)[0]
        T = transform_matrix(M)[0]  # migration to coalescence
        het = np.diag(T)  # take diagonal values (within pop coalesence time=heterozygosity)
        het = np.ndarray.tolist(het)
        het_list.append(het.copy())  # add another network step to the list

        F = transform_matrix(M)[1]  # migration to fst function
        fst_list.append(F.copy())  # add another network step to the list

    return het_list, fst_list


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


#
# def make_fst_stat(f: pd.DataFrame, ignore_isolated: bool) -> pd.DataFrame:
#     """
#      calculate the mean and median fst of each step excluding values of 1
#     :param f: dataframe of fst distribution in each step
#     :return: dataframe of average and median for each step
#     """
#     avg = []
#     med = []
#
#     if ignore_isolated == True:
#         for i in range(max(f['step']) + 1):  # We add 1 to include the max value in the range
#             fst_avg = f[(f['step'] == i) & (f['fst'] != 1)]['fst']
#             fst_med = f[(f['step'] == i) & (f['fst'] != 1)]['fst']
#             # we only append values if the series is not empty
#             if not fst_avg.empty:
#                 avg.append(mean(fst_avg))
#             if not fst_med.empty:
#                 med.append(median(fst_med))
#         step = list(range(max(f['step']) + 1))
#
#     if ignore_isolated == False:
#         for i in range(max(f['step']) + 1):
#             fst_avg = f[f['step'] == i]['fst']
#             avg.append(mean(fst_avg))
#             fst_med = f[f['step'] == i]['fst']
#             med.append(median(fst_med))
#             step = range(max(f['step']))
#     d = {'step': step, 'avg': avg, 'median': med}
#     df = pd.DataFrame(data=d)
#     return df


def make_fst_stat(f: pd.DataFrame, ignore_isolated: bool) -> pd.DataFrame:
    """
     calculate the mean and median fst of each step
     :param f: dataframe of fst distribution in each step
     :param ignore_isolated: if True, ignores rows where 'fst' value is 1
     :return: dataframe of average and median for each step
    """

    if ignore_isolated:
        # Ignore rows where 'fst' equals 1
        f = f[f['fst'] != 1]

    # Get unique steps
    unique_steps = f['step'].unique()

    # Calculate mean and median fst for each step
    avg = f.groupby('step')['fst'].mean().reindex(unique_steps).tolist()
    med = f.groupby('step')['fst'].median().reindex(unique_steps).tolist()

    d = {'step': unique_steps, 'avg': avg, 'median': med}
    df = pd.DataFrame(data=d)
    return df


def make_het_stat(f: pd.DataFrame, ignore_isolated: bool) -> pd.DataFrame:
    """
     calculate the mean and median heterozygosity of each step
    :param f: dataframe of heterozygosity distribution in each step
    :param ignore_isolated: if True, ignores rows where 'het' value is 1
    :return: dataframe of average and median for each step
    """

    if ignore_isolated:
        # Ignore rows where 'het' equals 1
        f = f[f['het'] != 1]

    # Get unique steps
    unique_steps = f['step'].unique()

    # Calculate mean and median heterozygosity for each step
    avg = f.groupby('step')['het'].mean().reindex(unique_steps).tolist()
    med = f.groupby('step')['het'].median().reindex(unique_steps).tolist()

    d = {'step': unique_steps, 'avg': avg, 'median': med}
    df = pd.DataFrame(data=d)
    return df


# create short name to call the desired function
function_mapping = {
    'rand': remove_edge_random,
    'cor': remove_edge_correlated,
    'dist': remove_edge_distance
}


def make_fragmentation(net: nx.Graph, frag_type: str, ignore_isolated: bool) -> tuple:
    """
    run the radom fragmentation pipeline to get heterozygosity data
    :param net: network
    :return: df with heterozygosity statistics
    """
    frag_type = function_mapping[frag_type]

    migration = frag_type(net=net)
    nets_number = len(migration)

    # migration2 = intervals(migration1) # take bins of the process

    genetics = calculate_genetics(migration_list=migration)

    # calculate heterozygosity
    het_dens = make_het_dist(genetics[0])
    het_stat = make_het_stat(het_dens, ignore_isolated=ignore_isolated)

    # calculate fst
    fst_dens = make_fst_dist(genetics[1])
    fst_stat = make_fst_stat(fst_dens, ignore_isolated=ignore_isolated)

    return nets_number, migration, het_dens, het_stat, fst_dens, fst_stat


def make_replicates(nets: list, frag_type: str, ignore_isolated: bool) -> tuple:
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

    for i in range(len(nets)):
        net = make_fragmentation(net=nets[i], frag_type=frag_type, ignore_isolated=ignore_isolated)
        nets_number.append(net[0])
        all_nets.append(net[1])
        het_dens.append(net[2])
        het_stat.append(net[3])
        fst_dens.append(net[4])
        fst_stat.append(net[5])

    # Combine the dataframes into a single dataframe
    nets_number = mean(nets_number)
    het_dens = pd.concat(het_dens)
    het_stat = pd.concat(het_stat)
    fst_dens = pd.concat(fst_dens)
    fst_stat = pd.concat(fst_stat)

    return nets_number, all_nets, het_dens, het_stat, fst_dens, fst_stat


# Function to apply to each network in the list
def apply_make_fragmentation(args):
    net, frag_type, ignore_isolated = args
    return make_fragmentation(net=net, frag_type=frag_type, ignore_isolated=ignore_isolated)

def make_replicates_new(nets: list, frag_type: str, ignore_isolated: bool) -> tuple:
    """
    run multiple iterations of the fragmentation process
    :param ignore_isolated: ignore isolated populations (1)
    :param frag_type: fragmentation type
    :param nets: list of networks
    :return: tuple of all values in dataframes
    """

    # Prepare arguments for the apply_make_fragmentation function
    args = [(net, frag_type, ignore_isolated) for net in nets]

    # Use a pool of workers
    with Pool() as p:
        results = p.map(apply_make_fragmentation, args)

    # Unpack the results
    nets_number, all_nets, het_dens, het_stat, fst_dens, fst_stat = zip(*results)

    # Combine the dataframes into a single dataframe
    nets_number = np.mean(nets_number)
    het_dens = pd.concat(het_dens)
    het_stat = pd.concat(het_stat)
    fst_dens = pd.concat(fst_dens)
    fst_stat = pd.concat(fst_stat)

    return nets_number, all_nets, het_dens, het_stat, fst_dens, fst_stat
