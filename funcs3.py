from typing import Tuple, List, Any

from Transformation import transform_matrix
from processes import remove_edge_random, remove_edge_correlated, remove_edge_distance
import numpy as np
import pandas as pd
from statistics import mean, median
from multiprocessing import Pool
import networkx as nx



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
    fst_list = []
    het_list = []

    for i in range(len(migration_list)):
        M = migration_list[i]
        M = normalize(M)
        T = transform_matrix(M)[0]  # migration to coalescence
        het = np.diag(T)  # take diagonal values (within pop coalesence time=heterozygosity)
        het = het/len(het)
        het = np.ndarray.tolist(het)
        het_list.append(het.copy())  # add another network step to the list

        F = transform_matrix(M)[1]  # migration to fst function
        fst_list.append(F.copy())  # add another network step to the list

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
        F_no_diag = f[i][~np.eye(len(f[i]), dtype=bool)]  # remove diagonals of zero
        if ignore:
            F_no_diag = F_no_diag[F_no_diag != 1]  # ignore all values of 1

        fst_values.extend(F_no_diag)  # extend the list with values
        steps.extend([i]*len(F_no_diag))  # extend the list with corresponding step

    # Create a DataFrame
    df = pd.DataFrame({
        'step': steps,
        'fst': fst_values
    })

    # Sort by 'step'
    df = df.sort_values(by='step')

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
    'int': remove_edge_intrusive,
    'reg': remove_edge_regressive,
    'div': remove_edge_divisive,
    'dist': remove_edge_distance,
}


def make_fragmentation(net: nx.Graph, frag_type: str, ignore: bool) -> tuple:
    """
    run the radom fragmentation pipeline to get genetic data
    :param net: network
    :return: df with heterozygosity statistics
    """
    frag_type = function_mapping[frag_type]

    migration = frag_type(net=net)
    nets_number = len(migration)

    # migration2 = intervals(migration1) # take bins of the process

    genetics = calculate_genetics(migration_list=migration)

    # calculate heterozygosity
    het_dens = make_het_dist(genetics[0], ignore=ignore)
    het_stat = make_het_stat(het_dens)

    # calculate fst
    fst_dens = make_fst_dist(genetics[1], ignore=ignore)
    fst_stat = make_fst_stat(fst_dens)

    return nets_number, migration, het_dens, het_stat, fst_dens, fst_stat


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

    for i in range(len(nets)):
        net = make_fragmentation(net=nets[i], frag_type=frag_type, ignore=ignore)
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
    net, frag_type, ignore = args
    return make_fragmentation(net=net, frag_type=frag_type, ignore=ignore)

def make_replicates_new(nets: list, frag_type: str, ignore: bool) -> tuple:
    """
    run multiple iterations of the fragmentation process
    :param ignore_isolated: ignore isolated populations (1)
    :param frag_type: fragmentation type
    :param nets: list of networks
    :return: tuple of all values in dataframes
    """

    # Prepare arguments for the apply_make_fragmentation function
    args = [(net, frag_type, ignore) for net in nets]

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

def make_networks(n_nets: int, n_nodes: int, net_type) -> list:
    """
    create a list of networks
    :param n_nets: number of networks
    :param n_nodes: number of nodes
    :param connectivity: degree of connectivity
    :param net_type: type of network: ER, RGG, or SF
    :return: list of networks
    """
    nets = []
    for net in range(n_nets):

        if net_type == 'ER':
            net = nx.erdos_renyi_graph(n=n_nodes, p=0.2)
            nets.append(net)
        if net_type == 'RGG':
            net = nx.random_geometric_graph(n=n_nodes, radius=0.3)
            nets.append(net)
        if net_type == 'AB':
            net = nx.barabasi_albert_graph(n=n_nodes, m=5)
            nets.append(net)
        if net_type == 'SW':
            net = nx.watts_strogatz_graph(n=n_nodes,k=9, p=0.1)
            nets.append(net)

    return nets