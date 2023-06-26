from typing import Tuple, List, Any

from Transformation import transform_matrix
from processes import remove_edge_random
import numpy as np
import networkx as nx


net = nx.random_geometric_graph(20,0.5)
mifration = remove_edge_random(net)


def calculate_genetics(migration_list: list) -> tuple[list[Any], list[Any]]:
    fst_list = []
    het_list = []

    for i in range(len(migration_list)):
        M = migration_list[i]
        T = transform_matrix(M)[0]  # migration to coalescence
        het = np.diag(T)  # take diagonal values (within pop coalesence time=heterozygosity)
        het = np.ndarray.tolist(het)
        het_list.append(het.copy())  # add another network step to the list

        F = transform_matrix(M)[1]  # migration to fst function
        fst_list.append(F.copy())  # add another network step to the list

    return fst_list, het_list


test = calculate_genetics(mifration)



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
