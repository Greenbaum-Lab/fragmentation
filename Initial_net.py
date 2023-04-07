import array
import random
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from Transformation import m_to_f
from Transformation import m_to_t


def remove_edge(m, n: int):
    """
    Remove a random edge from net m of type networkx
    :param m: initial migration net m
    :param n: no. of fragmentation steps
    :return: net after edge removal
    """

    for i in range(n):
        edges = list(nx.edges(m))
        edges_to_remove = (random.sample(edges, k=1))  # choose a random edge
        m.remove_edge(*(edges_to_remove[0]))

    print(f'I removed {n} edges')

    # nx.draw(m, with_labels=True)
    # plt.show()


n = 10  # no. of nodes
p = 1  # probability to connect nodes
net = nx.erdos_renyi_graph(n, p)  # create ER network
n_frag = 10  # no. of fragmentation steps


def calculate_fst(m) -> list:
    """
    Calculate Fst of M matrix after each random edge removal
    :param m: initial migration network M of networkx
    :return: list of lists of Fst matrices for each step of fragmentation
    """
    fst_list = []
    for i in range(n_frag):
        M = nx.attr_matrix(m)[0]  # take the matrix of the net
        F = m_to_f(M)  # migration to fst function
        remove_edge(m, 1)  # use the remove edge function
        fst_list.append(F)  # add another item (fragmentation step) to the list

    return fst_list


array = calculate_fst(net)
print(array)


def make_fst_data(f: list) -> pd.DataFrame:
    """
    take a list of F metrics and return a dataframe without diagonal values (zero)
    :param f: list of fst metrics
    :return: dataframe with colomm represent each matrix
    """
    fst_data = []
    for i in range(len(f)):
        F_no_diag = f[i][~np.eye(len(f[i]), dtype=bool)]  # remove diagonals of zero and concatante array
        F_no_diag = np.ndarray.tolist(F_no_diag)  # transform to list
        fst_data.append(F_no_diag)  # add another item (fragmentation step) to the list
        df = pd.DataFrame(fst_data)
        df = df.transpose()
    return df


fst_data = make_fst_data(array)
# make dataframe wist one colom of all values-need to remove "a"
fst_data = fst_data.stack().rename_axis(('a', 'step')).reset_index(name='fst')

from joypy import joyplot

plt.figure()

joyplot(
    data=fst_data[['fst', 'step']],
    by='step'
    , colormap=plt.cm.autumn, fade=True,
    figsize=(12, 8)
)
plt.title('pairwise Fst along fragmentation', fontsize=20)
plt.show()
plt.savefig('fst.png')
