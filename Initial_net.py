# define N=100 and u=0.000005
# M-create ER network-100 nodes and well-connected
# T-calculate the expected coalescence matrix
# F-calculate Fst matrix
# het= coalescence time of a pair of individuals sampled from the same population-distribution and avg
# fst= pairwise fst of population-distribution and avg
# fst= Tt-Ts/Tt= between pop coalescence time-within pop coalescence time
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

        print(f'I removed the edge: {edges_to_remove[0]} \n 'f'These are the remaining edges: {m.edges}')

       # nx.draw(m, with_labels=True)
       # plt.show()


n = 5  # no. of nodes
p = 0.7  # probability to connect nodes
net = nx.erdos_renyi_graph(n, p)  # create ER net
n_frag = 5  # no. of fragmentation steps


# M = nx.attr_matrix(net)[0]  # take the matrix of the net
# F = m_to_f(M)
# print(np.round(F, decimals=2))
#
# ####????when fully connected Fst is 0.0909 and not 0
# F_no_diag = F[~np.eye(len(F), dtype=bool)]  # remove diagonals of zero and concatante array
#
# # plot density plot of Fst
# g = sns.displot(data=F_no_diag, kind="kde")
# g.set_axis_labels("Fst", "Density")
# plt.show()


def calculate_fst(m):
    """
    Calculate Fst of M matrix after each random edge removal
    :param m: initial migration network M of networkx
    :return: list of lists of Fst values for each step of fragmentation
    """
    fst_list = []
    for i in range(n_frag):
        M = nx.attr_matrix(m)[0]  # take the matrix of the net
        F = m_to_f(M)
        F_no_diag = F[~np.eye(len(F), dtype=bool)]  # remove diagonals of zero and concatante array
        F_no_diag = np.ndarray.tolist(F_no_diag)
        remove_edge(m, 1)
        fst_list.append(F_no_diag)

    return fst_list


fst_dist = calculate_fst(net)

df = pd.DataFrame(fst_dist)
df = df.transpose()
print(df)

# plot density plot of Fst
g = sns.displot(data=df, kind="kde")
g.set_axis_labels("Fst", "Density")
plt.show()
