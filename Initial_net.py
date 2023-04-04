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

        nx.draw(m, with_labels=True)
        plt.show()


n = 20  # no. of nodes
p = 0.8  # probability to connect nodes
net = nx.erdos_renyi_graph(n, p)  # create ER net
n_frag = 5  # no. of fragmentation steps


# remove_edge(net, n_frag)
M = nx.attr_matrix(net)[0]
print(M)

# F_1 = m_to_f(M)
# print(np.round(F_1, decimals=2))
# nx.draw(net, with_labels=True)
# plt.show()


M = np.concatenate(M)
g = sns.displot(data=M, kind="kde")
g.set_axis_labels("Fst", "Density")

plt.show()



# data.speeding.plot.density(color='green')
# plt.title('Density plot for Speeding')
# plt.show()
# # loading the dataset
# # from seaborn library
#
# # viewing the dataset
# print(data.head(4))
#
# # plotting the density plot
# # for 'speeding' attribute
# # using plot.density()
# data.speeding.plot.density(color='green')
# plt.title('Density plot for Speeding')
# plt.show()
