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

n = 10  # 100 nodes
p = 0.4
M = nx.erdos_renyi_graph(n, p)


# nx.draw(M, with_labels=True)
# plt.show()
#
# print(M.number_of_edges())
# print(nx.average_clustering(M))
#
# edges = list(M.edges)
# print(edges)
#
# edges_to_remove = (random.sample(edges, k=1))  # choose a random edge
# # edges_to_remove = edges.index((edges_to_remove[0])) # to use index removal if needed
# edges.remove(edges_to_remove[0])  # edges after removal
# print(edges)
# M.remove_edge(*(edges_to_remove[0]))
# print(M.edges)
#
#
#
# edges_to_remove = (random.sample(edges, k=1))  # choose a random edge

def remove_edge(M, f: int):
    """
    remove a random edge from net M
    # :param matrix: initial migration net
    """

    for i in range(f):
        edges = list(nx.edges(M))
        edges_to_remove = (random.sample(edges, k=1))  # choose a random edge
        M.remove_edge(*(edges_to_remove[0]))

        print(f'I removed the edge: {edges_to_remove[0]} \n 'f'These are the remaining edges: {M.edges}')

        nx.draw(M, with_labels=True)
        plt.show()


n = 5  # no. of nodes
p = 0.8  # probability to connect nodes
net = nx.erdos_renyi_graph(n, p)
n_frag = 5
print(net.edges)

remove_edge(net, n_frag)
חטעךלויךלויךלו
