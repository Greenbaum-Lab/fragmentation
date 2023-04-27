import random
import networkx as nx
import matplotlib.pyplot as plt
from joypy import joyplot
from funcs import calculate_het
from funcs import make_het_dist
from funcs import make_het_stat
from funcs import calculate_fst
from funcs import make_fst_dist
from funcs import make_fst_stat

from fragmentation import remove_edge_random
from fragmentation import remove_edges_distance

from funcs import calculate_fst_and_plot

n = 10  # no. of nodes
p = 0.5  # probability to connect nodes
# seed = 666
net = nx.erdos_renyi_graph(n, p, seed=seed)  # create network
n_frag = 500  # no. of fragmentation steps
pos = nx.spring_layout(net, seed=55)  # set the fixed position for plotting the network
# random.seed(12)  # set random seed

# 50 nodes- 15 min per fst calculation

nx.draw(net, pos=pos, with_labels=True)
plt.show()

# fst = calculate_fst(m=net, frag_process=remove_edge_random, n=n_frag)
calculate_fst_and_plot(m=net, frag_process=remove_edge_random, n=n_frag)


# fst_dens = make_fst_dist(fst)
# fst_stat = make_fst_stat(fst_dens)



#
# heterozygosity = calculate_het(net, remove_edge_random)
# print(heterozygosity)
# heterozygosity_dens = make_het_dist(heterozygosity)
# print(heterozygosity_dens)
# het_stat = make_het_stat(heterozygosity_dens)
# print(het_stat)
# print()


# plotting avg and median
# plt.plot(fst_stat['step'], fst_stat['avg'], label="average")
# plt.plot(fst_stat['step'], fst_stat['median'], label="median")
#
# plt.xlabel('fragmentation process')
# plt.ylabel('pairwise fst')
# plt.legend()
# plt.show()
#

# # plotting avg and median
# plt.plot(het_stat['step'], het_stat['avg'], label="average")
# plt.plot(het_stat['step'], het_stat['median'], label="median")
#
# plt.xlabel('fragmentation process')
# plt.ylabel('unscaled heterozygosity')
# plt.legend()
# plt.show()
#
# # plot distribution of Fst values - ridge-lines
#
#
# # #
# # plt.figure()
# # joyplot(
# #     data=heterozygosity_dens[['het', 'step']],
# #     by='step',
# #     colormap=plt.cm.autumn, fade=True,
# #     figsize=(12, 8)
# # )
# # plt.title('pairwise Fst along fragmentation', fontsize=16)
# # plt.show()
#
#
# plt.figure()
# joyplot(
#     data=fst_dens[['fst', 'step']],
#     by='step',
#     colormap=plt.cm.autumn, fade=True,
#     figsize=(12, 8)
# )
# plt.title('pairwise Fst along fragmentation', fontsize=20)
# plt.show()
#
