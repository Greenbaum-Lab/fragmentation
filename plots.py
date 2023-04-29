import random
from statistics import mean
import networkx as nx
import matplotlib.pyplot as plt
from joypy import joyplot
import seaborn as sns
from fragmentation import remove_edge_random
from fragmentation import remove_edge_correlated
from fragmentation import remove_edges_distance
from funcs import make_fst_list
from funcs import make_fst_dist
from funcs import make_fst_stat
from funcs import make_het_list
from funcs import make_het_dist
from funcs import make_het_stat

n = 5  # no. of nodes
p = 0.8  # probability to connect nodes
seed = 666
net = nx.erdos_renyi_graph(n=n, p=p)  # create network
n_frag = 5000  # no. of fragmentation steps
pos = nx.spring_layout(net, seed=55)  # set the fixed position for plotting the network
# random.seed(12)  # set random seed


# analyze random fragmantation
# def frag_random(net, n_frag:int):
#     xxx = remove_edge_random(migration=net, n=n_frag)
#     fst = make_fst_list(migration_list=xxx)
#     fst_dens = make_fst_dist(fst)
#     fst_stat = make_fst_stat(fst_dens)
#     return fst_stat


# rand = frag_random(net=net, n_frag=10)
rand = remove_edge_random(migration=net, n=50)


fst = make_fst_list(migration_list=rand)
print(len(rand))




# def frag_cor(net, n_frag):
#     xxx = remove_edge_correlated(migration=net, n=n_frag)
#     fst = make_fst_list(migration_list=xxx)
#     fst_dens = make_fst_dist(fst)
#     fst_stat = make_fst_stat(fst_dens)
#     return fst_stat
#
# frag_cor(net=net, n_frag=n_frag)

# # #analyze correlated fragmantation
# migration_list_correlated = remove_edge_correlated(migration=net2, n=50)
# fst_correlated = make_fst_list(migration_list=migration_list_correlated)
# fst_dens_correlated = make_fst_dist(fst_correlated)
# fst_stat_cor = make_fst_stat(fst_dens_correlated)
#
# #analyze distance dependent fragmantation
# migration_list_distance = remove_edges_distance(migration=net3, n=40)
# fst_distance = make_fst_list(migration_list=migration_list_distance)
# fst_dens_distance = make_fst_dist(fst_distance)
# fst_stat_dist = make_fst_stat(fst_dens_distance)
#


# plotting avg and median
# plt.plot(fst_stat_rand['step'], fst_stat_rand['avg'], label="average", color="blue")
# plt.plot(fst_stat_rand['step'], fst_stat_rand['median'], label="median", color="blue", linestyle='dashed')
# plt.plot(fst_stat_cor['step'], fst_stat_cor['avg'], label="average", color="red")
# plt.plot(fst_stat_cor['step'], fst_stat_cor['median'], label="median", color="red", linestyle='dashed')
# plt.plot(fst_stat_dist['step'], fst_stat_dist['avg'], label="average", color="green")
# plt.plot(fst_stat_dist['step'], fst_stat_dist['median'], label="median", color="green", linestyle='dashed')

# plt.xlabel('fragmentation process')
# plt.ylabel('pairwise fst')
# plt.legend()
# plt.show()
#
#

# # # # plot distribution of Fst values - ridge-lines
# plt.figure()
# joyplot(
#     data=fst_dens[['fst', 'step']],
#     by='step',
#     colormap=plt.cm.autumn, fade=True,
#     figsize=(12, 8)
# )
# plt.title('pairwise Fst along fragmentation', fontsize=16)
# plt.show()
#


# het = make_het_list(migration_list)
# het_dens = make_het_dist(het)
# het_stat = make_het_stat(het_dens)


#
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
