import random
from statistics import mean
import networkx as nx
import matplotlib.pyplot as plt
from joypy import joyplot

import Transformation
from fragmentation import remove_edge_random
from fragmentation import remove_edge_correlated
from fragmentation import remove_edge_distance
from funcs import make_fst_list
from funcs import make_fst_dist
from funcs import make_fst_stat
from funcs import make_het_list
from funcs import make_het_dist
from funcs import make_het_stat
from funcs import frag_random
from funcs import frag_cor
from funcs import frag_dist
from funcs import calculate_centrality
n = 10  # no. of nodes
p = 0.8  # probability to connect nodes
seed = 666
net = nx.erdos_renyi_graph(n=n, p=p)  # create network
n_frag = 5000  # no. of fragmentation steps
pos = nx.spring_layout(net, seed=55)  # set the fixed position for plotting the network
# random.seed(12)  # set random seed

central_rand = calculate_centrality(remove_edge_random(net, n_frag))
central_cor = calculate_centrality(remove_edge_correlated(net, n_frag))
central_dist = calculate_centrality(remove_edge_distance(net, n_frag))

rand = frag_random(net=net, n_frag=n_frag)
cor = frag_cor(net=net, n_frag=n_frag)
dist = frag_dist(net=net, n_frag=n_frag)


# het = make_het_list(rand)
# het_dist = make_het_dist(het)
# het_stat = make_het_stat(het_dist)
# print(het_stat)

# plotting avg and median
plt.plot(rand['step'], rand['avg'], label="average", color="blue")
plt.plot(rand['step'], rand['median'], label="median", color="blue", linestyle='dashed')
plt.plot(cor['step'], cor['avg'], label="average", color="red")
plt.plot(cor['step'], cor['median'], label="median", color="red", linestyle='dashed')
plt.plot(dist['step'], dist['avg'], label="average", color="green")
plt.plot(dist['step'], dist['median'], label="median", color="green", linestyle='dashed')

plt.xlabel('fragmentation process')
plt.ylabel('pairwise fst')
plt.legend()
plt.show()

# plotting avg and median
plt.plot(central_rand['clustering'], rand['avg'], label="average", color="blue")
plt.plot(central_cor['clustering'], cor['avg'], label="average", color="red")
plt.plot(central_dist['clustering'], dist['avg'], label="average", color="green")

plt.xlabel('fragmentation process')
plt.ylabel('pairwise fst')
plt.legend()
plt.show()


#
# # # # plot distribution of Fst values - ridge-lines
# plt.figure()
# joyplot(
#     data=het_dist[['het', 'step']],
#     by='step',
#     colormap=plt.cm.autumn, fade=True,
#     figsize=(12, 8)
# )
# plt.title('pairwise Fst along fragmentation', fontsize=16)
# plt.show()



# het = make_het_list(migration_list)
# het_dens = make_het_dist(het)
# het_stat = make_het_stat(het_dens)


#

# plotting avg and median
# plt.plot(het_stat['step'], het_stat['avg'], label="average")
# plt.plot(het_stat['step'], het_stat['median'], label="median")
#
# plt.xlabel('fragmentation process')
# plt.ylabel('unscaled heterozygosity')
# plt.legend()
# plt.show()

