import random
from statistics import mean
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
from funcs import frag_random_giant_comp
from funcs import frag_cor_giant_comp
from funcs import frag_dist_giant_comp
from funcs import find_breaking_point




n = 10  # no. of nodes
p = 0.8  # probability to connect nodes
seed = 9
# net = nx.erdos_renyi_graph(n=n, p=p)  # create network
net = nx.random_geometric_graph(n=n, radius=0.8)
pos = nx.spring_layout(net, seed=98)  # set the fixed position for plotting the network
random.seed(65)  # set random seed

rand = frag_random_giant_comp(net=net)
cor = frag_cor_giant_comp(net=net)
dist = frag_dist_giant_comp(net=net)

central_rand = calculate_centrality(rand[2])
central_cor = calculate_centrality(cor[2])
central_dist = calculate_centrality(dist[2])




brk_rand = find_breaking_point(rand[2])
brk_cor = find_breaking_point(cor[2])
brk_dist = find_breaking_point(dist[2])



# plotting avg and median

color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette

plt.plot(rand[0]['step'], rand[0]['avg'], label="avg rand", color=color_palette(0))
plt.plot(rand[0]['step'], rand[0]['median'], label="med rand", color=color_palette(0), linestyle='dashed')
plt.plot(cor[0]['step'], cor[0]['avg'], label="med cor", color=color_palette(1))
plt.plot(cor[0]['step'], cor[0]['median'], label="med cor", color=color_palette(1), linestyle='dashed')
plt.plot(dist[0]['step'], dist[0]['avg'], label="avg dist", color=color_palette(2))
plt.plot(dist[0]['step'], dist[0]['median'], label="med dist", color=color_palette(2), linestyle='dashed')

plt.axvline(x=brk_rand, color=color_palette(0),ymax=0.1)
plt.axvline(x=brk_cor, color=color_palette(1),ymax=0.1)
plt.axvline(x=brk_dist, color=color_palette(2),ymax=0.1)
plt.xlabel('fragmentation process')
plt.ylabel('pairwise fst')
plt.legend()
plt.show()




# # merge centrality measures and fst
merged_rand = pd.merge(central_rand, rand[0], on='step')
merged_cor = pd.merge(central_cor, cor[0], on='step')
merged_dist = pd.merge(central_dist, dist[0], on='step')

plt.plot(merged_rand['clustering'], merged_rand['avg'], label="Random", color=color_palette(0))
plt.plot(merged_cor['clustering'], merged_cor['avg'], label="Correlated", color=color_palette(1))
plt.plot(merged_dist['clustering'], merged_dist['avg'], label="Distance", color=color_palette(2))

# plt.xlim(max(merged_rand['clustering']), 0.2)
plt.xlabel("Clustering")
plt.ylabel("Average Fst")
plt.legend()
plt.show()


plt.plot(merged_rand['betweenness'], merged_rand['avg'], label="Random", color=color_palette(0))
plt.plot(merged_cor['betweenness'], merged_cor['avg'], label="Correlated", color=color_palette(1))
plt.plot(merged_dist['betweenness'], merged_dist['avg'], label="Distance", color=color_palette(2))

plt.xlabel("betweenness")
plt.ylabel("Average Fst")
plt.legend()
plt.show()

#
# # # plot distribution of Fst values - ridge-lines
# plt.figure()
# joyplot(
#     data=rand[1][['fst', 'step']],
#     by='step',
#     colormap=plt.cm.autumn, fade=True,
#     figsize=(12, 8)
# )
#
# plt.title('pairwise Fst along random fragmentation', fontsize=16)
# plt.show()

#
# plt.figure()
# joyplot(
#     data=dist[1][['fst', 'step']],
#     by='step',
#     colormap=plt.cm.autumn, fade=True,
#     figsize=(12, 8)
# )
# plt.title('pairwise Fst along distance fragmentation', fontsize=16)
# plt.show()

#
# plt.figure()
# joyplot(
#     data=cor[1][['fst', 'step']],
#     by='step',
#     colormap=plt.cm.autumn, fade=True,
#     figsize=(12, 8)
# )
# plt.title('pairwise Fst along correlated fragmentation', fontsize=16)
# plt.show()
#


# het = make_het_list(migration_list)
# het_dens = make_het_dist(het)
# het_stat = make_het_stat(het_dens)


# plotting avg and median
# plt.plot(het_stat['step'], het_stat['avg'], label="average")
# plt.plot(het_stat['step'], het_stat['median'], label="median")
#
# plt.xlabel('fragmentation process')
# plt.ylabel('unscaled heterozygosity')
# plt.legend()
# plt.show()


# het = make_het_list(rand)
# het_dist = make_het_dist(het)
# het_stat = make_het_stat(het_dist)
# print(het_stat)
