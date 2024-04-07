import networkx as nx
import random
import statistics
from statistics import mean

import numpy as np
import seaborn as sns
from multiprocessing import Pool
from community import community_louvain
from joypy import joyplot
from matplotlib import pyplot as plt
import pandas as pd
from scipy import stats
import math
from funcs_analysis import load_data, plot_fragmentation, plot_data, giant_component_replicates, compute_mean_std, \
    plot_component_genetics, plot_centrality, plot_degree_distributions, plot_nodes_all, plot_het_central

import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pickle

fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)
print('I have finished loading. Now we start!')
############## plot snapshots of networks along steps of fragmentation
# plot_fragmentation(data)

####plot centrality vs fragmnetation
# plot_centrality(data,centrality='degree')

####plot centrality vs heterozygosity
plot_het_central(data,measure='transitivity',save=True)

####plot degree distribution of fragemtation types
# plot_degree_distributions(data)

###### plot individual nodes
# plot_nodes_all(data)













################ make distribution of breaking point
# # create list off nets
# nets = make_networks(n_nets=n_rep, n_nodes=n, connectivity=p, net_type='ER')

# def parallelize_list_comprehension(nets, function):
#     with Pool() as pool:
#         return pool.map(function, nets)
#
# breaking_point_rand = parallelize_list_comprehension(nets, remove_edge_random)
# breaking_point_rand = find_breakink_point_list(breaking_point_rand)
#
# file_name= "breaking_point_rand_ER.txt"
# with open(file_name, 'w') as file:
#     for item in breaking_point_rand:
#         file.write(str(item) + '\n')
#
# breaking_point_cor = parallelize_list_comprehension(nets, remove_edge_correlated)
# breaking_point_cor = find_breakink_point_list(breaking_point_cor)
#
# file_name= "breaking_point_cor_ER.txt"
# with open(file_name, 'w') as file:
#     for item in breaking_point_cor:
#         file.write(str(item) + '\n')
#
# breaking_point_dist = parallelize_list_comprehension(nets, remove_edge_distance)
# breaking_point_dist = find_breakink_point_list(breaking_point_dist)
# print("finish dist")
#
# file_name= "breaking_point_dist_ER.txt"
# with open(file_name, 'w') as file:
#     for item in breaking_point_dist:
#         file.write(str(item) + '\n')

######plot
# bins = 100
# sns.histplot(data=breaking_point_rand, bins=bins, kde=True, color='yellow', label='rand')
# sns.histplot(data=breaking_point_cor, bins=bins, kde=True, color='orange', label='cor')
# sns.histplot(data=breaking_point_dist, bins=bins, kde=True, color='grey', label='dist')
#
# plt.xlabel('Value')
# plt.ylabel('Count')
# plt.title('Histogram')
# plt.legend()
# plt.savefig("breaking.png", format="png")
#
# plt.show()



########################### plot  centrality vs genetics
# central = {}
# for name, label in zip(names, labels):
#     central[name] = list(mean_centrality[name][centrality])
#
#     plt.plot(mean_values[name], central[name], label=label)
#     plt.fill_between(central[name], mean_centrality[name][centrality] - std_centrality[name][centrality],
#                      mean_centrality[name][centrality] + std_centrality[name][centrality], alpha=0.2)
#
# plt.gca().invert_xaxis()
# plt.xlabel('Heterozygosity')
# plt.ylabel(centrality)
# plt.title(f'heterozygosity vs. {centrality}')
# # plt.legend()
# plt.savefig("clust het.svg", format="svg")
# plt.show()
#
# for name, label in zip(names, labels):
#     plt.plot(mean_values[name], central[name], label=label)
#     # Assuming mean_values[name] is correctly ordered and corresponds to central[name]
#     # And assuming you have a way to calculate confidence intervals for heterozygosity vs. centrality
#     lower_bound = central[name] - std_deviation  # You need to calculate or retrieve std_deviation
#     upper_bound = central[name] + std_deviation  # Similarly, calculate or retrieve std_deviation
#     plt.fill_between(mean_values[name], lower_bound, upper_bound, alpha=0.2)
