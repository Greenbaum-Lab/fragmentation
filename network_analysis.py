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
    plot_component_genetics, plot_centrality, plot_degree_distributions, plot_nodes_all, plot_het_central, \
    get_distance_matrix, \
    get_euclidean_matrix, plot_matrix_relationship, plot_node_centrality, plot_network_stacked, measure_network_metrics, \
    calculate_centrality, compute_modularity, measure_giant_component

import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pickle

fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
# fragmentation_types = ['opt']
net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)
print('I have finished loading. Now we start!')


def calculate_centrality(all_nets: list,
                         measure: list = ['clustering', 'degree', 'component',
                                          'modularity', 'transitivity',
                                          'connectivity', 'connect']) -> (
        pd.DataFrame, pd.DataFrame):
    """
    Calculate specified centrality measures of networks over multiple replicates.

    :param all_nets: list of lists of migration networks
    :param measures: list of centrality measures to compute ('clustering', 'path', 'degree' or any combination)

    :return: two dataframes - one with the average values for the specified centrality measures at each step
             and the other with the standard deviations of these values.
    """
    data = []
    for replica, nets in enumerate(all_nets):
        for step, net in enumerate(nets):
            record = {'replica': replica, 'step': step}

            if 'clustering' in measure:
                record['clustering'] = nx.average_clustering(net)

            if 'transitivity' in measure:
                record['transitivity'] = nx.transitivity(net)

            if 'degree' in measure:
                degree = sum(nx.degree_centrality(net).values()) / len(net.nodes)
                record['degree'] = degree

            if 'connect' in measure:
                record['connect'] = nx.average_node_connectivity(net)

            if 'modularity' in measure:
                # partition = community_louvain.best_partition(net, resolution=1)
                # record['modularity'] = community_louvain.modularity(partition, net)
                record['modularity'] = compute_modularity(net)

            if 'connectivity' in measure:
                record['connectivity'] = weighted_algebraic_connectivity(net)

            if 'component' in measure:
                record['component'] = measure_giant_component(net)

            data.append(record)

    df = pd.DataFrame(data)
    return df



def plot_het_central(data: dict, measure: str, save=bool):
    fragmentation_types = list(data.keys())
    plt.figure()
    color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette

    for i, frag_type in enumerate(data.keys()):
        het = data[frag_type][3]
        central = calculate_centrality(data[frag_type][1], measure=measure)
        merged = pd.merge(het, central, how='outer')
        merged = merged[merged[measure] != 0]

        sns.regplot(x='component', y='avg', data=merged, fit_reg=True, order=2,
                    truncate=True, scatter_kws={'s': 50, 'alpha': 0.01, 'color': color_palette(i)},
                    line_kws={'lw': 2, 'label': frag_type})

    if measure == 'modularity':
        plt.gca().invert_xaxis()
        plt.ylim(-0.1, 1.1)

    if measure == 'component':
        max_val = max(plt.gca().get_xlim()[1], plt.gca().get_ylim()[1])
        min_val = min(plt.gca().get_xlim()[1], plt.gca().get_ylim()[1])
        plt.plot([0.05, max_val], [0, max_val], linestyle='--', color='black')

    plt.xlabel('Fraction of nodes in the largest component', fontsize=16)
    plt.ylabel('Heterozygosity', fontsize=16)
    plt.legend()

    if save:
        plt.savefig(f'./figs/het_{measure}.jpg', format="jpg")
    plt.show()

# frag = 'div'
# with open(f'RGG, {frag}_ignore_False.pickle', 'rb') as file:
#     rand = pickle.load(file)


############## plot snapshots of networks along steps of fragmentation
# plot_fragmentation(data)

####plot centrality vs fragmnetation
# plot_centrality(data,centrality='degree')

####plot centrality vs heterozygosity
plot_het_central(data, measure='component', save=True)

####plot degree distribution of fragemtation types
# plot_degree_distributions(data)

###### plot individual nodes
# plot_nodes_all(data)

##### plot heterozygisuty vs. node centrality
# plot_node_centrality(rand,step=0,centrality='degree',log=False,frag=frag)

### plot correlation between centrality and heterozygosity for all processes
# results = compute_correlation_all(data,centrality='degree')
# plot_mean_with_ci(results)


##plot fst-distance relationship
# matrix = rand[7][0][0]
# net = rand[1][0][0]
# # distance_matrix = get_euclidean_matrix(net)
# distance_matrix = get_distance_matrix(net)
# plot_matrix_relationship(distance_matrix=distance_matrix,fst_matrix=matrix)

# make barchart to show the proprtion of strcutures in the network
# x = measure_network_metrics(rand)
# plot_network_stacked(x,frag=frag)


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
