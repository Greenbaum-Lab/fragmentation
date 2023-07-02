import networkx
import random
import statistics
from statistics import mean
import time
import seaborn as sns
from multiprocessing import Pool

import networkx as nx
from joypy import joyplot
from matplotlib import pyplot as plt
import pandas as pd

from funcs3 import make_replicates_new
from processes import find_breaking_point, find_breakink_point_list, remove_edge_random, remove_edge_correlated, \
    remove_edge_distance
from funcs2 import frag_rand, frag_cor, frag_dist, het_rand, het_cor, het_dist, make_networks

# n = 50  # no. of nodes
# p = 0.4  # probability to connect nodes
# n_rep = 200
#
# # Record the starting time
# start_time = time.time()
#
# # create list off nets
# nets = make_networks(n_nets=n_rep, n_nodes=n, connectivity=p, net_type='RGG')
# print("finish nets")
#
#
# def parallelize_list_comprehension(nets, function):
#     with Pool() as pool:
#         return pool.map(function, nets)
#
#
# breaking_point_rand = parallelize_list_comprehension(nets, remove_edge_random)
# breaking_point_rand = find_breakink_point_list(breaking_point_rand)
#
# breaking_point_cor = parallelize_list_comprehension(nets, remove_edge_correlated)
# breaking_point_cor = find_breakink_point_list(breaking_point_cor)
#
# breaking_point_dist = parallelize_list_comprehension(nets, remove_edge_distance)
# breaking_point_dist = find_breakink_point_list(breaking_point_dist)
#
#
# running_time = time.time() - start_time
# print("Running time:", running_time, "seconds")
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

#######add parameters of percolation
def calculate_centrality_single(net: list) -> pd.DataFrame:
    """
    Calculate the degree of network
    :param net: list of migration networks
    :return: dataframe of degree and clustering for each step
    """
    m = net.copy()
    clustering = list(map(lambda x: nx.average_clustering(x), m))
    betweenness = list(map(lambda x: sum(nx.betweenness_centrality(x).values()) / len(x), m))
    step = range(len(m))
    d = {'step': step, 'clustering': clustering, 'betweenness': betweenness}
    df = pd.DataFrame(data=d)
    return df



n = 20  # no. of nodes
p = 0.4  # probability to connect nodes
n_rep = 20

color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette

# # create list off nets
nets = make_networks(n_nets=n_rep, n_nodes=n, connectivity=p, net_type='RGG')
#
# run the pipeline for all fragmentation types
rand = make_replicates_new(nets=nets, frag_type='dist', ignore_isolated=False)

all_nets = rand[1]
def measure_giant_component(network: nx.Graph):
    largest_component = max(nx.connected_components(network), key=len)
    return len(largest_component)


def giant_component_replicates(all_nets: list) -> pd.DataFrame:
    data = []
    for i, networks_list in enumerate(all_nets):
        for step, network in enumerate(networks_list[1:]):  # Starting from the second network
            size_giant_component = measure_giant_component(network)
            data.append({'replicate': i, 'step': step, 'avg': size_giant_component})

    df = pd.DataFrame(data)
    return df

#calculate giant component measures
giant_component = giant_component_replicates(all_nets)

#######heterozygosity
# Calculate mean and std deviation of GC
mean_giant_component = giant_component.groupby('step')['avg'].mean()
confidence_giant_component = giant_component.groupby('step')['avg'].std()

# Calculate mean and std deviation of Het
mean_rand = rand[3].groupby('step')['avg'].mean()
confidence_rand = rand[3].groupby('step')['avg'].std()

# Plotting the relationship
plt.plot(mean_giant_component, label='Mean')
plt.plot(mean_rand, label='Het')

plt.fill_between(mean_giant_component.index, mean_giant_component - confidence_giant_component, mean_giant_component + confidence_giant_component, alpha=0.2)
plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)

plt.xlabel('Fragmentation')
plt.ylabel('Mean Number of Nodes in Giant Component')
plt.title('Mean Number of Nodes in Giant Component vs. Step')
plt.legend()
plt.show()

#######fst
# Calculate mean and std deviation of GC
mean_giant_component = giant_component.groupby('step')['avg'].mean()/20
confidence_giant_component = giant_component.groupby('step')['avg'].std()/20

# Calculate mean and std deviation of Het
mean_rand = rand[5].groupby('step')['avg'].mean()
confidence_rand = rand[5].groupby('step')['avg'].std()

# Plotting the relationship
plt.plot(mean_giant_component, label='component')
plt.plot(mean_rand, label='Het')

plt.fill_between(mean_giant_component.index, mean_giant_component - confidence_giant_component, mean_giant_component + confidence_giant_component, alpha=0.2)
plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)

plt.xlabel('Fragmentation')
plt.ylabel('Mean Number of Nodes in Giant Component')
plt.title('Mean Number of Nodes in Giant Component vs. Step')
plt.legend()
plt.show()




def calculate_centrality(all_nets: list) -> (pd.DataFrame, pd.DataFrame):
    """
    Calculate the average clustering and betweenness centrality of networks over multiple replicates.
    :param all_nets: list of lists of migration networks
    :return: two dataframes - one with the average values for clustering and betweenness centrality at each step
             and the other with the standard deviations of these values.
    """
    data = []
    for i, nets in enumerate(all_nets):
        for step, net in enumerate(nets):
            clustering = nx.average_clustering(net)
            betweenness = sum(nx.betweenness_centrality(net).values()) / len(net)
            data.append({'replicate': i, 'step': step, 'clustering': clustering, 'betweenness': betweenness})

    df = pd.DataFrame(data)
    mean_centrality = df.groupby('step').mean().drop(columns='replicate')
    std_centrality = df.groupby('step').std().drop(columns='replicate')
    mean_centrality.columns = ['clustering', 'betweenness']

    return mean_centrality, std_centrality


mean_centrality, std_centrality = calculate_centrality(all_nets)

# plot Clustering
plt.plot(mean_centrality['clustering'], label='Clustering')
plt.fill_between(mean_centrality.index, mean_centrality['clustering'] - std_centrality['clustering'],
                 mean_centrality['clustering'] + std_centrality['clustering'], alpha=0.2)

plt.xlabel('Step')
plt.ylabel('Clustering')
plt.legend()
plt.show()

# plot Betweenness
plt.plot(mean_centrality['betweenness'], label='Betweenness')
plt.fill_between(mean_centrality.index, mean_centrality['betweenness'] - std_centrality['betweenness'],
                 mean_centrality['betweenness'] + std_centrality['betweenness'], alpha=0.2)

plt.xlabel('Step')
plt.ylabel('Betweenness')
plt.legend()
plt.show()




clust = list(mean_centrality['clustering'])

# Plotting the relationship
plt.plot(clust, mean_rand)
plt.fill_between(clust, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)

plt.xlabel('Clustering')
plt.ylabel('Mean Fst')
plt.title('Fst across clustering')
plt.legend()
plt.show()


clust = list(mean_centrality['betweenness'])

# Plotting the relationship
plt.plot(clust, mean_rand)
plt.fill_between(clust, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)

plt.xlabel('Betweenness')
plt.ylabel('Mean Fst')
plt.title('Fst across betweenness')
plt.legend()
plt.show()
