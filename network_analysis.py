import networkx
import random
import statistics
from statistics import mean
import time
import seaborn as sns
from multiprocessing import Pool

import networkx as nx
from community import community_louvain
from joypy import joyplot
from matplotlib import pyplot as plt
import pandas as pd

from processes import find_breaking_point, find_breakink_point_list, remove_edge_random, remove_edge_correlated, \
    remove_edge_distance

import pickle


n = 50  # no. of nodes
p = 0.4  # probability to connect nodes
n_rep = 10

# Record the starting time
start_time = time.time()
#
# # create list off nets
# nets = make_networks(n_nets=n_rep, n_nodes=n, connectivity=p, net_type='ER')
# print("finish nets")

#
# def parallelize_list_comprehension(nets, function):
#     with Pool() as pool:
#         return pool.map(function, nets)
#
#
# breaking_point_rand = parallelize_list_comprehension(nets, remove_edge_random)
# breaking_point_rand = find_breakink_point_list(breaking_point_rand)
# print("finish rand")
#
# file_name= "breaking_point_rand_ER.txt"
# with open(file_name, 'w') as file:
#     for item in breaking_point_rand:
#         file.write(str(item) + '\n')
#
# breaking_point_cor = parallelize_list_comprehension(nets, remove_edge_correlated)
# breaking_point_cor = find_breakink_point_list(breaking_point_cor)
# print("finish cor")
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

#### upload files
# breaking_point_rand = []
# with open('/run/user/1000/gvfs/sftp:host=132.64.61.5,user=lab-heavy/home/lab-heavy/PycharmProjects/fragmentation/breaking_point_rand.txt', 'r') as file:
#     for line in file:
#         item = line.strip()
#         breaking_point_rand.append(item)
#
# breaking_point_cor = []
# with open('/run/user/1000/gvfs/sftp:host=132.64.61.5,user=lab-heavy/home/lab-heavy/PycharmProjects/fragmentation/breaking_point_cor.txt', 'r') as file:
#     for line in file:
#         item = line.strip()
#         breaking_point_cor.append(item)
#
#
# breaking_point_dist = []
# with open('/run/user/1000/gvfs/sftp:host=132.64.61.5,user=lab-heavy/home/lab-heavy/PycharmProjects/fragmentation/breaking_point_dist.txt', 'r') as file:
#     for line in file:
#         item = line.strip()
#         breaking_point_dist.append(item)
#
# breaking_point_rand = [int(element) for element in breaking_point_rand]
# breaking_point_cor = [int(element) for element in breaking_point_cor]
# breaking_point_dist = [int(element) for element in breaking_point_dist]

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
#
# running_time = time.time() - start_time
# print("Running time:", running_time, "seconds")


#######add parameters of percolation
def calculate_centrality_single(net: list) -> pd.DataFrame:
    """
    Calculate the degree of network
    :param net: list of migration networks
    :return: dataframe of degree and clustering for each step
    """
    m = net.copy()
    clustering = list(map(lambda x: nx.average_clustering(x), m))
    betweenness = list(map(lambda x: sum(nx.edge_betweenness_centrality(x).values()) / len(x), m))
    step = range(len(m))
    d = {'step': step, 'clustering': clustering, 'betweenness': betweenness}
    df = pd.DataFrame(data=d)
    return df

def measure_giant_component(network: nx.Graph):
    """
    measure the no. of nodes in the giant component
    :param network:
    :return:
    """
    largest_component = max(nx.connected_components(network), key=len)
    return len(largest_component)


def giant_component_replicates(all_nets: list) -> pd.DataFrame:
    data = []
    for i, networks_list in enumerate(all_nets):
        for step, network in enumerate(networks_list):
            size_giant_component = measure_giant_component(network)
            data.append({'replicate': i, 'step': step, 'avg': size_giant_component})

    df = pd.DataFrame(data)
    return df




# def calculate_centrality(all_nets: list) -> (pd.DataFrame, pd.DataFrame):
#     """
#     Calculate centrality meausres of networks over multiple replicates.
#     :param all_nets: list of lists of migration networks
#     :return: three dataframes - one with the average values for clustering, path and degree centrality at each step
#              and the other with the standard deviations of these values.
#     """
#     data = []
#     for i, nets in enumerate(all_nets):
#         for step, net in enumerate(nets):
#             clustering = nx.average_clustering(net)
#             degree = nx.average_degree_connectivity(net)
#             data.append({'replicate': i, 'step': step, 'clustering': clustering})
#
#     df = pd.DataFrame(data)
#     mean_centrality = df.groupby('step').mean().drop(columns='replicate')
#     std_centrality = df.groupby('step').std().drop(columns='replicate')
#     mean_centrality.columns = ['clustering', 'path', 'degree']
#
#     return mean_centrality, std_centrality
#

def calculate_centrality(all_nets: list, measures: list = ['clustering', 'modularity']) -> (
pd.DataFrame, pd.DataFrame):
    """
    Calculate specified centrality measures of networks over multiple replicates.

    :param all_nets: list of lists of migration networks
    :param measures: list of centrality measures to compute ('clustering', 'path', 'degree' or any combination)

    :return: two dataframes - one with the average values for the specified centrality measures at each step
             and the other with the standard deviations of these values.
    """
    data = []
    for i, nets in enumerate(all_nets):
        for step, net in enumerate(nets):
            record = {'replicate': i, 'step': step}

            if 'clustering' in measures:
                record['clustering'] = nx.transitivity(net)

            if 'modularity' in measures:
                partition = community_louvain.best_partition(net)
                record['modularity'] = community_louvain.modularity(partition, net)

            if 'algebric' in measures:
                record['clustering'] = nx.algebraic_connectivity(net)
            data.append(record)

    df = pd.DataFrame(data)

    # Calculate the means and standard deviations for the specified centrality measures
    mean_centrality = df.groupby('step').mean().drop(columns='replicate')
    std_centrality = df.groupby('step').std().drop(columns='replicate')

    return mean_centrality, std_centrality


# # # Load the tuple using pickle

with open('RGG, rand_ignore_False.pickle', 'rb') as file:
    rand = pickle.load(file)

with open('RGG, cor_ignore_False.pickle', 'rb') as file:
    cor = pickle.load(file)

with open('RGG, dist_ignore_False.pickle', 'rb') as file:
    dist = pickle.load(file)

with open('RGG, int_ignore_False.pickle', 'rb') as file:
    int = pickle.load(file)

with open('RGG, reg_ignore_False.pickle', 'rb') as file:
    reg = pickle.load(file)

with open('RGG, div_ignore_False.pickle', 'rb') as file:
    div = pickle.load(file)

color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette
print("finish load !!!")

all_nets_rand = rand[1]
all_nets_cor = cor[1]
all_nets_dist = dist[1]
all_nets_reg = reg[1]
all_nets_int = int[1]
all_nets_div = div[1]



# #calculate giant component measures
# giant_component = giant_component_replicates(all_nets)
#
# #######heterozygosity
# # Calculate mean and std deviation of GC
# mean_giant_component = giant_component.groupby('step')['avg'].mean()
# confidence_giant_component = giant_component.groupby('step')['avg'].std()
#
# # Calculate mean and std deviation of Het
# mean_rand = rand[3].groupby('step')['avg'].mean()
# confidence_rand = rand[3].groupby('step')['avg'].std()
#
# # Plotting the relationship
# plt.plot(mean_giant_component, label='Mean')
# plt.plot(mean_rand, label='Het')
#
# plt.fill_between(mean_giant_component.index, mean_giant_component - confidence_giant_component, mean_giant_component + confidence_giant_component, alpha=0.2)
# plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
#
# plt.xlabel('Fragmentation')
# plt.ylabel('Mean Number of Nodes in Giant Component')
# plt.title('Mean Number of Nodes in Giant Component vs. Step')
# plt.legend()
# plt.show()

#
# # #calculate giant component measures for all processes
# giant_component_rand = giant_component_replicates(all_nets_rand)
# giant_component_cor = giant_component_replicates(all_nets_cor)
# giant_component_dist = giant_component_replicates(all_nets_dist)
# print("giant comp")
# # Calculate mean and std deviation of GC
# mean_giant_component_rand = giant_component_rand.groupby('step')['avg'].mean()
# confidence_giant_component_rand = giant_component_rand.groupby('step')['avg'].std()
# print("mean of stuff")
#
# mean_giant_component_cor = giant_component_cor.groupby('step')['avg'].mean()
# confidence_giant_component_cor = giant_component_cor.groupby('step')['avg'].std()
#
# mean_giant_component_dist = giant_component_dist.groupby('step')['avg'].mean()
# confidence_giant_component_dist = giant_component_dist.groupby('step')['avg'].std()

# # Calculate mean and std deviation of Het
# mean_rand = rand[3].groupby('step')['avg'].mean()
# confidence_rand = rand[3].groupby('step')['avg'].std()
# mean_cor = cor[3].groupby('step')['avg'].mean()
# confidence_cor = cor[3].groupby('step')['avg'].std()
# mean_dist = dist[3].groupby('step')['avg'].mean()
# confidence_dist = dist[3].groupby('step')['avg'].std()
# mean_int = int[3].groupby('step')['avg'].mean()
# confidence_int = int[3].groupby('step')['avg'].std()
# mean_reg = reg[3].groupby('step')['avg'].mean()
# confidence_reg = reg[3].groupby('step')['avg'].std()
# mean_div = div[3].groupby('step')['avg'].mean()
# confidence_div = div[3].groupby('step')['avg'].std()
# print("mean of het")
#
# mean_centrality_rand, std_centrality_rand = calculate_centrality(all_nets_rand,measures='clustering')
# mean_centrality_cor, std_centrality_cor = calculate_centrality(all_nets_cor,measures='clustering')
# mean_centrality_dist, std_centrality_dist = calculate_centrality(all_nets_dist,measures='clustering')
# mean_centrality_int, std_centrality_int = calculate_centrality(all_nets_int,measures='clustering')
# mean_centrality_reg, std_centrality_reg = calculate_centrality(all_nets_reg,measures='clustering')
# mean_centrality_div, std_centrality_div = calculate_centrality(all_nets_div,measures='clustering')
#
# # plot Clustering
# plt.plot(mean_centrality_rand['clustering'], label='Random')
# plt.plot(mean_centrality_cor['clustering'], label='Correlated')
# plt.plot(mean_centrality_int['clustering'], label='Patchy')
# plt.plot(mean_centrality_reg['clustering'], label='Regressive')
# plt.plot(mean_centrality_div['clustering'], label='Divisive')
# plt.plot(mean_centrality_dist['clustering'], label='Distance-dependent')
# print("centrality!")
#
# plt.fill_between(mean_centrality_rand.index, mean_centrality_rand['clustering'] - std_centrality_rand['clustering'],
#                  mean_centrality_rand['clustering'] + std_centrality_rand['clustering'], alpha=0.2)
# plt.fill_between(mean_centrality_cor.index, mean_centrality_cor['clustering'] - std_centrality_cor['clustering'],
#                  mean_centrality_cor['clustering'] + std_centrality_cor['clustering'], alpha=0.2)
# plt.fill_between(mean_centrality_dist.index, mean_centrality_dist['clustering'] - std_centrality_dist['clustering'],
#                  mean_centrality_dist['clustering'] + std_centrality_dist['clustering'], alpha=0.2)
# plt.fill_between(mean_centrality_int.index, mean_centrality_int['clustering'] - std_centrality_int['clustering'],
#                  mean_centrality_int['clustering'] + std_centrality_int['clustering'], alpha=0.2)
# plt.fill_between(mean_centrality_reg.index, mean_centrality_reg['clustering'] - std_centrality_reg['clustering'],
#                  mean_centrality_reg['clustering'] + std_centrality_reg['clustering'], alpha=0.2)
# plt.fill_between(mean_centrality_div.index, mean_centrality_div['clustering'] - std_centrality_div['clustering'],
#                  mean_centrality_div['clustering'] + std_centrality_div['clustering'], alpha=0.2)
#
# plt.xlabel('Step')
# plt.ylabel('Clustering')
# plt.legend()
# plt.savefig("clust rgg.png", format="png")
# plt.show()





names = ['rand', 'cor', 'dist', 'int', 'reg', 'div']
labels = ['Random', 'Correlated', 'Distance', 'Patchy', 'Regressive', 'Divisive']

# Dictionary to store mean and confidence values
mean_values = {}
confidence_values = {}

# Calculate mean and std deviation for each name
for name in names:
    data = locals()[name][3].groupby('step')['avg']
    mean_values[name] = data.mean()
    confidence_values[name] = data.std()

# Dictionary to store centrality values
mean_centrality = {}
std_centrality = {}

# Calculate centrality for each name
for name, label in zip(names, labels):
    mean_centrality[name], std_centrality[name] = calculate_centrality(locals()['all_nets_' + name], measures='clustering')

# Plot the centralities and fill between the confidence intervals
for name, label in zip(names, labels):
    plt.plot(mean_centrality[name]['clustering'], label=label)
    plt.fill_between(mean_centrality[name].index,
                     mean_centrality[name]['clustering'] - std_centrality[name]['clustering'],
                     mean_centrality[name]['clustering'] + std_centrality[name]['clustering'],
                     alpha=0.2)

# Setting labels and legend
plt.xlabel('Step')
plt.ylabel('Clustering')
plt.legend()
plt.savefig("clust rgg.png", format="png")
plt.show()


####relationship  centrality vs genetics

clust = {}
for name, label in zip(names, labels):
    clust[name] = list(mean_centrality[name]['clustering'])

    plt.plot(clust[name], mean_values[name], label=label)
    plt.fill_between(clust[name], mean_values[name] - confidence_values[name],
                     mean_values[name] + confidence_values[name], alpha=0.2)

plt.gca().invert_xaxis()
plt.xlabel('Clustering')
plt.ylabel('Heterozygosity')
plt.title('Heterozygosity vs. clustering')
plt.legend()
plt.savefig("clust het.png", format="png")
plt.show()


#
# clust_rand = list(mean_centrality_rand['betweenness'])
# clust_cor = list(mean_centrality_cor['betweenness'])
# clust_dist = list(mean_centrality_dist['betweenness'])
#
# # Plotting the relationship
# plt.plot(clust_rand, mean_rand, label='random')
# plt.plot(clust_cor, mean_cor, label='correlated')
# plt.plot(clust_dist, mean_dist, label='distance')
#
# plt.fill_between(clust_rand, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
# plt.fill_between(clust_cor, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
# plt.fill_between(clust_dist, mean_dist - confidence_dist, mean_dist + confidence_dist, alpha=0.2)
#
#
# plt.gca().invert_xaxis()
# plt.xlabel('Betweenness')
# plt.ylabel('Heterozygosity')
# plt.title('Heterozygosity vs. betweenness')
# plt.legend()
# plt.savefig("betweenness het.png", format="png")
# plt.show()

