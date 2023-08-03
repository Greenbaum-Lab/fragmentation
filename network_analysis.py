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
    betweenness = list(map(lambda x: sum(nx.betweenness_centrality(x).values()) / len(x), m))
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




# # # Load the tuple using pickle
# with open('rand_include.pickle', 'rb') as file:
#     rand = pickle.load(file)
#
# with open('dist_include.pickle', 'rb') as file:
#     rand = pickle.load(file)
#
with open('rand_include.pickle', 'rb') as file:
    rand = pickle.load(file)

with open('cor_include.pickle', 'rb') as file:
    cor = pickle.load(file)

with open('dist_include.pickle', 'rb') as file:
    dist = pickle.load(file)

color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette
print("finsh whhoppu!!!")


all_nets_rand = rand[1]
all_nets_cor = cor[1]
all_nets_dist = dist[1]

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




# #calculate giant component measures for all processes
giant_component_rand = giant_component_replicates(all_nets_rand)
giant_component_cor = giant_component_replicates(all_nets_cor)
giant_component_dist = giant_component_replicates(all_nets_dist)

# Calculate mean and std deviation of GC
mean_giant_component_rand = giant_component_rand.groupby('step')['avg'].mean()
confidence_giant_component_rand = giant_component_rand.groupby('step')['avg'].std()

mean_giant_component_cor = giant_component_cor.groupby('step')['avg'].mean()
confidence_giant_component_cor = giant_component_cor.groupby('step')['avg'].std()

mean_giant_component_dist = giant_component_dist.groupby('step')['avg'].mean()
confidence_giant_component_dist = giant_component_dist.groupby('step')['avg'].std()

# Calculate mean and std deviation of Het
mean_rand = rand[3].groupby('step')['avg'].mean()
confidence_rand = rand[3].groupby('step')['avg'].std()
mean_cor = cor[3].groupby('step')['avg'].mean()
confidence_cor = cor[3].groupby('step')['avg'].std()
mean_dist = dist[3].groupby('step')['avg'].mean()
confidence_dist = dist[3].groupby('step')['avg'].std()

mean_centrality_rand, std_centrality_rand = calculate_centrality(all_nets_rand)
mean_centrality_cor, std_centrality_cor = calculate_centrality(all_nets_cor)
mean_centrality_dist, std_centrality_dist = calculate_centrality(all_nets_dist)

# plot Clustering
plt.plot(mean_centrality_rand['clustering'], label='random')
plt.plot(mean_centrality_cor['clustering'], label='correlated')
plt.plot(mean_centrality_dist['clustering'], label='distance')

plt.fill_between(mean_centrality_rand.index, mean_centrality_rand['clustering'] - std_centrality_rand['clustering'],
                 mean_centrality_rand['clustering'] + std_centrality_rand['clustering'], alpha=0.2)
plt.fill_between(mean_centrality_cor.index, mean_centrality_cor['clustering'] - std_centrality_cor['clustering'],
                 mean_centrality_cor['clustering'] + std_centrality_cor['clustering'], alpha=0.2)
plt.fill_between(mean_centrality_dist.index, mean_centrality_dist['clustering'] - std_centrality_dist['clustering'],
                 mean_centrality_dist['clustering'] + std_centrality_dist['clustering'], alpha=0.2)

plt.xlabel('Step')
plt.ylabel('Clustering')
plt.legend()
plt.savefig("clustering rgg.png", format="png")
plt.show()

# plot Betweenness
plt.plot(mean_centrality_rand['betweenness'], label='random')
plt.plot(mean_centrality_cor['betweenness'], label='correlated')
plt.plot(mean_centrality_dist['betweenness'], label='distance')

plt.fill_between(mean_centrality_rand.index, mean_centrality_rand['betweenness'] - std_centrality_rand['betweenness'],
                 mean_centrality_rand['betweenness'] + std_centrality_rand['betweenness'], alpha=0.2)
plt.fill_between(mean_centrality_cor.index, mean_centrality_cor['betweenness'] - std_centrality_cor['betweenness'],
                 mean_centrality_cor['betweenness'] + std_centrality_cor['betweenness'], alpha=0.2)
plt.fill_between(mean_centrality_dist.index, mean_centrality_dist['betweenness'] - std_centrality_dist['betweenness'],
                 mean_centrality_dist['betweenness'] + std_centrality_dist['betweenness'], alpha=0.2)

plt.xlabel('Step')
plt.ylabel('Betweenness')
plt.legend()
plt.savefig("betweenness rgg.png", format="png")
plt.show()


####relationship  centrality vs genetics

clust_rand = list(mean_centrality_rand['clustering'])
clust_cor = list(mean_centrality_cor['clustering'])
clust_dist = list(mean_centrality_dist['clustering'])

# Plotting the relationship
plt.plot(clust_rand, mean_rand, label='random')
plt.plot(clust_cor, mean_cor, label='correlated')
plt.plot(clust_dist, mean_dist, label='distance')

plt.fill_between(clust_rand, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
plt.fill_between(clust_cor, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
plt.fill_between(clust_dist, mean_dist - confidence_dist, mean_dist + confidence_dist, alpha=0.2)


plt.gca().invert_xaxis()
plt.xlabel('Clustering')
plt.ylabel('Heterozygosity')
plt.title('Heterozygosity vs. clustering')
plt.legend()
plt.savefig("clustering het.png", format="png")
plt.show()



clust_rand = list(mean_centrality_rand['betweenness'])
clust_cor = list(mean_centrality_cor['betweenness'])
clust_dist = list(mean_centrality_dist['betweenness'])

# Plotting the relationship
plt.plot(clust_rand, mean_rand, label='random')
plt.plot(clust_cor, mean_cor, label='correlated')
plt.plot(clust_dist, mean_dist, label='distance')

plt.fill_between(clust_rand, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
plt.fill_between(clust_cor, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
plt.fill_between(clust_dist, mean_dist - confidence_dist, mean_dist + confidence_dist, alpha=0.2)


plt.gca().invert_xaxis()
plt.xlabel('Betweenness')
plt.ylabel('Heterozygosity')
plt.title('Heterozygosity vs. betweenness')
plt.legend()
plt.savefig("betweenness het.png", format="png")
plt.show()

