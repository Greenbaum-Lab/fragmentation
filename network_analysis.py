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
from scipy import stats

from processes import find_breaking_point, find_breakink_point_list, remove_edge_random, remove_edge_correlated, \
    remove_edge_distance

import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pickle


# # create list off nets
# nets = make_networks(n_nets=n_rep, n_nodes=n, connectivity=p, net_type='ER')
# print("finish nets")


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


#######FUNCTIONS
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
    :return: length of giant components
    """
    largest_component = max(nx.connected_components(network), key=len)
    return len(largest_component) / len(network)


def giant_component_replicates(all_nets: list) -> pd.DataFrame:
    """
    measure the no. of nodes in the giant component for a list of networks
    :param all_nets: list of networks
    :return: dataframe
    """
    data = []
    for i, networks_list in enumerate(all_nets):
        for step, network in enumerate(networks_list):
            size_giant_component = measure_giant_component(network)
            data.append({'replicate': i, 'step': step, 'avg': size_giant_component})

    df = pd.DataFrame(data)
    return df


def calculate_centrality(all_nets: list, measures: list = ['clustering', 'degree', 'modularity','connect']) -> (
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
                # record['clustering'] = nx.average_clustering(net)

            if 'degree' in measures:
                degree = sum(nx.degree_centrality(net).values()) / len(net.nodes)
                record['degree'] = degree

            if 'connect' in measures:
                record['connect'] = nx.average_node_connectivity(net)

            if 'modularity' in measures:
                partition = community_louvain.best_partition(net)
                record['modularity'] = community_louvain.modularity(partition, net)

            # if 'algebric' in measures:
            #     record['clustering'] = nx.algebraic_connectivity(net)
            data.append(record)

    df = pd.DataFrame(data)

    # Calculate the means and standard deviations for the specified centrality measures
    mean_centrality = df.groupby('step').mean().drop(columns='replicate')
    std_centrality = df.groupby('step').std().drop(columns='replicate')

    return mean_centrality, std_centrality


def compute_mean_std(data):
    """
    Helper function to compute mean and standard deviation for given data.
    """
    mean = data.groupby('step')['avg'].mean()
    confidence = data.groupby('step')['avg'].std()
    return mean, confidence



# Load the tuple using pickle

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

print("finish load !!!")

################# plot giant component vs heterozygosity
#choose data
frag = rand
# Get giant component measures for all replicates
giant_component_rand = giant_component_replicates(frag[1])

# Calculate mean and std deviation for all replicates
mean_gc_rand, conf_gc_rand = compute_mean_std(giant_component_rand)
mean_het_rand, conf_het_rand = compute_mean_std(frag[3])

# Het - Assuming you want to plot heterozygosity only for rand as shown in your provided code
plt.plot(mean_het_rand, label='Heterozygosity')
plt.plot(mean_gc_rand, label='Giant component')

plt.fill_between(mean_het_rand.index, mean_het_rand - conf_het_rand, mean_het_rand + conf_het_rand, alpha=0.2)
plt.fill_between(mean_gc_rand.index, mean_gc_rand - conf_gc_rand, mean_gc_rand + conf_gc_rand, alpha=0.2)

plt.xlabel('Fragmentation')
plt.ylabel('Fraction of Nodes')
plt.legend()
plt.show()


######################### plot network centrality vs fragmentaion steps
load and store list of networks
names = ['rand', 'cor', 'dist', 'int', 'reg', 'div']
labels = ['Random', 'Correlated', 'Distance', 'Patchy', 'Regressive', 'Divisive']

centrality = 'modularity'

all_nets = {
    'rand': rand[1],  # Extracting the network data, which is the second element of the tuple
    'cor': cor[1],
    'dist': dist[1],
    'int': int[1],
    'reg': reg[1],
    'div': div[1],
}

# Dictionary to store mean and confidence values
mean_values = {}
confidence_values = {}

# Calculate mean and std deviation for each fragmentation
for name in names:
    data = locals()[name][3].groupby('step')['avg']
    mean_values[name] = data.mean()
    confidence_values[name] = data.std()

# Dictionary to store centrality values
mean_centrality = {}
std_centrality = {}

# Calculate centrality for each name
for name, label in zip(names, labels):
    mean_centrality[name], std_centrality[name] = calculate_centrality(all_nets[name], measures=[centrality])

# Plot centrality and fill between the confidence intervals
for name, label in zip(names, labels):
    plt.plot(mean_centrality[name][centrality], label=label)
    plt.fill_between(mean_centrality[name].index,
                     mean_centrality[name][centrality] - std_centrality[name][centrality],
                     mean_centrality[name][centrality] + std_centrality[name][centrality],
                     alpha=0.2)

# Setting labels and legend
plt.xlabel('Step')
plt.ylabel('modularity')
plt.legend()
plt.show()

########################### plot  centrality vs genetics
central = {}
for name, label in zip(names, labels):
    central[name] = list(mean_centrality[name]['modularity'])

    plt.plot(central[name], mean_values[name], label=label)
    plt.fill_between(central[name], mean_values[name] - confidence_values[name],
                     mean_values[name] + confidence_values[name], alpha=0.2)

plt.gca().invert_xaxis()
plt.xlabel('modularity')
plt.ylabel('Heterozygosity')
plt.title('Heterozygosity vs. modularity')
plt.legend()
plt.savefig("clust het.png", format="png")
plt.show()
