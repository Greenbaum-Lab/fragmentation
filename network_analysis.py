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
    return len(largest_component)/len(network)


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


def calculate_centrality(all_nets: list, measures: list = ['clustering','degree', 'modularity']) -> (
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

            if 'degree' in measures:
                degree = sum(nx.degree_centrality(net).values()) / len(net.nodes)
                record['degree'] = degree

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


def plot_network_realization(data: tuple, step: int, save: bool = False):
    """
    Visualize a network based on betweenness centrality and heterozygosity for a specified step.

    Parameters:
    - data (tuple): A tuple containing relevant data.
                    data[1] is a list of replicates, with each replica being a list of networks.
                    data[2] is a DataFrame with columns 'step' and 'het' (heterozygosity).
    - step (int): The step at which to visualize the network.
    - save (bool, optional): Whether to save the figure. Default is False.

    Returns:
    None
    """

    # Extract relevant network for the specified step
    nets = data[1]
    x = random.choice(range(len(nets)))  # Choose a random replicate from the list
    net = nets[x][step]

    # Calculate betweenness centrality for the selected network
    bet = nx.betweenness_centrality(net)

    # Extract corresponding heterozygosity data
    het_df = data[2]
    het_df = het_df.loc[het_df['step'] == step]
    het_df = het_df.iloc[x*50: (x*50)+50]
    het = {index: value for index, value in enumerate(het_df['het'])}

    # Normalize betweenness centrality values for color mapping
    max_betweenness = max(bet.values())
    normalized_betweenness = {node: value / max_betweenness for node, value in bet.items()}
    colors = plt.cm.Reds([normalized_betweenness[node] for node in net.nodes()])

    # Layout configuration for network visualization
    pos = nx.spring_layout(net, k=0.4, iterations=50)

    # Visualize the network
    nx.draw_networkx_nodes(net, node_color=colors, pos=pos, edgecolors='black',
                           linewidths=1.5, node_size=[v * 300 for v in het.values()])
    nx.draw_networkx_edges(net, pos=pos, edge_color='gray')
    plt.title("Color: Betweenness Centrality | Size: Heterozygosity")

    # Save the figure if specified
    if save:
        plt.savefig("net_betweenness.svg", format='svg')

    # Display the plot
    plt.show()


def centrality_correlation(data: tuple, step: int) -> None:
    """
    Compute and visualize the correlation between betweenness centrality
    and a heterozygosity for nodes in a network at a specified step.

    Parameters:
    - data (tuple): A tuple containing relevant data.
                    data[1] is a list of replicates with each replica is a list of networks
                    data[2] is a DataFrame with columns 'step' and 'het' (heterozygosity)
    - step (int): The step at which to compute the correlation.

    Returns:
    None
    """

    # Extract relevant network step
    nets = data[1]
    nets_step = [net[step - 1] for net in nets]

    # Calculate betweenness centrality for each network in 'result'
    bet = [nx.betweenness_centrality(net) for net in nets_step]

    # Extract relevant heterozygosity values based on the given step
    het_df = data[2]
    het_filtered = het_df.loc[het_df['step'] == step]
    het = het_filtered['het'].tolist()

    # Flatten the 'bet' list of dictionaries into a single list
    flattened_bet = [value for d in bet for value in d.values()]

    # Create a dataframe for correlation analysis
    df = pd.DataFrame({'het': het, 'bet': flattened_bet})

    # Plot regression line
    plot_regression(df)

def plot_regression(df: pd.DataFrame, save=bool) -> None:
    """
    Plot a regression between heterozygosity and betweenness attributes.

    Parameters:
    - df (pd.DataFrame): DataFrame containing 'het' and 'bet' columns.

    Returns:
    None
    """

    plt.figure(figsize=(8, 6))
    sns.regplot(x='bet', y='het', data=df)

    # Compute correlation coefficient and p-value
    r, p = stats.pearsonr(df['bet'], df['het'])
    r_squared = r**2

    # Annotate the plot with r^2 and p-value
    plt.annotate(f'$R^2$ = {r_squared:.3f}', xy=(0.1, 0.9), xycoords='axes fraction', fontsize=14)
    plt.annotate(f'p-value = {p:.5f}', xy=(0.1, 0.85), xycoords='axes fraction', fontsize=14)

    plt.xlabel('Beteeness', fontsize=14)
    plt.ylabel('Heterozygosity', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    if save:
        plt.savefig("cor.svg", format='svg')

    plt.show()




# # Load the tuple using pickle

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




# ################# plot giant component vs heterozygosity
# #choose data
# frag = rand
# # Get giant component measures for all replicates
# giant_component_rand = giant_component_replicates(frag[1])
#
# # Calculate mean and std deviation for all replicates
# mean_gc_rand, conf_gc_rand = compute_mean_std(giant_component_rand)
# mean_het_rand, conf_het_rand = compute_mean_std(frag[3])
#
# # Het - Assuming you want to plot heterozygosity only for rand as shown in your provided code
# plt.plot(mean_het_rand, label='Heterozygosity')
# plt.plot(mean_gc_rand, label='Giant component')
#
# plt.fill_between(mean_het_rand.index, mean_het_rand - conf_het_rand, mean_het_rand + conf_het_rand, alpha=0.2)
# plt.fill_between(mean_gc_rand.index, mean_gc_rand - conf_gc_rand, mean_gc_rand + conf_gc_rand, alpha=0.2)
#
# plt.xlabel('Fragmentation')
# plt.ylabel('Fraction of Nodes')
# plt.legend()
# plt.show()


########################## plot network centrality vs fragmentaion steps
# load and store list of networks
names = ['rand', 'cor', 'dist', 'int', 'reg', 'div']
labels = ['Random', 'Correlated', 'Distance', 'Patchy', 'Regressive', 'Divisive']

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
    mean_centrality[name], std_centrality[name] = calculate_centrality(all_nets[name], measures=['degree'])

# Plot centrality and fill between the confidence intervals
for name, label in zip(names, labels):
    plt.plot(mean_centrality[name]['degree'], label=label)
    plt.fill_between(mean_centrality[name].index,
                     mean_centrality[name]['degree'] - std_centrality[name]['degree'],
                     mean_centrality[name]['degree'] + std_centrality[name]['degree'],
                     alpha=0.2)

# Setting labels and legend
plt.xlabel('Step')
plt.ylabel('Clustering')
plt.legend()
plt.show()


# ########################### plot  centrality vs genetics
# clust = {}
# for name, label in zip(names, labels):
#     clust[name] = list(mean_centrality[name]['clustering'])
#
#     plt.plot(clust[name], mean_values[name], label=label)
#     plt.fill_between(clust[name], mean_values[name] - confidence_values[name],
#                      mean_values[name] + confidence_values[name], alpha=0.2)
#
# plt.gca().invert_xaxis()
# plt.xlabel('Clustering')
# plt.ylabel('Heterozygosity')
# plt.title('Heterozygosity vs. clustering')
# plt.legend()
# plt.savefig("clust het.png", format="png")
# plt.show()


