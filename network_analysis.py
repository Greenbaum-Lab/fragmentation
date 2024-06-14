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
    get_euclidean_matrix, plot_matrix_relationship, calculate_centrality, compute_modularity

import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pickle

fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
frag = 'dist'

fragmentation_types = [frag]
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

####plot centrality vs heterozygosity
# plot_het_central(data, measure='component', save=True)


################################
################################ stack plot
def measure_giant_component(network: nx.Graph, min_size: int = 4):
    """
    measure the no. of nodes in the giant component
    :param network:
    :return: length of giant components
    """
    largest_component = max(nx.connected_components(network), key=len)
    if len(largest_component) <= min_size:
        return 0
    return len(largest_component) / len(network)


def measure_isolated_nodes(network: nx.Graph) -> int:
    """
    Measure the number of isolated nodes in the network.
    :param network: NetworkX graph
    :return: Number of isolated nodes
    """
    isolated_nodes = list(nx.isolates(network))
    return len(isolated_nodes) / len(network)

def measure_components(network: nx.Graph, min_size: int = 4) -> int:
    """
    Measure the number of components with a size greater than or equal to a given threshold,
    excluding the giant component.
    :param network: NetworkX graph
    :param min_size: Minimum size of components to be counted
    :return: Number of nodes in large components excluding the giant component
    """
    largest_component = max(nx.connected_components(network), key=len)

    components = [
        comp for comp in nx.connected_components(network)
        if (comp != largest_component or len(comp) == min_size) and len(comp) >= min_size
    ]

    return sum(len(comp) for comp in components) / len(network)


def measure_waste(network: nx.Graph, max_size: int = 3, min_size: int = 2) -> int:
    """
    Measure the number of components with a size greater than or equal to a given threshold,
    excluding the giant component.
    :param network: NetworkX graph
    :param min_size: Minimum size of components to be counted
    :return: Number of nodes in large components excluding the giant component
    """
    components = [comp for comp in nx.connected_components(network) if min_size <= len(comp) <= max_size]
    num_nodes_in_medium_components = sum(len(comp) for comp in components)
    return num_nodes_in_medium_components / len(network)


def measure_network_metrics(networks: list) -> pd.DataFrame:
    """
    Measure various metrics of the networks and return them as a DataFrame:
    - Size of the giant component
    - Number of isolated nodes
    - Number of components with 4 or more nodes excluding the giant component
    :param networks: List of NetworkX graphs
    :return: DataFrame with metrics for each network
    """
    metrics = []

    for step, network in enumerate(networks):
        giant_component= measure_giant_component(network)
        isolated_nodes = measure_isolated_nodes(network)
        components = measure_components(network)
        waste = measure_waste(network)

        scaled_metrics = {
            "step": step,
            "giant": giant_component ,
            "isolated": isolated_nodes,
            "components": components,
            "waste": waste,
        }
        total = giant_component+isolated_nodes+components+waste

        if total != 1:
            print(f"sum of values higher than it should be in {step}")

        metrics.append(scaled_metrics)

    return pd.DataFrame(metrics)



def measure_network_metrics_replicas(replicas: list) -> pd.DataFrame:
    """
    Measure metrics for a list of lists of networks (replicas) and return a DataFrame
    including a column for the replica index.
    :param replicas: List of lists of NetworkX graphs
    :return: DataFrame with metrics for each network and replica
    """
    all_metrics = []

    for replica_index, networks in enumerate(replicas):
        replica_metrics = measure_network_metrics(networks)
        replica_metrics['replica'] = replica_index
        all_metrics.append(replica_metrics)

    return pd.concat(all_metrics, ignore_index=True)


def plot_network_stacked(df: pd.DataFrame, frag):
    """
    Plot the metrics as stacked bar charts with error bars.
    :param df: DataFrame containing the metrics to plot
    """
    # Ensure the DataFrame is sorted by 'step'
    df = df.sort_values(by='step')

    # Plot stacked bar chart
    fig, ax = plt.subplots(figsize=(40, 30))
    bottom = np.zeros(len(df))

    # Columns to plot (without '_ci' suffix)
    columns = ['waste', 'isolated', 'components', 'giant']

    # Use a color palette from matplotlib
    colors = plt.cm.Dark2.colors[:len(columns)]

    for col, color in zip(columns, colors):
        mean_values = df[col].values
        ax.bar(df['step'], mean_values, bottom=bottom, label=col,
               color=color, edgecolor='black',linewidth=0)
        bottom += mean_values

    ax.set_xlabel('Step', fontsize=24)
    ax.set_ylabel('Proportion of the network (%)', fontsize=24)
    plt.yticks(fontsize=24)
    plt.xticks(fontsize=24)
    ax.set_title(frag)
    plt.savefig(f'./figs/stack_proportion_{frag}.jpg')
    plt.show()


def calculate_statistics(df):
    """Calculate mean and 95% confidence interval for all columns in the dataframe."""
    result = []

    # Select all columns except 'step' and 'replica'
    columns_to_analyze = df.columns.difference(['step', 'replica'])

    for column in columns_to_analyze:
        mean_values = df.groupby('step')[column].mean()
        sem = df.groupby('step')[column].sem()  # Standard error of the mean
        confidence_interval = 1.96 * sem  # 95% confidence interval

        # Create a DataFrame for this column's statistics
        column_stats = pd.DataFrame({
            'step': mean_values.index,
            f'{column}': mean_values.values,
            f'{column}_ci': confidence_interval.values
        })

        result.append(column_stats)

    # Concatenate all column statistics DataFrames along the 'step' index
    result_df = pd.concat(result, axis=1)

    # Remove duplicate 'step' columns
    result_df = result_df.loc[:, ~result_df.columns.duplicated()]

    # Fill NaN values in confidence intervals with zeros
    for column in columns_to_analyze:
        result_df[f'{column}_ci'] = result_df[f'{column}_ci'].fillna(0)

    return result_df

pd.set_option('display.max_rows', None)

# make barchart to show the proprtion of strcutures in the network
data = data[frag]
x= data[2]
# print(max(x['step']))

def last_step(df):
    last_steps = x.groupby('replica')['step'].max()
    last_step = last_steps.median()
    return last_step
last_steps = x.groupby('replica')['step'].max()
print(last_steps.mean(), last_steps.median())

# Plot the distribution of the last steps
plt.hist(last_steps, bins=range(min(last_steps), max(last_steps) + 2), edgecolor='black')
plt.xlabel('Step')
plt.ylabel('Frequency')
plt.title('Distribution of Last Steps in Each Replica')
plt.xticks(range(min(last_steps), max(last_steps) + 1))  # Set x-ticks to be integer steps
plt.show()
print(x)


na_percentage = x.groupby('step')['het'].apply(lambda x: x.isna().mean() * 100)
print(na_percentage)
# networks = data[1]
# matrices = measure_network_metrics_replicas(networks)
# print(matrices)
# stats = calculate_statistics(matrices)
# print(stats)
#
#
# plot_network_stacked(matrices,frag=frag)
# #









############## plot snapshots of networks along steps of fragmentation
# plot_fragmentation(data)

####plot centrality vs fragmnetation
# plot_centrality(data,centrality='degree')

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

