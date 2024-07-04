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
from funcs import load_data, calculate_statistics

import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pickle

fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
fragmentation_types = ['rand', 'cor']

frag = 'div'

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

        if measure == 'component':
            sns.regplot(x='component', y='avg', data=merged, fit_reg=True, order=2,
                        truncate=True, scatter_kws={'s': 50, 'alpha': 0.01, 'color': color_palette(i)},
                        line_kws={'lw': 2, 'label': frag_type})

            max_val = max(plt.gca().get_xlim()[1], plt.gca().get_ylim()[1])
            min_val = min(plt.gca().get_xlim()[1], plt.gca().get_ylim()[1])
            plt.plot([0.05, max_val], [0, max_val], linestyle='--', color='black')
            plt.xlabel('Fraction of nodes in the largest component', fontsize=16)


    if measure == 'modularity':
        sns.regplot(x='modularity', y='avg', data=merged, fit_reg=True, order=2,
                    truncate=True, scatter_kws={'s': 50, 'alpha': 0.01, 'color': color_palette(i)},
                    line_kws={'lw': 2, 'label': frag_type})

        plt.gca().invert_xaxis()
        plt.ylim(-0.1, 1.1)
        plt.xlabel('Community structure', fontsize=16)


    plt.ylabel('Heterozygosity', fontsize=16)
    plt.legend()

    if save:
        plt.savefig(f'./figs/het_{measure}.jpg', format="jpg")
    plt.show()


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
        giant_component = measure_giant_component(network)
        isolated_nodes = measure_isolated_nodes(network)
        components = measure_components(network)
        waste = measure_waste(network)

        total = giant_component + isolated_nodes + components + waste

        # Round the first three metrics
        giant = round(giant_component / total, 2)
        isolated = round(isolated_nodes / total, 2)
        components = round(components / total, 2)

        # Adjust the last metric so the total sums up to 1
        waste = 1 - giant - isolated - components

        scaled_metrics = {
            "step": step,
            "giant": giant,
            "isolated": isolated,
            "components": components,
            "waste": waste,
        }

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



def calculate_statistics(df):
    """Calculate mean and 95% confidence interval for all columns in the dataframe."""
    result = []

    # Select all columns except 'step' and 'replica'
    columns_to_analyze = df.columns.difference(['step', 'replica'])

    for column in columns_to_analyze:
        mean_values = round(df.groupby('step')[column].mean(),2)
        # sem = df.groupby('step')[column].sem()  # Standard error of the mean
        # confidence_interval = 1.96 * sem  # 95% confidence interval

        # Create a DataFrame for this column's statistics
        column_stats = pd.DataFrame({
            'step': mean_values.index,
            f'{column}': mean_values.values,
            # f'{column}_ci': confidence_interval.values
        })

        result.append(column_stats)

    # Concatenate all column statistics DataFrames along the 'step' index
    result_df = pd.concat(result, axis=1)

    # Remove duplicate 'step' columns
    result_df = result_df.loc[:, ~result_df.columns.duplicated()]
    result_df['waste'] = 1 - result_df['giant'] - result_df['isolated'] - result_df['components']

    # Fill NaN values in confidence intervals with zeros
    # for column in columns_to_analyze:
    #     result_df[f'{column}_ci'] = result_df[f'{column}_ci'].fillna(0)

    return result_df


def plot_network_stacked_area(df: pd.DataFrame, frag: str):
    """
    Plot the metrics as stacked area charts.
    :param df: DataFrame containing the metrics to plot
    :param frag: Fragmentation type
    """
    # Ensure the DataFrame is sorted by 'step'
    df = df.sort_values(by='step')

    # Create a new figure and axes with a specific size
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define the columns to plot and the colors to use
    columns = ['waste', 'isolated', 'components', 'giant']
    colors = plt.cm.Dark2.colors[:len(columns)]

    # Prepare the data for the stackplot
    x_values = df['step'].values
    y_values = [df[col].values for col in columns]

    # Create the stackplot
    ax.stackplot(x_values, y_values, labels=columns, colors=colors, alpha=0.8)

    # Set the labels and title
    ax.set_xlabel('Step', fontsize=24)
    ax.set_ylabel('Proportion of the network (%)', fontsize=24)
    ax.set_title(frag)
    plt.ylim(0, 1)
    # Add a legend
    # plt.legend(loc='upper left')

    # Save the figure
    plt.savefig(f'./figs/stack_proportion_{frag}.jpg')

    # Display the plot
    plt.show()


#####plot stacks
# networks = data[1]
# matrices = measure_network_metrics_replicas(networks)
# stats = calculate_statistics(matrices)
# plot_network_stacked_area(stats,frag=frag)



####plot centrality vs heterozygosity
# plot_het_central(data, measure='modularity', save=True)


####plot centrality vs fragmnetation
# plot_centrality(data,centrality='connectivity')


##plot fst-distance relationship
# data = data[frag]
# matrix = data[7][0][20]
# net = data[1][0][20]
# distance_matrix = get_euclidean_matrix(net)
# perform_mantel_test(matrix,distance_matrix)

# distance_matrix = get_distance_matrix(net)
# plot_matrix_relationship(distance_matrix=distance_matrix,fst_matrix=matrix)

def calculate_mantel_correlation(data):
    # Initialize a list to store the results
    mantel_results = []

    # Get the number of replicas
    num_replicas = len(data[7])

    # Loop over all replicas
    for rep in range(num_replicas):
        # Get the number of steps for the current replica
        num_steps = len(data[7][rep])
        num_steps = 10

        # Loop over all steps
        for step in range(num_steps):
            # Get the FST matrix and network for the current step
            matrix = data[7][rep][step]
            net = data[1][rep][step]

            distance_matrix = get_euclidean_matrix(net)

            correlation = perform_mantel_test(matrix, distance_matrix, perms=999,
                                              method='pearson',print=False)[0]

            mantel_results.append({
                'replica': rep,
                'step': step,
                'correlation': correlation
            })

    return pd.DataFrame(mantel_results)


def plot_step_vs_correlation(df):
    plt.figure(figsize=(10, 6))
    plt.plot(df['step'], df['correlation'])
    plt.xlabel('Step')
    plt.ylabel('Correlation')
    plt.show()


from scipy import stats

def calculate_statistics(df):
    # Group the data by 'step'
    grouped = df.groupby('step')

    # Calculate the mean correlation for each step
    mean_correlation = grouped['correlation'].mean()

    # Calculate the standard error of the mean for each step
    sem_correlation = grouped['correlation'].sem()

    ci_correlation = 1.96 * sem_correlation  # 95% confidence interval

    # Create a new DataFrame for the results
    results = pd.DataFrame({
        'step': mean_correlation.index,
        'correlation_mean': mean_correlation.values,
        'correlation_ci': ci_correlation.values
    })

    return results

def plot_correlation_with_ci(df):
    plt.figure(figsize=(10, 6))
    sns.lineplot(x='step', y='correlation_mean', data=df)
    plt.fill_between(df['step'], df['correlation_mean'] - df['correlation_ci'], df['correlation_mean'] + df['correlation_ci'], color='b', alpha=0.1)
    plt.xlabel('Step')
    plt.ylabel('Correlation')
    plt.show()


def plot_correlation_with_ci(df, fragmentation_types):
    plt.figure(figsize=(10, 6))

    # If a single fragmentation type is provided, convert it to a list
    if isinstance(fragmentation_types, str):
        fragmentation_types = [fragmentation_types]

    # Loop over all fragmentation types
    for frag_type in fragmentation_types:
        # Filter the data for the current fragmentation type
        df_filtered = df[df['fragmentation_type'] == frag_type]

        # Plot the data for the current fragmentation type
        sns.lineplot(x='step', y='correlation_mean', data=df_filtered, label=frag_type)
        plt.fill_between(df_filtered['step'], df_filtered['correlation_mean'] - df_filtered['correlation_ci'], df_filtered['correlation_mean'] + df_filtered['correlation_ci'], alpha=0.1)

    plt.xlabel('Step')
    plt.ylabel('Correlation')
    plt.legend()
    plt.show()


def process_multiple_fragmentation_types(data, fragmentation_types):
    for frag_type in fragmentation_types:
        frag_data = data[frag_type]
        correlations = calculate_mantel_correlation(frag_data)
        stats = calculate_statistics(correlations)
        plot_correlation_with_ci(stats,fragmentation_types=frag_type)


# x=calculate_mantel_correlation(data)
# print(x)
# stats= calculate_statistics(x)
# print(stats)
# plot_correlation_with_ci(stats, fragmentation_types)



# fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
# fragmentation_types = ['rand', 'cor']
#
# net = 'RGG'
# ignore = False
# data = load_data(fragmentation_types, net, ignore)

