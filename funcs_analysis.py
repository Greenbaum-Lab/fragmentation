import math
import pickle
from statistics import mean

import networkx as nx
import numpy as np
import pandas as pd
from infomap import Infomap
from joypy import joyplot
from matplotlib import pyplot as plt
import seaborn as sns
from mantel import test
from scipy.stats import pearsonr
from scipy.stats import norm


# pd.set_option('display.max_rows', None)

def load_data(fragmentation_types, net, ignore):
    data = {}
    for frag_type in fragmentation_types:
        filename = f'RGG, {frag_type}_ignore_{ignore}.pickle'
        with open(filename, 'rb') as file:
            data[frag_type] = pickle.load(file)
    print("I finished loading!")
    return data


def find_breaking_point(networks):
    """
    find the index of the list where the network is no longer connected
    """
    for index, network in enumerate(networks):
        if not nx.is_connected(network):
            return index
    return None


def find_breakink_point_list(networks: list):
    breaking_point = []
    for net in networks:
        x = find_breaking_point(net)
        breaking_point.append(x)
    return breaking_point

def calculate_statistics(df, index):
    """Calculate mean and 95% confidence interval."""
    mean_values = df[index].groupby('step')['avg'].mean()
    sem = df[index].groupby('step')['avg'].sem()  # Standard error of the mean
    confidence_interval = 1.96 * sem  # 95% confidence interval
    return mean_values, confidence_interval


def plot_data(data, index, ylabel, measure, save=False):
    """Plot data with mean and 95% confidence interval."""
    color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette
    plt.figure()

    # Plot each dataset's mean and confidence interval
    for frag_type, datasets in data.items():
        mean_values, confidence_interval = calculate_statistics(datasets, index)
        plt.plot(mean_values, label=frag_type.capitalize())
        plt.fill_between(mean_values.index, mean_values - confidence_interval, mean_values + confidence_interval,
                         alpha=0.2)

    # Add breaking points and other plot details
    for i, (frag_type, datasets) in enumerate(data.items()):
        breaking_point = mean(find_breakink_point_list(datasets[1]))
        plt.axvline(x=breaking_point, color=color_palette(i), ymax=0.1)

    plt.xlabel('Fragmentation step', fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.legend()
    if save:
        plt.savefig(f'./figs/genetics_general_{measure}.jpg', format="jpg")
    plt.show()


def filter_intervals(df, interval_percentage=15):
    """
    Filter the DataFrame to include only specific intervals of steps.

    Args:
    df (pd.DataFrame): The original DataFrame with 'step' and 'replica' columns.
    interval_percentage (int): The percentage interval for filtering steps.

    Returns:
    pd.DataFrame: Filtered DataFrame.
    """
    # Determine the maximum step value
    max_step = df['step'].max()

    # Calculate interval step based on the percentage
    interval_step = max_step * interval_percentage // 100

    # Create a list of steps to include
    steps_to_include = list(range(0, max_step - 50, interval_step))

    # Filter the DataFrame to include only these steps
    filtered_df = df[df['step'].isin(steps_to_include)]
    return filtered_df


def intervals(lst):
    """
    take snapshots of the process
    :param lst:
    :return:
    """
    if len(lst) <= 50:
        return lst
    n = 19  # number of bins (-1)
    interval = max((len(lst) - 1) // n, 1)
    return lst[:n * interval:interval] + [lst[-1]]



def plot_distribution(df, name, type):
    """
    Histogram ridgeline plots for 'fst' or 'het' data.
    :param df: DataFrame containing the data to plot.
    :param name: Name of the value to plot ('fst' or 'het').
    :param type: The type of analysis or data grouping.
    """

    plt.figure(figsize=(12, 12))
    # Using joyplot to create histograms
    fig, axes = joyplot(
        data=df[[name, 'step']],
        by='step',
        hist=True,
        bins=100,
        overlap=0.2,
        colormap=plt.cm.viridis,
        fade=False,
        range_style='all',
        linecolor="black",
        linewidth=0.1,
        normalize=True
    )

    # Set plot title
    title = f'{name.capitalize()} - {type} fragmentation'
    fig.suptitle(title, fontsize=18)
    plt.ylabel('Step', fontsize=18)
    plt.xlabel(name, fontsize=18)

    for ax in axes:
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlim(-0.1, 1.3)  # Set appropriate limits

    plt.savefig('./figs/''dist' + '_' + name + '_' + type + '.png')
    plt.show()


def plot_all_distributions(data, types=['fst', 'het']):
    """
    Plots distributions for all fragmentation types and data types ('fst', 'het').

    :param data: Dictionary containing loaded data for each fragmentation type.
    :param types: List of data types to plot ('fst', 'het').
    """
    for frag_type in data:
        for data_type in types:
            df_index = 4 if data_type == 'fst' else 2
            df = filter_intervals(data[frag_type][df_index])

            plot_distribution(df, data_type, frag_type)

            print(f"Plot generated for {frag_type} - {data_type}")


#


def plot_fragmentation(data):
    """
    Plots network snapshot across fragmentation processes.
    each fragmentation in its own row.

    :param data: A dictionary of loaded network data, keyed by fragmentation type.
    """
    steps = [0, 50, 100, 150, 200, 250]
    fragmentation_types = list(data.keys())
    pos = nx.spring_layout(data['rand'][1][10][0], k=0.2, iterations=20, seed=50)
    num_rows = len(fragmentation_types)
    num_cols = len(steps)
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 4 * num_rows))

    for row_idx, frag_type in enumerate(fragmentation_types):
        net_data = data[frag_type]

        for col_idx, step in enumerate(steps):
            net = net_data[1][10][step]

            ax = axes[row_idx, col_idx] if num_rows > 1 else axes[col_idx]
            nx.draw_networkx(net, pos=pos, ax=ax, node_size=20, with_labels=False)
            ax.set_title(f"{step}" if row_idx == 0 else "", fontsize=22)  # Only set step number for the first column
            if col_idx == 0:
                # Label the rows with the fragmentation type
                ax.set_ylabel(frag_type, fontsize=36)

    plt.savefig("./figs/fragmentation processes.png")
    plt.tight_layout()
    plt.show()


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


def compute_mean_std(data):
    """
    Helper function to compute mean and standard deviation for given data.
    """
    mean = data.groupby('step')['avg'].mean()
    confidence = data.groupby('step')['avg'].std()
    return mean, confidence


def plot_component_genetics(data):
    """
    Plot the fraction of nodes in a giant component with heterozygosity
    along fragmentation. Do it for all fragmentation types.
    :param data:
    :return:
    """
    fragmentation_types = list(data.keys())
    # to allow plotting any number of frag types
    num_rows = math.ceil(len(fragmentation_types) / 3)

    fig, axes = plt.subplots(num_rows, 3, figsize=(20, 4 * num_rows))
    axes = axes.flatten()  # Flatten the axes array for easy indexing

    for idx, frag_type in enumerate(fragmentation_types):
        data_frag = data[frag_type]

        giant_component = giant_component_replicates(data_frag[1])
        mean_gc_rand, conf_gc_rand = compute_mean_std(giant_component)
        mean_het_rand, conf_het_rand = compute_mean_std(data_frag[3])

        ax = axes[idx]

        ax.plot(mean_het_rand, label='Heterozygosity')
        ax.plot(mean_gc_rand, label='Giant component')

        ax.fill_between(mean_het_rand.index, mean_het_rand - conf_het_rand,
                        mean_het_rand + conf_het_rand, alpha=0.2)
        ax.fill_between(mean_gc_rand.index, mean_gc_rand - conf_gc_rand,
                        mean_gc_rand + conf_gc_rand, alpha=0.2)

        ax.set_xlabel('Step', fontsize=20)
        ax.set_ylabel('GC/Heterozygosity', fontsize=20)
        ax.set_title(frag_type, fontsize=20, ha='left', loc='left')

        ax.legend()
    plt.savefig('./figs/giant_component.jpg')
    plt.show()


def compute_modularity(net):

    im = Infomap(silent=True, markov_time=1, variable_markov_time=True,flow_model='undirected',num_trials=1)

    # Add edges to the Infomap instance
    for edge in net.edges():
        im.add_link(*edge)
    im.run()

    return im.codelength




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
    for i, nets in enumerate(all_nets):
        for step, net in enumerate(nets):
            record = {'replicate': i, 'step': step}

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

    # Calculate the means and standard deviations for the specified centrality measures
    mean_centrality = df.groupby('step').mean().drop(columns='replicate')
    std_centrality = df.groupby('step').std().drop(columns='replicate')

    return mean_centrality, std_centrality


def plot_centrality(data, centrality='modularity'):
    """
    Plots centrality measures and their confidence
    intervals against fragmentation steps.

    :param data: Dictionary containing network data for all fragmentation types.
    :param centrality: The centrality measure to plot.
    """
    names = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
    labels = ['Random', 'Correlated', 'Intrusive', 'Distance', 'Regressive', 'Divisive', 'Optimal']

    plt.figure(figsize=(10, 6))

    for name, label in zip(names, labels):
        # Calculate centrality and its standard deviation using your function
        mean_centrality, std_centrality = calculate_centrality(data[name][1], measure=[centrality])

        steps = mean_centrality.index

        # Plotting the centrality measure for the current fragmentation type
        plt.plot(steps, mean_centrality[centrality], label=label)

        lower_bound = mean_centrality[centrality] - std_centrality[centrality]
        upper_bound = mean_centrality[centrality] + std_centrality[centrality]

        # Plotting the confidence interval as a shaded area
        plt.fill_between(steps, lower_bound, upper_bound, alpha=0.2)

    plt.xlabel('Step', fontsize=22)
    plt.ylabel(centrality.capitalize(), fontsize=18)
    plt.title(f'{centrality.capitalize()} along Fragmentation', fontsize=22)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./figs/{centrality}.jpg')
    plt.show()


def weighted_algebraic_connectivity(G):
    """Calculate the weighted algebraic connectivity of a graph with disconnected components."""
    components = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    total_weight = sum(len(comp) for comp in components)  # Total weight based on the number of nodes in each component
    weighted_connectivity = sum(len(comp) * nx.algebraic_connectivity(comp) for comp in components) / total_weight
    return weighted_connectivity


def weighted_algebraic_connectivity1(net):
    """Calculate the weighted algebraic connectivity of a graph with disconnected components."""
    if nx.is_connected(net):
        return nx.algebraic_connectivity(net)

    components = list(nx.connected_components(net))
    total_weight = nx.number_of_nodes(net)

    connectivity_sum = 0
    for comp in components:
        comp_size = len(comp)
        comp_subgraph = net.subgraph(comp)
        # Avoid calculation for single-node components as algebraic connectivity would be 0
        if comp_size > 1:
            comp_connectivity = nx.algebraic_connectivity(comp_subgraph)
            connectivity_sum += comp_size * comp_connectivity
        # For a single-node component, you could decide to add or not add to the sum, depending on interpretation
        # In this context, skipping as algebraic connectivity is not defined for single nodes in a meaningful way

    weighted_connectivity = connectivity_sum / total_weight
    return weighted_connectivity


def compute_degree_distributions(data, frag_type, step):
    """
    Computes pooled degree distribution for a specified step across all iterations for a given fragmentation type.

    :param data: Dictionary containing the loaded network data, keyed by fragmentation type.
    :param frag_type: The fragmentation type to analyze.
    :param step: Fragmentation step to analyze.
    :return: A tuple containing the degrees and their counts.
    """
    all_degrees = []  # List to collect all degrees across all iterations at the specified step

    for iteration_networks in data[frag_type][1]:
        # Check if the current iteration has the specified step
        if step < len(iteration_networks):
            # Access the network at the specified step within this iteration
            network_at_step = iteration_networks[step]
            # Extend the collected degrees with degrees from this network
            all_degrees.extend([deg for _, deg in network_at_step.degree()])

    # Compute the degree distribution from the pooled degrees
    degree_counts = np.bincount(all_degrees)
    deg = np.arange(len(degree_counts))

    return deg, degree_counts


def plot_degree_distributions(data):
    """
    Plots degree distributions across fragmentation processes for specified steps,
    with each fragmentation type in its own row.
    """
    steps = [0, 50, 100, 150, 200, 250]
    fragmentation_types = list(data.keys())
    num_rows = len(fragmentation_types)
    num_cols = len(steps)
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 4 * num_rows), constrained_layout=True)

    for row_idx, frag_type in enumerate(fragmentation_types):
        for col_idx, step in enumerate(steps):
            deg, degree_counts = compute_degree_distributions(data, frag_type, step)
            ax = axes[row_idx, col_idx]
            ax.bar(deg, degree_counts, color='grey', alpha=0.7)
            ax.set_title(f"Step {step}" if row_idx == 0 else "", fontsize=18)
            ax.set_xlabel('Degree', fontsize=14)
            ax.set_ylabel('Count', fontsize=14)
            ax.set_xlim(0, max(deg) + 1)
            plt.tick_params(axis='both', which='major', labelsize=22)  # Increase tick labels font size

            if col_idx == 0:
                ax.set_ylabel('')
                # Label the rows with the fragmentation type
                ax.text(-0.1, 0.5, frag_type, fontsize=18, ha='right', va='center', transform=ax.transAxes, rotation=90)
    plt.savefig("./figs/degree_distributions.png")
    plt.show()


def select_nodes(df, num_nodes=1):
    """
    Selects a specified number of random node indices for each replica.

    :param df: DataFrame containing the heterozygosity data (dat[2]).
    :param num_nodes: Number of nodes to select per replica.
    :return: A dictionary with replicas as keys and lists of selected node indices as values.
    """
    selection_dict = {}
    for replica in df['replica'].unique():
        df_replica = df[df['replica'] == replica]
        nodes_per_replica = np.argmax(df_replica['step'].to_numpy()[1:] != df_replica['step'].to_numpy()[:-1]) + 1
        random_indices = np.random.choice(nodes_per_replica, min(num_nodes, nodes_per_replica), replace=False)
        selection_dict[replica] = random_indices
    return selection_dict


def extract_selected_nodes(df):
    """
    Extracts rows for the selected nodes across all steps for each replica.

    :param df: DataFrame containing heterozygosity data.
    :param selection_dict: A dictionary with replicas as keys and lists of selected node indices as values.
    :return: DataFrame with the extracted rows, including a node_number column.
    """

    selection_dict = select_nodes(df, num_nodes=1)
    selected_rows_across_replicas = pd.DataFrame()
    for replica, indices in selection_dict.items():
        df_replica = df[df['replica'] == replica]
        nodes_per_replica = np.argmax(df_replica['step'].to_numpy()[1:] != df_replica['step'].to_numpy()[:-1]) + 1
        for index in indices:
            node_rows = df_replica.iloc[index::nodes_per_replica].copy()
            node_rows['node_number'] = index + 1  # Assign node number
            selected_rows_across_replicas = pd.concat([selected_rows_across_replicas, node_rows], ignore_index=True)
    return selected_rows_across_replicas


# def plot_nodes(df, frag_type):
#     """
#     Plots the heterozygosity ('het') values for each node across steps using a pivot table approach.
#
#     :param df: DataFrame with 'het' values, 'step', 'node_number', and 'replica'.
#     """
#     # Create a unique identifier for each node across replicas
#     df['node_replica_id'] = df['node_number'].astype(str) + '_replica_' + df['replica'].astype(str)
#
#     # Pivot the DataFrame
#     pivot_df = df.pivot_table(index='step', columns='node_replica_id', values='het')
#
#     plt.figure(figsize=(10, 6))
#     # Plotting each column in the pivot table
#     for column in pivot_df.columns:
#         plt.plot(pivot_df.index, pivot_df[column], color='grey', alpha=0.2)
#
#     plt.xlabel('Step', fontsize=18)
#     plt.ylabel('Heterozygosity', fontsize=18)
#     plt.title(f'{frag_type} fragmentation', fontsize=20)
#     plt.tick_params(axis='both', which='major', labelsize=18)  # Increase tick labels font size
#     plt.tight_layout()
#     plt.show()

def plot_nodes(df, frag_type):
    """
    Plots the heterozygosity ('het') values for each node across steps using a pivot table approach.

    :param df: DataFrame with 'het' values, 'step', 'node_number', and 'replica'.
    """
    # Create a unique identifier for each node across replicas
    df['node_replica_id'] = df['node_number'].astype(str) + '_replica_' + df['replica'].astype(str)

    # Pivot the DataFrame
    pivot_df = df.pivot_table(index='step', columns='node_replica_id', values='het')

    plt.figure(figsize=(10, 6))
    # Plotting each column in the pivot table
    for column in pivot_df.columns:
        plt.plot(pivot_df.index, pivot_df[column], color='grey', alpha=0.2)

    plt.xlabel('Step', fontsize=18)
    plt.ylabel('Heterozygosity', fontsize=18)
    plt.title(f'{frag_type} fragmentation', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tight_layout()
    plt.show()


def plot_nodes_all(data):
    fragmentation_types = list(data.keys())

    for frag_type in fragmentation_types:
        het = data[frag_type][2]
        het_nodes = extract_selected_nodes(het)

        plot_nodes(het_nodes, frag_type)


def plot_het_central(data: dict, measure: str, save=bool):
    fragmentation_types = list(data.keys())

    plt.figure()

    for frag_type in fragmentation_types:
        het = compute_mean_std(data[frag_type][3])[0]
        central = calculate_centrality(data[frag_type][1], measure=measure)[0]

        plt.plot(het, central, label=frag_type.capitalize())

    plt.xlabel('Heterozygosity', fontsize=16)
    plt.ylabel('Transitivity', fontsize=16)
    plt.legend()
    plt.gca().invert_xaxis()

    if save == True:
        plt.savefig(f'./figs/het_{measure}.jpg', format="jpg")
    plt.show()


def prepare_het_df(dat, step: int):
    """
    Process the input DataFrame to add a 'node' column for each replica,
    considering up to the first 50 rows for each replica.

    Parameters:
    - df: pandas.DataFrame with columns ['replica', 'step', 'het']

    Returns:
    - het_df: pandas.DataFrame, the processed DataFrame with an additional 'node' column
    """
    all_het = []
    nodes = 50
    df = dat[2]

    # Iterate over each unique replica
    for replica in df['replica'].unique():
        # Select rows based on the 'step' value and up to the first 50 rows for the current replica
        replica_net = df[(df['replica'] == replica) & (df['step'] == step)].head(nodes)
        # Generate a sequence for the 'node' column within each replica slice
        replica_net['node'] = range(replica_net.shape[0])

        all_het.append(replica_net)

    # Concatenate all processed slices into a single DataFrame and reset the index
    het_df = pd.concat(all_het).reset_index(drop=True)

    return het_df


def calculate_node_centrality(dat, step: int, centrality: str):
    """
    Calculates the specified centrality for the first network in each replica
    and organizes the results into a DataFrame, ensuring that the specified step
    index is available to avoid IndexError.

    Parameters:
    - dat: A nested list where dat[1] contains replicas, and each replica contains networks.
    - step: The step index to look for in each replica.
    - centrality: The type of centrality to calculate ('betweenness' or 'degree').

    Returns:
    - central_df: A DataFrame with columns 'replica', 'node', and 'central' for the specified centrality.
    """
    central_dict = {}

    # Iterate over each replica in dat[1]
    for rep_index in range(len(dat[1])):
        # Ensure the step index is within the bounds of the list for this replica
        if step > len(dat[1][rep_index]):
            print(f"Skipping replica {rep_index}: step {step} is out of range.")
            continue

        # Safely get the network at the given step
        net = dat[1][rep_index][step]

        if centrality == 'betweenness':
            central = nx.betweenness_centrality(net)
        elif centrality == 'degree':
            central = nx.degree_centrality(net)
        else:
            raise ValueError("Unsupported centrality type. Use 'betweenness' or 'degree'.")

        central_dict[rep_index] = central

    # Convert the dictionary to a DataFrame
    central_df = pd.DataFrame([
        {'replica': rep, 'node': node, 'central': centrality}
        for rep, centrality_dict in central_dict.items()
        for node, centrality in centrality_dict.items()
    ])

    return central_df




def calculate_node_centrality(dat, step: int, centrality: str):
    """
    Calculates the betweenness centrality for the first network in each replica
    and organizes the results into a DataFrame, ensuring that the specified step
    index is available to avoid IndexError.

    Parameters:
    - dat: A nested list where dat[1] contains replicas, and each replica contains networks.
    - step: The step index to look for in each replica.

    Returns:
    - central_df: A DataFrame with columns 'replica', 'node', and 'central' for betweenness centrality.
    """
    central_dict = {}

    # Iterate over each replica in dat[1]
    for rep_index in range(len(dat[1])):
        # Ensure the step index is within the bounds of the list for this replica
        if step < len(dat[1][rep_index]):
            net = dat[1][rep_index][step]  # Safely get the network at the given step

            if centrality == 'betweenness':
                central = nx.betweenness_centrality(net)

            if centrality == 'degree':
                central = nx.degree_centrality(net)
            central_dict[rep_index] = central

        else:
            print(f"Skipping replica {rep_index}: step {step} is out of range.")
            continue

    # Convert the dictionary to a DataFrame
    central_df = pd.DataFrame([
        {'replica': rep, 'node': node, 'central': centrality}
        for rep, centrality_dict in central_dict.items()
        for node, centrality in centrality_dict.items()
    ])

    return central_df

def merge_het_central(dat, step: int, centrality: str, frag: str, log=bool):
    het = prepare_het_df(dat, step)
    central = calculate_node_centrality(dat, step, centrality)
    # remove zero values
    central = central[central['central'] != 0]

    if log == True:
        central['central'] = np.log10(central['central'])

    final_df = pd.merge(het, central)
    return final_df

#
# def merge_het_central(dat, step: int, centrality: str, frag: str, log=bool):
#     """
#     Merge heterogeneity and centrality data, removing rows with zero centrality.
#     If only one row remains after removing zero centrality rows, discard the entire DataFrame.
#
#     Parameters:
#     - dat: Input data.
#     - step: The step index to process.
#     - centrality: The type of centrality to calculate ('betweenness' or 'degree').
#     - frag: Fragmentation strategy or other parameter.
#     - log: Whether to apply log transformation to the centrality values.
#
#     Returns:
#     - final_df: The merged DataFrame after processing.
#     """
#     het = prepare_het_df(dat, step)
#     central = calculate_node_centrality(dat, step, centrality)
#
#     # Remove zero values
#     central = central[central['central'] != 0]
#
#     if log:
#         central['central'] = np.log10(central['central'])
#
#     final_df = pd.merge(het, central, on='node')  # Assuming 'node' is the common column for merging
#     return final_df
#


def plot_node_centrality(dat, step: int, centrality: str, frag: str, log=bool):

    final_df = merge_het_central(dat, step, centrality, frag, log)
    sns.regplot(x='central', y='het', data=final_df, fit_reg=True, order=2,
                scatter_kws={'s': 50, 'alpha': 0.1, 'color': 'grey'})
    plt.ylabel("Heterozygosity", fontsize=18)
    plt.xlabel(centrality.capitalize(), fontsize=18)
    plt.ylim(0, 1.8)
    plt.savefig(f'./figs/node_{centrality}_{step}_{frag}.jpg', format="jpg")
    plt.show()




def het_central_process_level(data, frag: str,centrality: str):
    """
    Create a DataFrame for all steps and replicas

    Parameters:
        data (pd.DataFrame): Input data.
        frag (str): Fragmentation strategy or other parameter used in merge_het_central.

    Returns:
        pd.DataFrame: Concatenated DataFrame for all steps and replicas.
    """
    # Initialize a list to store the DataFrames for each step
    df_list = []
    steps = 300
    # Iterate over the range of steps
    for step in range(0, steps):
        # Generate the DataFrame for the current step
        step_df = merge_het_central(data, step, centrality, frag, False)
        # Append the DataFrame to the list
        df_list.append(step_df)

    # Concatenate all the DataFrames in the list into a single DataFrame
    results_df = pd.concat(df_list, ignore_index=True)

    return results_df




def compute_correlation(data, frag: str,centrality:str):
    """
    Calculate the Pearson correlation coefficient between 'het' and 'central'
    for each combination of 'step' and 'replica' in the DataFrame.

    Parameters:
        df (pd.DataFrame): Input DataFrame containing 'step', 'het', 'replica', and 'central' columns.

    Returns:
        pd.DataFrame: DataFrame with columns 'step', 'replica', and 'cor', containing the correlation values.
    """
    # Initialize a list to store the results
    results = []

    df = het_central_process_level(data, centrality, frag)


    # Group by 'step' and 'replica'
    grouped = df.groupby(['step', 'replica'])

    # Iterate over each group
    for (step, replica), group in grouped:

        # Check if the group has at least 2 rows
        if len(group) <= 2:
            continue

        # Calculate the Pearson correlation coefficient
        correlation = pearsonr(x=group['central'],y=group['het'])[0]
        # Append the results as a dictionary
        results.append({'step': step, 'replica': replica, 'cor': correlation})

    # Convert the results list to a DataFrame
    results_df = pd.DataFrame(results)

    return results_df




def calculate_statistics_cor(correlation_df,frag_type):
    """
    Calculate the mean and 95% confidence interval for each step across all replicas for a single frag_type.

    Parameters:
        correlation_df (pd.DataFrame): DataFrame containing the correlation results with columns 'step', 'replica', and 'cor'.

    Returns:
        pd.DataFrame: DataFrame with columns 'step', 'mean_cor', 'ci_lower', 'ci_upper', containing the mean correlation and confidence interval bounds for each step.
    """
    # Initialize a list to store the statistics
    statistics = []

    # Group by 'step'
    grouped = correlation_df.groupby('step')

    # Iterate over each group
    for step, group in grouped:
        # Calculate mean correlation
        mean_cor = group['cor'].mean()

        # Calculate standard error
        se = group['cor'].std() / np.sqrt(len(group['cor']))

        # Calculate the 95% confidence interval
        ci_lower, ci_upper = norm.interval(0.95, loc=mean_cor, scale=se)

        # Append the results
        statistics.append({
            'step': step,
            'mean_cor': mean_cor,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'frag': frag_type
        })

    # Convert the statistics list to a DataFrame
    statistics_df = pd.DataFrame(statistics)

    return statistics_df


def compute_correlation_all(data, centrality:str):

    all_correlations = []
    for frag_type, dataset in data.items():
        cor_frag = compute_correlation(dataset, centrality, frag_type)
        cor_frag = calculate_statistics_cor(cor_frag,frag_type)
        all_correlations.append(cor_frag)

    return all_correlations



def plot_mean_with_ci(data_list):
    """
    Plot the mean correlation with shaded 95% confidence intervals for each 'frag' type.

    Parameters:
        data_list (list): List of DataFrames, each containing 'step', 'mean_cor', 'ci_lower', 'ci_upper', and 'frag' columns.
    """

    # Iterate over the list of DataFrames
    for df in data_list:
        # Extract the frag type (assuming it's the same for all rows in the df)
        frag = df['frag'].iloc[0]
        plt.plot(df['step'], df['mean_cor'], label=frag)
        plt.fill_between(df['step'], df['ci_lower'], df['ci_upper'], alpha=0.2)

    plt.xlabel('Step', fontsize=20)
    plt.ylabel('Correlation (r)', fontsize=20)
    plt.legend(title='Fragmentation Type')
    plt.savefig(f'./figs/cor_central.jpg', format="jpg")
    plt.show()


fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
# fragmentation_types = ['rand', 'cor']

net = 'RGG'
ignore = False
# data = load_data(fragmentation_types, net, ignore)



frag = 'rand'
with open(f'RGG, {frag}_ignore_False.pickle', 'rb') as file:
    rand = pickle.load(file)





def make_het_dist(het_list: list, ignore: bool=False) -> pd.DataFrame:
    """
    Takes a list of heterozygosity values and returns a DataFrame.
    If ignore_ones is set to True, it will ignore all values of 1.

    Args:
    het_list : List of heterozygosity vectors
    ignore_ones : If True, ignores all values of 1

    Returns:
    DataFrame with a column of all the heterozygosity values and the corresponding fragmentation step.
    """
    df = pd.DataFrame(het_list)

    if ignore:
        df = df.replace(0.02, np.nan)  # replace all 1s with NaN

    df = df.stack().rename_axis(('step', 'delete')).reset_index(name='het')
    df = df.drop(columns=['delete'])

    return df











# get df of het of a single node
# rand = rand[2]
# df = extract_selected_nodes(rand)
# df = df[df['replica'] == 1]
# df = df['het']
# df = pd.DataFrame(df).reset_index(drop=True)
# print(df)
# df.to_csv('df.csv')


def get_distance_matrix(net, default_distance=50):
    nodes = list(net.nodes())
    n = len(nodes)
    distance_matrix = np.full((n, n), default_distance)  # Initialize matrix with default distance
    node_index = {node: idx for idx, node in enumerate(nodes)}  # Map nodes to indices

    # Calculate shortest paths using Floyd-Warshall algorithm
    # This considers all path lengths and sets distances for all connected pairs
    path_lengths = dict(nx.all_pairs_dijkstra_path_length(net))

    for i, distances in path_lengths.items():
        for j, dist in distances.items():
            distance_matrix[node_index[i]][node_index[j]] = dist

    return distance_matrix



def plot_matrix_relationship(distance_matrix, fst_matrix, method='pearson', perms=999):
    # Convert all zeros to NaN in both matrices
    #convert 50 to NaN. 50 is the defeault value for isolated nodes
    distance_matrix = np.where((distance_matrix == 0) | (distance_matrix == 50), np.nan, distance_matrix)
    fst_matrix = np.where(fst_matrix == 0, np.nan, fst_matrix)

    # Perform Mantel test, expecting a dictionary as a return value
    result = test(fst_matrix, distance_matrix, perms=perms, method=method, ignore_nans=True)
    print(f"Correlation: {result[0]}")
    print(f"P-value: {result[1]}")

    # Flatten the matrices for plotting, ignoring NaN values
    flat_matrix1 = distance_matrix.flatten()
    flat_matrix2 = fst_matrix.flatten()

    valid_indices = ~np.isnan(flat_matrix1) & ~np.isnan(flat_matrix2)
    flat_matrix1 = flat_matrix1[valid_indices]
    flat_matrix2 = flat_matrix2[valid_indices]

    # Create scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(flat_matrix1, flat_matrix2, color='blue', edgecolor='k', alpha=0.7)

    # Add labels and title
    plt.xlabel('Euclidean distance')
    plt.ylabel('Fst')
    plt.title('Random fragmentation in the 150th step')

    # Add a line of best fit
    m, b = np.polyfit(flat_matrix1, flat_matrix2, 1)
    plt.plot(flat_matrix1, m * flat_matrix1 + b, color='red')
    #
    # coeffs = np.polyfit(flat_matrix1, flat_matrix2, 2)  # Quadratic fit
    # p = np.poly1d(coeffs)  # Create polynomial function
    # t = np.linspace(min(flat_matrix1), max(flat_matrix1), 500)
    # plt.plot(t, p(t), color='red')

    plt.show()


def get_euclidean_matrix(net):
    # Extract node positions into a numpy array
    nodes = list(net.nodes())
    positions = np.array([net.nodes[node]['pos'] for node in nodes])

    # Calculate the difference matrix for each dimension
    diff = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]

    # Compute the Euclidean distance matrix
    distance_matrix = np.linalg.norm(diff, axis=-1)

    return distance_matrix


# [1-all networks][replica no.][step number]
# [2-all heterozygosity][replica no.][step number]




