import math
import pickle
from statistics import mean

import networkx as nx
import numpy as np
import pandas as pd
from infomap import Infomap
from joypy import joyplot
from matplotlib import pyplot as plt


def load_data(fragmentation_types, net, ignore):
    data = {}
    for frag_type in fragmentation_types:
        filename = f'RGG, {frag_type}_ignore_{ignore}.pickle'
        with open(filename, 'rb') as file:
            data[frag_type] = pickle.load(file)
    print("I finished loading!")
    return data


def plot_data(data, index, ylabel, measure, save=bool):
    color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette
    plt.figure()
    for frag_type, datasets in data.items():
        mean_values = datasets[index].groupby('step')['avg'].mean()
        confidence = datasets[index].groupby('step')['avg'].std()
        plt.plot(mean_values, label=frag_type.capitalize())
        plt.fill_between(mean_values.index, mean_values - confidence, mean_values + confidence, alpha=0.2)

    # Add breaking points and other plot details
    for i, (frag_type, datasets) in enumerate(data.items()):
        breaking_point = mean(find_breakink_point_list(datasets[1]))
        plt.axvline(x=breaking_point, color=color_palette(i), ymax=0.1)

    plt.xlabel('Fragmentation step', fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.legend()
    if save == True:
        plt.savefig(f'./figs/genetics general {measure}.jpg', format="jpg")
    plt.show()


def filter_intervals(df, interval_percentage=10):
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
    steps_to_include = list(range(0, max_step - 40, interval_step))

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


def plot_distribution(df, name, type):
    """
    Desnsity plots for 'fst' or 'het' data.
    :param df: DataFrame containing the data to plot.
    :param name: Name of the value to plot ('fst' or 'het').
    """

    # Plotting
    plt.figure(figsize=(8, 8))
    fig, axes = joyplot(
        data=df[[name, 'step']],
        by='step', overlap=3,
        colormap=plt.cm.viridis, fade=True, range_style='all',
        linecolor="black", linewidth=0.1
    )

    # Set plot title
    title = f'{name.capitalize()} - {type} Fragmentation'
    fig.suptitle(title, fontsize=18)
    plt.ylabel('Step', fontsize=18)
    plt.xlabel(type, fontsize=18)

    for ax in axes:
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlim(-0.1, 1.3)

    plt.savefig('./figs/' + name +'_'+ type + '.png')
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
    im = Infomap(silent=True, markov_time=1, variable_markov_time=True)

    # Add edges to the Infomap instance
    for edge in net.edges():
        im.addLink(*edge)
    im.run()

    return im.codelength


def calculate_centrality(all_nets: list,
                         measure: list = ['clustering', 'degree','component',
                                           'modularity','transitivity',
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
            ax.set_xlim(0, max(deg)+1)
            plt.tick_params(axis='both', which='major', labelsize=22)  # Increase tick labels font size

            if col_idx == 0:
                ax.set_ylabel('')
                # Label the rows with the fragmentation type
                ax.text(-0.1, 0.5, frag_type, fontsize=18, ha='right', va='center', transform=ax.transAxes, rotation=90)
    plt.savefig("./figs/degree_distributions.png")
    plt.show()


def select_nodes(df, num_nodes=5):
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

    selection_dict = select_nodes(df, num_nodes=5)
    selected_rows_across_replicas = pd.DataFrame()
    for replica, indices in selection_dict.items():
        df_replica = df[df['replica'] == replica]
        nodes_per_replica = np.argmax(df_replica['step'].to_numpy()[1:] != df_replica['step'].to_numpy()[:-1]) + 1
        for index in indices:
            node_rows = df_replica.iloc[index::nodes_per_replica].copy()
            node_rows['node_number'] = index + 1  # Assign node number
            selected_rows_across_replicas = pd.concat([selected_rows_across_replicas, node_rows], ignore_index=True)
    return selected_rows_across_replicas


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
    plt.tick_params(axis='both', which='major', labelsize=18)  # Increase tick labels font size
    plt.tight_layout()
    plt.show()

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


fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
net = 'RGG'
ignore = False
# data = load_data(fragmentation_types, net, ignore)

# with open('RGG, rand_ignore_True.pickle', 'rb') as file:
#     rand = pickle.load(file)
#
# with open('RGG, div_ignore_True.pickle', 'rb') as file:
#     div = pickle.load(file)

print("finsh load!!!!")


def plot_het_central(data:dict,measure:str,save=bool):

    fragmentation_types = list(data.keys())

    plt.figure()

    for frag_type in fragmentation_types:
        het = compute_mean_std(data[frag_type][3])[0]
        central = calculate_centrality(data[frag_type][1], measure=measure)[0]

        plt.plot(het,central, label=frag_type.capitalize())

    plt.xlabel('Heterozygosity', fontsize=16)
    plt.ylabel('Average node connectivity', fontsize=16)
    plt.legend()
    if save == True:
        plt.savefig(f'./figs/het_{measure}.jpg', format="jpg")
    plt.show()


# plot_het_central(data,measure='modularity',save=True)

net = nx.random_geometric_graph(100,0.1)
nx.draw_networkx(net)
plt.show()
net2 = nx.random_geometric_graph(100,0.15)
nx.draw_networkx(net)
plt.show()
print(weighted_algebraic_connectivity1(net))
print(weighted_algebraic_connectivity1(net))

print(weighted_algebraic_connectivity(net2))
