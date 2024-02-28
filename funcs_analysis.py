import math
import pickle
from statistics import mean

import networkx as nx
import numpy as np
import pandas as pd
from infomap import Infomap
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
        plt.close()
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
    steps_to_include = list(range(0, max_step, interval_step))

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
                         measures: list = ['clustering', 'degree', 'modularity','transitivity'
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

            if 'clustering' in measures:
                record['clustering'] = nx.average_clustering(net)

            if 'transitivity' in measures:
                record['transitivity'] = nx.transitivity(net)

            if 'degree' in measures:
                degree = sum(nx.degree_centrality(net).values()) / len(net.nodes)
                record['degree'] = degree

            if 'connect' in measures:
                record['connect'] = nx.average_node_connectivity(net)

            if 'modularity' in measures:
                # partition = community_louvain.best_partition(net, resolution=1)
                # record['modularity'] = community_louvain.modularity(partition, net)
                record['modularity'] = compute_modularity(net)

            if 'connectivity' in measures:
                record['connectivity'] = weighted_algebraic_connectivity(net)

            if 'component' in measures:
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
        mean_centrality, std_centrality = calculate_centrality(data[name][1], measures=[centrality])

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


def weighted_algebraic_connectivity(net):
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


def extract_nodes(df):
    """
    tracking each node separately and extract its heterozygosity
    :param df: df of heterozygous of all nodes for all replicates (het dens)
    :return: df for each node along the fragmentation for each replica
    """
    df.reset_index(drop=True, inplace=True)

    # count how many nodes are in a network
    nodes = np.argmax(df['step'] != 0)

    # Randomly select 5 unique nodes from the network
    random_indices = np.random.choice(nodes, 5, replace=False)

    # Initialize an empty DataFrame to store the selected rows
    selected_rows = pd.DataFrame()

    # Iterate over each randomly selected index and get the corresponding rows
    for index in random_indices:
        node_rows = df.iloc[index::nodes].copy()
        selected_rows = pd.concat([selected_rows, node_rows])

    return selected_rows





def plot_heterozygosity(data):
    """
    Plots heterozygosity along fragmentation steps for selected nodes in a single fragmentation type.

    :param data: DataFrame returned by extract_nodes function.
    """
    plt.figure(figsize=(10, 6))

    # Assuming 'data' already contains only the rows for the nodes of interest
    # Group by 'replica' to plot each replica's path separately
    for replica, group in data.groupby('replica'):
        plt.plot(group['step'], group['het'], label=f'Replica {replica}', alpha=0.5)

    plt.xlabel('Step')
    plt.ylabel('Heterozygosity')
    plt.title('Heterozygosity along fragmentation')
    plt.legend()
    plt.show()


def plot_heterozygosity_for_nodes(df, num_nodes=5):
    """
    Plots heterozygosity for a subset of nodes across steps.

    :param df: DataFrame containing heterozygosity data.
    :param num_nodes: Number of nodes to plot.
    """
    # Determine the total number of nodes per step in the first replica
    nodes_per_step = df[df['replica'] == df['replica'].unique()[0]]['step'].value_counts().iloc[0]

    # Select a subset of nodes to track
    node_indices = np.linspace(0, nodes_per_step - 1, num=num_nodes, dtype=int)

    plt.figure(figsize=(10, 6))

    # Plot heterozygosity for each selected node
    for node_idx in node_indices:
        # Extract rows corresponding to this node across all steps and replicas
        node_data = df.iloc[node_idx::nodes_per_step].copy()

        # Assuming steps are consistent across replicas, use the first replica's steps for x-axis
        steps = df[df['replica'] == df['replica'].unique()[0]]['step'].unique()

        plt.plot(steps, node_data['het'], label=f'Node {node_idx}', alpha=0.5)

    plt.xlabel('Step', fontsize=12)
    plt.ylabel('Heterozygosity', fontsize=12)
    plt.title('Heterozygosity for Selected Nodes Across Steps', fontsize=14)
    plt.legend()
    plt.tight_layout()
    plt.show()

fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)
pd.set_option('display.max_rows', None)

print(data['rand'][2])
def track_sampled_nodes_heterozygosity(df):
    """
    Samples 5 random nodes from step 0 of each replica and tracks their heterozygosity across all steps.

    :param df: DataFrame with columns ['step', 'het', 'replica'], where the index corresponds to nodes.
    :return: A dictionary with keys as replicas and values as DataFrames of tracked heterozygosity for sampled nodes.
    """
    sampled_nodes_data = {}  # Dictionary to store DataFrames for sampled nodes in each replica
    replicas = df['replica'].unique()

    for replica in replicas:
        # Filter for step 0 in the current replica
        step_0_data = df[(df['step'] == 0) & (df['replica'] == replica)]

        # Sample 5 random nodes
        sampled_nodes = step_0_data.sample(n=5, random_state=1).index

        # Track these nodes across all steps for the current replica
        node_data = pd.DataFrame()
        for node_idx in sampled_nodes:
            # Assuming consistent ordering of nodes across steps within a replica
            node_df = df[df.index % (df['step'].nunique() * 50) == node_idx]
            node_data = pd.concat([node_data, node_df], ignore_index=True)

        sampled_nodes_data[replica] = node_data

    return sampled_nodes_data


def plot_sampled_nodes_heterozygosity(sampled_nodes_data, replica):
    """
    Plots the heterozygosity of sampled nodes across steps for a specific replica.

    :param sampled_nodes_data: Dictionary with DataFrames of tracked heterozygosity for sampled nodes.
    :param replica: The replica to plot data for.
    """
    node_data = sampled_nodes_data[replica]
    plt.figure(figsize=(10, 6))

    # Plot heterozygosity for each sampled node
    for node_idx in node_data['index'].unique():
        node_df = node_data[node_data['index'] == node_idx]
        plt.plot(node_df['step'], node_df['het'], label=f'Node {node_idx}')

    plt.xlabel('Step')
    plt.ylabel('Heterozygosity')
    plt.title(f'Heterozygosity Across Steps for Sampled Nodes in Replica {replica}')
    plt.legend()
    plt.show()


# Assuming 'het_data' is your DataFrame containing the heterozygosity data for 'rand' fragmentation type
het_data = data['rand'][2]

# Track the heterozygosity for sampled nodes
sampled_nodes_data = track_sampled_nodes_heterozygosity(het_data)

# Plot the heterozygosity for sampled nodes in a specific replica, for example, replica 0
plot_sampled_nodes_heterozygosity(sampled_nodes_data, 0)




# Create a pivot DataFrame
pivot_df = het.pivot(columns='replica', index='step', values='het')

# Plot using the pivot DataFrame
plt.figure(figsize=(10, 6))
plt.plot(pivot_df.index, pivot_df, color='grey', alpha=0.1)

plt.xlabel('Step')
plt.ylabel('Heterozygosity')
plt.title('Heterozygosity along RANDOM fragmentation')
plt.show()