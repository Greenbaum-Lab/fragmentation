from statistics import mean

import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from funcs import calculate_statistics, access_het_dist, access_fst_dist
from funcs_analysis import load_data


#
#
#
# def make_networks(n_nets: int, n_nodes: int, net_type) -> list:
#     """
#     create a list of networks
#     :param n_nets: number of networks
#     :param n_nodes: number of nodes
#     :param connectivity: degree of connectivity
#     :param net_type: type of network: ER, RGG, or SF
#     :return: list of networks
#     """
#     nets = []
#     for net in range(n_nets):
#
#         if net_type == 'ER':
#             net = nx.erdos_renyi_graph(n=n_nodes, p=0.2)
#             nets.append(net)
#         if net_type == 'RGG':
#             net = nx.random_geometric_graph(n=n_nodes, radius=0.3)
#             nets.append(net)
#         if net_type == 'AB':
#             net = nx.barabasi_albert_graph(n=n_nodes, m=5)
#             nets.append(net)
#         if net_type == 'SW':
#             net = nx.watts_strogatz_graph(n=n_nodes,k=9, p=0.1)
#             nets.append(net)
#
#     return nets
#
# nets = make_networks(100, 50, 'RGG')
#
# num_edges = [net.number_of_edges() for net in nets]
#
# # Create a histogram of the number of edges
# plt.hist(num_edges, bins=50)
#
# plt.title('Distribution of Number of Edges')
# plt.xlabel('Number of Edges')
# plt.ylabel('Frequency')
#
# plt.show()


def find_breaking_point(networks):
    """
    find the index of the list where the network is no longer connected
    param: networks: list of networks across fragmentation for single replica
    """

    num_edges = networks[0].number_of_edges()

    for index, network in enumerate(networks):
        if not nx.is_connected(network):
            return index/num_edges * 100
    return None


def find_breakink_point_list(networks_list: list):
    """get the breaking point for all replicas"""
    breaking_point = []
    for net_list in networks_list:
        x = find_breaking_point(net_list)
        breaking_point.append(x)
    return breaking_point



def normalize_steps(data):
    """
    Normalize the step index of a Pandas Series or the 'step' column of a DataFrame
    to a percentage scale from 0 to 100.

    Parameters:
    data (pd.Series or pd.DataFrame): A Pandas Series with step indices or a DataFrame with a 'step' column.

    Returns:
    pd.Series or pd.DataFrame: A new Pandas Series with normalized step indices or a DataFrame with a normalized 'step' column.
    """
    if isinstance(data, pd.Series):
        # Normalize step values to percentage for Series
        normalized_index = (data.index / data.index.max()) * 100
        normalized_series = pd.Series(data.values, index=normalized_index)
        return normalized_series

    elif isinstance(data, pd.DataFrame) and 'step' in data.columns:
        # Normalize step values to percentage for DataFrame
        data['step'] = (data['step'] / data['step'].max()) * 100
        return data
    else:
        raise ValueError("Input must be a Pandas Series or a DataFrame with a 'step' column.")

def plot_data(data, index, ylabel, measure):
    """Plot data with mean and 95% confidence interval."""
    color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette
    plt.figure()

    # Plot each dataset's mean and confidence interval
    for frag_type, datasets in data.items():
        mean_values, confidence_interval = calculate_statistics(datasets, index)
        mean_values = normalize_steps(mean_values)
        confidence_interval = normalize_steps(confidence_interval)
        plt.plot(mean_values, label=frag_type.capitalize())
        plt.fill_between(mean_values.index, mean_values - confidence_interval, mean_values + confidence_interval,
                         alpha=0.2)
    # Add breaking points and other plot details
    for i, (frag_type, datasets) in enumerate(data.items()):
        breaking_point = mean(find_breakink_point_list(datasets[1]))
        plt.axvline(x=breaking_point, color=color_palette(i), ymax=0.05,linewidth=4)

    plt.xlabel('Fragmentation (%)', fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.legend()
    plt.savefig(f'./figs/genetics_general_{measure}.jpg', format="jpg")
    plt.show()


################### distribution ####################
def filter_intervals(data, interval_percentage=25,measure='het' or 'fst'):
    """
    Filter the DataFrame to include only specific intervals of steps.

    Args:
    data: dict of fragmetatation type.
    interval_percentage (int): The percentage interval for filtering steps.

    Returns:
    pd.DataFrame: Filtered DataFrame.
    """
    if measure == 'fst':
        df = access_fst_dist(data)
    else:
        df = access_het_dist(data)
    # Determine the maximum step value
    max_step = df['step'].max()

    # Calculate interval step based on the percentage
    interval_step = max_step * interval_percentage // 100

    # Create a list of steps to include
    steps_to_include = list(range(0, max_step, interval_step))
    steps_to_include = steps_to_include[:4]

    # Filter the DataFrame to include only these steps
    filtered_df = df[df['step'].isin(steps_to_include)]
    return filtered_df


def plot_distribution(df, measure='het' or 'fst', type=str):
    """ Plot the distribution of a given measure across different steps.
    """
    # Create a figure and axes
    fig, ax = plt.subplots()
    # Get unique steps
    unique_steps = df['step'].unique()

    # Generate reversed color gradient
    colors = plt.cm.YlGnBu(np.linspace(0, 1, len(unique_steps)))[::-1]

    # Plot histogram for each step with increasing alpha
    for i, step in enumerate(unique_steps):
        if measure == 'fst':
            values = df[df['step'] == step]['fst']
        if measure == 'het':
            values = df[df['step'] == step]['het']

        ax.hist(values, bins=40, alpha=0.4, label=f'Step {step}', density=True,
                color=colors[i], edgecolor='black')

    # Set titles and labels
    ax.set_xlabel('Fst' if measure == 'fst' else 'Heterozygosity', fontsize=20)
    ax.set_ylabel('Density (%) ',fontsize=20)
    ax.legend()

    # Optional: set x and y limits
    if measure == 'het':
        ax.set_xlim(0, 1.4)
    ax.set_ylim(0, 15)

    # Show the plot
    plt.savefig(f'./figs/dist_{measure}_{type}.jpg', format="jpg")
    plt.show()

def plot_distribution_het(df, type=str):
    """ Plot the distribution of heterozygosity across different steps as a ridgeline plot. """
    # Create a figure and axes
    fig, ax = plt.subplots()
    # Get unique steps
    unique_steps = df['step'].unique()

    # Generate reversed color gradient
    colors = plt.cm.YlGnBu(np.linspace(0, 1, len(unique_steps)))[::-1]

    # Plot histogram for each step with increasing alpha
    for i, step in enumerate(unique_steps):
        values = df[df['step'] == step]['het']
        ax.hist(values, bins=40, alpha=0.4, label=f'Step {step}', color=colors[i], edgecolor='black', density=True)

        if i < len(unique_steps) - 1:
        # Offset each step's histogram by a certain amount
            for rect in ax.patches:
                rect.set_y(rect.get_y() + 6)  # Adjust this value to change the vertical spacing between histograms
                ax.axhline(y=rect.get_y(), color='black', linewidth=0.5)  # Add a line below the histogram

            else:
                rect.set_y(0)

    ax.set_xlabel('Heterozygosity', fontsize=20)
    ax.set_ylabel('Density', fontsize=20)
    ax.set_yticks([])
    ax.legend()

    ax.set_xlim(0, 1.4)
    ax.set_ylim(0, (4 + len(unique_steps) * 6)) # Adjust this value to match the vertical spacing between histograms

    # Remove the box (rectangle) around the plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.savefig(f'./figs/dist_het_{type}.jpg', format="jpg")
    plt.show()





def plot_distribution_fst(df, type=str):
    """ Plot the distribution of heterozygosity across different steps as a ridgeline plot. """
    # Create a figure and axes
    fig, ax = plt.subplots()
    # Get unique steps
    unique_steps = df['step'].unique()

    # Generate reversed color gradient
    colors = plt.cm.YlGnBu(np.linspace(0, 1, len(unique_steps)))[::-1]

    # Plot histogram for each step with increasing alpha
    for i, step in enumerate(unique_steps):
        values = df[df['step'] == step]['fst']
        ax.hist(values, bins=40, alpha=0.4, label=f'Step {step}', color=colors[i], edgecolor='black', density=True)

        if i < len(unique_steps) - 1:
        # Offset each step's histogram by a certain amount
            for rect in ax.patches:
                rect.set_y(rect.get_y() + 6)  # Adjust this value to change the vertical spacing between histograms
                ax.axhline(y=rect.get_y(), color='black', linewidth=0.5)  # Add a line below the histogram
            else:
                ax.axhline(y=rect.get_y(), color='black', linewidth=0.5)  # Add a line below the histogram

    ax.set_xlabel('Fst', fontsize=20)
    ax.set_ylabel('Density', fontsize=20)
    ax.set_yticks([])
    ax.legend()

    ax.set_ylim(0, (4 + len(unique_steps) * 6)) # Adjust this value to match the vertical spacing between histograms

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.savefig(f'./figs/dist_fst_{type}.jpg', format="jpg")
    plt.show()





############### individual nodes ####################
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


#######################
####################### plot data
# fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt']
# # fragmentation_types = ['rand']
# net = 'RGG'
# ignore = False
# data = load_data(fragmentation_types, net, ignore)
#
# # Plot fst and het along fragmentation
# plot_data(data, 5, 'Pairwise Fst', measure='fst')
# plot_data(data, 3, 'Heterozygosity',measure='heterozygosity')
#

##############plot distributions
##### one frag type each time
# fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt']
fragmentation_types = ['div']
frag_type = fragmentation_types[0]

net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)

df = filter_intervals(data[frag_type],measure='het')
plot_distribution_het(df,type=frag_type)
df = filter_intervals(data[frag_type],measure='fst')
plot_distribution_fst(df,type=frag_type)


####################### plot individual nodes
################
# fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt']
# # fragmentation_types = ['rand']
# net = 'RGG'
# ignore = False
# data = load_data(fragmentation_types, net, ignore)
#
# plot_nodes_all(data)
