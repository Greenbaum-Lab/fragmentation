import pickle
from statistics import mean

import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from funcs import calculate_statistics, access_het_dist, access_fst_dist, normalize_steps, load_data, access_fst_mean, \
    access_het_mean, access_networks

# pd.set_option('display.max_rows', None)

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


# all data is a dictionary with keys as fragmentation types and values as lists of dataframes\lists
def process_data_for_single_type(frag_data, measure: str):
    """process data for plotting of single fragmentation type.
     calculate the mean and 95% confidence interval and breaking point of a network.
     Args:
        frag_data: dictionary of all data, keyed by fragmentation type.
        measure: The measure to plot ('fst' or 'het')."""

    if measure == 'fst':
        df = access_fst_mean(frag_data)
    if measure == 'het':
        df = access_het_mean(frag_data)

    df_stat = calculate_statistics(df)

    df_res = normalize_steps(df_stat)
    # breaking_point = mean(find_breakink_point_list(access_networks(frag_data)))

    return df_res #,breaking_point

def plot_all_fragmentation_types(data, measure: str):
    """Plot data for all fragmentation types with mean and 95% confidence interval and breaking point.
    Args:
        data: A dictionary of all data, keyed by fragmentation type.
        measure: The measure to plot ('fst' or 'het')."""

    color_palette = plt.get_cmap('tab10')
    plt.figure(figsize=(10, 6))

    # i-index; frag_type-string of fragmentation type; df-dataframe of avg fst\het
    for i, (frag_type, df) in enumerate(data.items()):

        color = color_palette(i)

        df_stat = process_data_for_single_type(data[frag_type], measure)

        plt.plot(df_stat['mean'], label=frag_type, color=color)
        plt.fill_between(df_stat.index, df_stat['mean'] - df_stat['sd'], df_stat['mean'] + df_stat['sd'],
                         alpha=0.2)
        # plt.axvline(x=breaking_point, color=color, ymax=0.05, linewidth=4)

    plt.xlabel('Fragmentation (%)', fontsize=34)
    plt.ylabel(measure, fontsize=34)
    # plt.legend()
    plt.ylim(-0.05, 1.05)
    plt.tick_params(axis='both', which='major', labelsize=25)
    plt.savefig(f'./figs/genetics_general_{measure}.svg', format="svg")
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
    plt.savefig(f'./figs/dist_{measure}_{type}.svg', format="svg")
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
        ax.hist(values, bins=40, alpha=0.6, label=f'Step {step}', color=colors[i], edgecolor='black', density=True)

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
    # Increase tick labels size
    ax.tick_params(axis='x', labelsize=16)  # Increase x-axis tick labels size
    ax.tick_params(axis='y', labelsize=16)  # Increase y-axis tick labels size

    # Remove the box (rectangle) around the plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.12)

    plt.savefig(f'./figs/dist_het_{type}.svg', format="svg")
    plt.show()

def plot_distribution_fst(df, type=str):
    """ Plot the distribution of heterozygosity across different steps as a ridgeline plot. """
    # Create a figure and axes
    fig, ax = plt.subplots()
    # Get unique steps
    unique_steps = df['step'].unique()

    # Generate reversed color gradient
    colors = plt.cm.YlOrRd(np.linspace(0, 1, len(unique_steps)))[::-1]

    # Plot histogram for each step with increasing alpha
    for i, step in enumerate(unique_steps):
        values = df[df['step'] == step]['fst']
        ax.hist(values, bins=40, alpha=0.8, label=f'Step {step}', color=colors[i], edgecolor='black', density=True)

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

    ax.set_ylim(0, (4 + len(unique_steps) * 6 + 4))  # vertical spacing between histograms
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.12)
    # Increase tick labels size
    ax.tick_params(axis='x', labelsize=16)  # Increase x-axis tick labels size
    ax.tick_params(axis='y', labelsize=16)  # Increase y-axis tick labels size

    # Remove the box (rectangle) around the plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.savefig(f'./figs/dist_fst_{type}.svg', format="svg")
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
    num_plots = len(fragmentation_types)
    num_cols = 2
    num_rows = (num_plots + num_cols - 1) // num_cols

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, num_rows * 5))
    axes = axes.flatten()

    for i, frag_type in enumerate(fragmentation_types):
        het = data[frag_type][2]
        het_nodes = extract_selected_nodes(het)

        # Create a unique identifier for each node across replicas
        het_nodes['node_replica_id'] = het_nodes['node_number'].astype(str) + '_replica_' + het_nodes['replica'].astype(str)
        het_nodes = normalize_steps(het_nodes)

        # Pivot the DataFrame
        pivot_df = het_nodes.pivot_table(index='step', columns='node_replica_id', values='het')
        pivot_df = normalize_steps(pivot_df)

        ax = axes[i]
        for column in pivot_df.columns:
            ax.plot(pivot_df.index, pivot_df[column], color='grey', alpha=0.2)

        ax.set_title(f'{frag_type} fragmentation', fontsize=20)
        ax.set_ylabel('Heterozygosity', fontsize=22)
        ax.tick_params(axis='both', which='major', labelsize=18)

    # Remove any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig('./figs/SUP_ind.svg')
    plt.show()


############### variance ####################
def calculate_variance(data, fragmentation_types: list):
    """
    This function calculates the variance of heterozygosity for each node in the network
    across all steps for each fragmentation type. The calculation is performed separately
    for each replica, and then the mean variance is calculated across replicas for each step.
    """
    all_data = []

    for frag_type in fragmentation_types:
        frag_data = access_het_dist(data[f'{frag_type}'])

        steps = frag_data['step'].unique()
        for step in steps:
            all_replicas = []
            for replica in frag_data['replica'].unique():
                replica_data = frag_data[(frag_data['replica'] == replica) & (frag_data['step'] == step)]
                variance = replica_data['het'].var()
                all_replicas.append(variance)
            mean_variance = np.mean(all_replicas)
            sd = np.std(all_replicas)
            all_data.append({'fragmentation_type': frag_type, 'step': step, 'variance': mean_variance, 'sd': sd})

    df = pd.DataFrame(all_data)
    df.to_csv('./variance.csv', index=False)
    return df


def plot_variance(df):
    """
    Plot the variance of heterozygosity for each fragmentation type across all steps.
    """
    color_palette = plt.get_cmap('tab10')
    plt.figure(figsize=(10, 6))

    fragmentation_types = df['fragmentation_type'].unique()
    for i, frag_type in enumerate(fragmentation_types):
        color = color_palette(i)
        frag_df = df[df['fragmentation_type'] == frag_type]
        frag_df['step'] = ((frag_df['step'] - frag_df['step'].min()) /
                           (frag_df['step'].max() - frag_df['step'].min()) * 100)
        plt.plot(frag_df['step'], frag_df['variance'], label=frag_type, color=color)
        plt.fill_between(frag_df['step'], frag_df['variance'] - frag_df['sd'], frag_df['variance'] + frag_df['sd'],
                         alpha=0.2, color=color)

    plt.xlabel('Fragmentation (%)', fontsize=20)
    plt.ylabel('Variance', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.savefig('./figs/paper figs/SUP_variance.svg', format='svg')
    plt.show()


#######################
####################### plot data
fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt','wrst']
# fragmentation_types = ['cor']
data = load_data(fragmentation_types)

plot_all_fragmentation_types(data, measure='fst')
plot_all_fragmentation_types(data, measure='het')


##############plot distributions
##### one frag type each time
# fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt','wrst']
# fragmentation_types = ['reg']
# frag_type = fragmentation_types[0]
#
# data = load_data(fragmentation_types)
#
# df = filter_intervals(data[frag_type],measure='het')
# plot_distribution_het(df,type=frag_type)
# df = filter_intervals(data[frag_type],measure='fst')
# plot_distribution_fst(df,type=frag_type)


####################### plot individual nodes
################
# fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt', 'wrst']
# # fragmentation_types = ['rand']
# data = load_data(fragmentation_types)
# plot_nodes_all(data)


############# calculate and plot variance across nodes in the network
# fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt','wrst']
# data = load_data(fragmentation_types)
# calculate_variance(data, fragmentation_types)
# df = pd.read_csv('./variance.csv')
# plot_variance(df)