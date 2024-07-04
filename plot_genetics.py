from statistics import mean

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt

from funcs_analysis import load_data



def make_networks(n_nets: int, n_nodes: int, net_type) -> list:
    """
    create a list of networks
    :param n_nets: number of networks
    :param n_nodes: number of nodes
    :param connectivity: degree of connectivity
    :param net_type: type of network: ER, RGG, or SF
    :return: list of networks
    """
    nets = []
    for net in range(n_nets):

        if net_type == 'ER':
            net = nx.erdos_renyi_graph(n=n_nodes, p=0.2)
            nets.append(net)
        if net_type == 'RGG':
            net = nx.random_geometric_graph(n=n_nodes, radius=0.3)
            nets.append(net)
        if net_type == 'AB':
            net = nx.barabasi_albert_graph(n=n_nodes, m=5)
            nets.append(net)
        if net_type == 'SW':
            net = nx.watts_strogatz_graph(n=n_nodes,k=9, p=0.1)
            nets.append(net)

    return nets



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


def plot_data(data, index, ylabel, measure):
    """Plot data with mean and 95% confidence interval."""
    color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette
    plt.figure()
    # last_step = last_step(data['rand'][2])
    print(f"Last step for : {last_step}")

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
    # plt.xlim(None,265)
    plt.legend()
    plt.savefig(f'./figs/genetics_general_{measure}.jpg', format="jpg")
    plt.show()



def filter_intervals(df, interval_percentage=25):
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
    steps_to_include = steps_to_include[:4]

    # Filter the DataFrame to include only these steps
    filtered_df = df[df['step'].isin(steps_to_include)]
    return filtered_df




def plot_distribution(df, measure='het' or 'fst',type=str):
    """ Plot the distribution of a given measure across different steps.
    """
    # Create a figure and axes
    fig, ax = plt.subplots()
    # Get unique steps
    unique_steps = df['step'].unique()

    # Generate reversed color gradient
    colors = plt.cm.YlOrRd(np.linspace(0.3, 1, len(unique_steps)))[::-1]

    # Plot histogram for each step with increasing alpha
    for i, step in enumerate(unique_steps):
        if measure == 'fst':
            values = df[df['step'] == step]['fst']
        if measure == 'het':
            values = df[df['step'] == step]['het']
        ax.hist(values, bins=40, alpha=0.7, label=f'Step {step}', density=True,
                color=colors[i], edgecolor='black')

    # Set titles and labels
    ax.set_xlabel('Fst' if measure == 'fst' else 'Heterozygosity')
    ax.set_ylabel('Density (%) ')
    ax.legend()

    # Optional: set x and y limits
    # ax.set_xlim(0, 1.4)
    ax.set_ylim(0, 20)

    # Show the plot
    plt.savefig(f'./figs/dist_{measure}_{type}.jpg', format="jpg")
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
fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt', 'opt2', 'wrst']
# fragmentation_types = ['int']
net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)


# Plot fst and het along fragmentation
plot_data(data, 5, 'Pairwise Fst',measure='fst')
plot_data(data, 3, 'Heterozygosity',measure='heterozygosity')



##############plot distributions
##### one frag type each time
fragmentation_types = ['int']
net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)


df = filter_intervals(data[fragmentation_types][2])
plot_distribution(df,measure='het',type=fragmentation_types)



####################### plot individual nodes
################
# plot_nodes_all(data)
