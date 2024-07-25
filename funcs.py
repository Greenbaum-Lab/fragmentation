import pickle

import networkx as nx
import numpy as np
import pandas as pd
from infomap import Infomap
from matplotlib import pyplot as plt


########## general functions

def load_data(fragmentation_types, net, ignore):
    data = {}
    for frag_type in fragmentation_types:
        filename = f'RGG, {frag_type}_ignore_{ignore}.pickle'
        with open(filename, 'rb') as file:
            data[frag_type] = pickle.load(file)
    print("I finished loading!")
    return data


def calculate_statistics(df):
    """Calculate mean and 95% confidence interval."""
    column = df.columns.difference(['replica', 'step']).values.tolist()
    if len(column) > 1:
        column = column[0]
    mean_values = df.groupby('step')[column].mean()
    sem = df.groupby('step')[column].sem()  # Standard error of the mean
    confidence_interval = 1.96 * sem  # 95% confidence interval
    return mean_values, confidence_interval

def access_het_dist(frag_data:list):
    return frag_data[2]

def access_fst_dist(frag_data:list):
    return frag_data[4]

def access_het_mean(frag_data:list):
    return frag_data[3]

def access_fst_mean(frag_data:list):
    return frag_data[5]

def access_networks(frag_data:list):
    return frag_data[1]


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

def compute_modularity(net):

    im = Infomap(silent=True, markov_time=1, variable_markov_time=True,flow_model='undirected',num_trials=10)

    # Add edges to the Infomap instance
    for edge in net.edges():
        im.add_link(*edge)
    im.run()

    return im.codelength

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

                record['modularity'] = compute_modularity(net)

            # if 'connectivity' in measure:
            #     record['connectivity'] = weighted_algebraic_connectivity(net)

            if 'component' in measure:
                record['component'] = measure_giant_component(net)

            data.append(record)

    df = pd.DataFrame(data)

    return df

# # Calculate the means and standard deviations for the specified centrality measures
# mean_centrality = df.groupby('step').mean().drop(columns='replicate')
# std_centrality = df.groupby('step').std().drop(columns='replicate')
#
# return mean_centrality, std_centrality





########### extra functions
def plot_fragmentation(data, replica: int):
    """
    Plots network snapshot across fragmentation processes.
    each fragmentation in its own row.

    :param data: A dictionary of loaded network data, keyed by fragmentation type.
    """
    steps = [50, 100, 150, 200]
    fragmentation_types = list(data.keys())
    pos = nx.spring_layout(data['rand'][1][replica][0], k=0.2, iterations=20, seed=50)
    num_rows = len(fragmentation_types)
    num_cols = len(steps)
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 4 * num_rows))

    for row_idx, frag_type in enumerate(fragmentation_types):
        net_data = data[frag_type]

        for col_idx, step in enumerate(steps):
            net = net_data[1][replica][step]

            ax = axes[row_idx, col_idx] if num_rows > 1 else axes[col_idx]
            nx.draw_networkx(net, pos=pos, ax=ax, node_size=20, with_labels=False)
            ax.set_title(f"{step}" if row_idx == 0 else "", fontsize=22)  # Only set step number for the first column
            if col_idx == 0:
                # Label the rows with the fragmentation type
                ax.set_ylabel(frag_type, fontsize=36)

    plt.savefig("./figs/fragmentation processes.svg")
    plt.tight_layout()
    plt.show()



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


######## plot fragmentaion snapshots
# fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt']
# net = 'RGG'
# ignore = False
# data = load_data(fragmentation_types, net, ignore)
# # plot_fragmentation(data,replica=83)
