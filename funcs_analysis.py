import pickle
from statistics import mean

import networkx as nx
import pandas as pd
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
        plt.savefig(f'genetics general {measure}.jpg', format="jpg")
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

    plt.savefig("fragmentation processes.png")
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
    num_rows = math.ceil(len(fragmentation_types)/3)

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

        ax.set_xlabel('Step',fontsize=20)
        ax.set_ylabel('GC/Heterozygosity',fontsize=20)
        ax.set_title(frag_type,fontsize=20,ha='left',loc='left')

        ax.legend()
    plt.savefig('giant_component.jpg')
    plt.show()





from infomap import Infomap


def compute_modularity(net):
    im = Infomap(silent=True, markov_time=1, variable_markov_time=True)

    # Add edges to the Infomap instance
    for edge in net.edges():
        im.addLink(*edge)
    im.run()

    return im.codelength


def calculate_centrality(all_nets: list, measures: list = ['clustering', 'degree', 'modularity', 'connect']) -> (
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

            if 'trans' in measures:
                record['trans'] = nx.transitivity(net)

            if 'degree' in measures:
                degree = sum(nx.degree_centrality(net).values()) / len(net.nodes)
                record['degree'] = degree

            if 'connect' in measures:
                record['connect'] = nx.average_node_connectivity(net)

            if 'modularity' in measures:
                # partition = community_louvain.best_partition(net, resolution=1)
                # record['modularity'] = community_louvain.modularity(partition, net)
                record['modularity'] = compute_modularity(net)

            # if 'algebric' in measures:
            #     record['clustering'] = nx.algebraic_connectivity(net)
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
    plt.savefig(f'{centrality}.jpg')
    plt.show()