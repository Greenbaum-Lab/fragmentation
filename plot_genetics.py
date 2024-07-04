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



# Generate 100 RGG networks with exactly 250 edges each
networks = make_rgg(n_nets=100, n_nodes=50, target_edges=250)


num_edges = [net.number_of_edges() for net in networks]
print(num_edges)
# Plot the histogram of the number of edges
plt.figure(figsize=(10, 6))
plt.hist(num_edges, bins=100, edgecolor='black')
plt.title('Distribution of Number of Edges in 100 RGG Networks')
plt.xlabel('Number of Edges')
plt.ylabel('Frequency')
plt.show()


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