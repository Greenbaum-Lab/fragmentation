import random
import networkx as nx
import pickle

import numpy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
random.seed(45)

def add_replica(df, step_column='step'):
    """
    add a column with the number of replica for each row based on 'step'
    """
    # Find where the step resets
    df.loc[:, 'replica'] = (df[step_column] < df[step_column].shift(1)).astype(int)

    # Compute the cumulative sum to get the replica number
    df.loc[:, 'replica'] = df['replica'].cumsum()

    return df
def extract_nodes(df):
    """
    tracking each node separately and extract its heterozygosity
    :param df: df of heterozygous of all nodes for all replicates (het dens)
    :return: df for each node along the fragmentation for each replica
    """
    # count how many nodes are in a network
    nodes = np.argmax(df['step'] != 0)

    # Randomly select 5 unique nodes from the network
    random_indices = np.random.choice(nodes, 5, replace=False,seed=4)

    # Initialize an empty DataFrame to store the selected rows
    selected_rows = pd.DataFrame()

    # Iterate over each randomly selected index and get the corresponding rows
    for index in random_indices:
        node_rows = df.iloc[index::nodes].copy()
        selected_rows = pd.concat([selected_rows, node_rows])

    return selected_rows

def plot_network_realization(data: tuple, step: int, save: bool = False):
    """
    Visualize a network based on betweenness centrality and heterozygosity for a specified step.

    Parameters:
    - data (tuple): A tuple containing relevant data.
                    data[1] is a list of replicates, with each replica being a list of networks.
                    data[2] is a DataFrame with columns 'step' and 'het' (heterozygosity).
    - step (int): The step at which to visualize the network.
    - save (bool, optional): Whether to save the figure. Default is False.

    Returns:
    None
    """

    # Extract relevant network for the specified step
    nets = data[1]
    replica = random.choice(range(len(nets)))  # Choose a random replicate from the list
    net = nets[replica][step]

    # Extract corresponding heterozygosity data
    het_df = add_replica(data[2])
    het_df = het_df.loc[(het_df['step'] == step) & (het_df['replica'] == replica)]
    het = {index: value for index, value in enumerate(het_df['het'])}

    # Calculate and normalize betweenness centrality for the selected network
    bet = nx.betweenness_centrality(net)
    max_betweenness = max(bet.values())
    normalized_betweenness = {node: value / max_betweenness for node, value in bet.items()}
    colors = plt.cm.Reds([normalized_betweenness[node] for node in net.nodes()])

    # Layout and visualization
    pos = nx.spring_layout(net, k=0.2, iterations=20)
    nx.draw_networkx_nodes(net, node_color=colors, pos=pos, edgecolors='black',
                           linewidths=1.5, node_size=[v * 300 for v in het.values()])
    nx.draw_networkx_edges(net, pos=pos, edge_color='gray')
    plt.title("Color: Betweenness Centrality | Size: Heterozygosity")

    # Save the figure if specified
    if save:
        plt.savefig("net_betweenness.svg", format='svg')

    # Display the plot
    plt.show()



def plot_network_realization(data: tuple, step: int, centrality,save: bool = False):
    """
    Visualize a network based on betweenness centrality and heterozygosity for a specified step.

    Parameters:
    - data (tuple): A tuple containing relevant data.
                    data[1] is a list of replicates, with each replica being a list of networks.
                    data[2] is a DataFrame with columns 'step' and 'het' (heterozygosity).
    - step (int): The step at which to visualize the network.
    - save (bool, optional): Whether to save the figure. Default is False.

    Returns:
    None
    """

    # Extract relevant network for the specified step
    nets = data[1]
    replica = random.choice(range(len(nets)))  # Choose a random replicate from the list
    net = nets[replica][step]

    # Extract corresponding heterozygosity data
    het_df = add_replica(data[2])
    het_df = het_df.loc[(het_df['step'] == step) & (het_df['replica'] == replica)]
    het = {index: value for index, value in enumerate(het_df['het'])}

    # Calculate and normalize betweenness centrality for the selected network
    if centrality == "bet":
        bet = nx.betweenness_centrality(net)
    if centrality == "clos":
        bet = nx.closeness_centrality(net)
    max_betweenness = max(bet.values())
    normalized_betweenness = {node: value / max_betweenness for node, value in bet.items()}
    colors = plt.cm.Reds([normalized_betweenness[node] for node in net.nodes()])

    # Layout and visualization
    pos = nx.spring_layout(net, k=0.2, iterations=20)
    nx.draw_networkx_nodes(net, node_color=colors, pos=pos, edgecolors='black',
                           linewidths=1.5, node_size=[v * 300 for v in het.values()])
    nx.draw_networkx_edges(net, pos=pos, edge_color='gray')
    plt.title("Color: Betweenness Centrality | Size: Heterozygosity")

    # Save the figure if specified
    if save:
        plt.savefig("net_betweenness.svg", format='svg')

    # Display the plot
    plt.show()




def plot_regression(df: pd.DataFrame, save=True) -> None:
    """
    Plot a regression between heterozygosity and betweenness attributes.

    Parameters:
    - df (pd.DataFrame): DataFrame containing 'het' and 'bet' columns.

    Returns:
    None
    """

    plt.figure(figsize=(8, 6))
    sns.regplot(x='bet', y='het', data=df, order=2, scatter_kws={'s':50,'alpha':0.1,'color':'grey'})

    # # Compute correlation coefficient and p-value
    # r, p = stats.pearsonr(df['bet'], df['het'])
    # r_squared = r ** 2
    #
    # # Annotate the plot with r^2 and p-value
    # plt.annotate(f'$R^2$ = {r_squared:.3f}', xy=(0.1, 0.9), xycoords='axes fraction', fontsize=14)
    # plt.annotate(f'p-value = {p:.5f}', xy=(0.1, 0.85), xycoords='axes fraction', fontsize=14)

    plt.xlabel('Betweenness', fontsize=14)
    plt.ylabel('Heterozygosity', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    if save:
        plt.savefig("cor.svg", format='svg')

    plt.show()



def centrality_correlation(data: tuple,centrality, step: int) -> None:
    """
    Compute and visualize the correlation between betweenness centrality
    and heterozygosity for nodes in a network at a specified step.

    Parameters:
    - data (tuple): A tuple containing relevant data.
                    data[1] is a list of replicates with each replica being a list of networks.
                    data[2] is a DataFrame with columns 'step' and 'het' (heterozygosity).
    - step (int): The step at which to compute the correlation.

    Returns:
    None
    """

    # 1. Extract relevant network and compute betweenness centrality
    nets = data[1]

    # if the net is broken before the step-skip it
    if centrality == "bet":
        bet = [value for net in nets if step < len(net) for value in nx.betweenness_centrality(net[step]).values()]

    if centrality == "clos":
        bet = [value for net in nets if step < len(net) for value in nx.closeness_centrality(net[step]).values()]

    # 2. Extract heterozygosity values for the given step
    het_df = data[2]
    het_values = het_df[het_df['step'] == step]['het'].tolist()

    # 3. Combine heterozygosity and betweenness centrality values for correlation analysis
    df = pd.DataFrame({'het': het_values, 'bet': bet})

    # 4. Visualize correlation
    plot_regression(df)




with open('RGG, dist_ignore_False.pickle', 'rb') as file:
    rand = pickle.load(file)
#
# het = rand[2]
# het.reset_index(drop=True, inplace=True)
# het_values = extract_nodes(het)
# het_values = add_replica(het_values)
#
# # Create a pivot DataFrame
# pivot_df = het_values.pivot(columns='replica', index='step', values='het')
#
# # Plot using the pivot DataFrame
# plt.figure(figsize=(10, 6))
# plt.plot(pivot_df.index, pivot_df, color='grey', alpha=0.1)
#
# plt.xlabel('Step')
# plt.ylabel('Heterozygosity')
# plt.title('Heterozygosity along DIVISIVE fragmentation')
# plt.show()
#
#
plot_network_realization(rand, step=150, centrality="clos")
centrality_correlation(rand, step=150, centrality="bet")

