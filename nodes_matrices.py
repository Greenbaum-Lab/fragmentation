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

from funcs import load_data, access_networks, access_het_dist


def add_nodes(data, nodes=50):
    """
    Add a 'node' column to the DataFrame for each replica, considering up to the first 50 rows for each replica.

    Parameters:
    - data (pd.DataFrame): DataFrame with columns ['replica', 'step', 'het'].
    - nodes (int): Number of nodes to consider for each replica. Default is 50.

    Returns:
    - pd.DataFrame: Processed DataFrame with an additional 'node' column.
    """

    all_het = []
    df = access_het_dist(data)
    steps = df['step'].max()

    for replica in df['replica'].unique():
        for step in range(steps):
            # Select rows based on the 'step' value and up to the first 50 rows for the current replica
            replica_net = df[(df['replica'] == replica) & (df['step'] == step)].head(nodes)
            # Generate a sequence for the 'node' column within each replica slice
            replica_net['node'] = range(replica_net.shape[0])
            all_het.append(replica_net)
    # Concatenate all processed slices into a single DataFrame and reset the index
    het_df = pd.concat(all_het).reset_index(drop=True)
    return het_df

def calculate_node_centrality(data, centrality: str):
    """
    Calculate the specified centrality for each network in each replica and organize the results into a DataFrame.

    Parameters:
    - data (list): A nested list where data[1] contains replicas, and each replica contains networks.
    - centrality (str): The type of centrality to calculate. Can be 'bet' for betweenness or 'degree' for degree centrality.

    Returns:
    - pd.DataFrame: DataFrame with columns ['node', 'central', 'step', 'replica'] containing the centrality values.
    """
    centrality_funcs = {
        'bet': nx.betweenness_centrality,
        'degree': lambda net: dict(nx.degree(net))
    }

    nets = access_networks(data)
    central_list = []

    for replica in range(len(nets)):
        for step in range(len(nets[0])):
            net = nets[replica][step]
            central = centrality_funcs[centrality](net)
            central_df = pd.DataFrame.from_dict(central, orient='index', columns=['central'])
            central_df['step'] = step
            central_df['replica'] = replica
            central_df = central_df.reset_index().rename(columns={'index': 'node'})
            central_list.append(central_df)

    return pd.concat(central_list)


def make_het_central(data, centrality: str,):
    """
    merge heterozygosity and centrality data into a single DataFrame.
    remove centrality-zero values
    """
    het = add_nodes(data)
    central = calculate_node_centrality(data, centrality)
    final_df = pd.merge(het, central, on=['node', 'step', 'replica'])
    final_df = final_df[final_df['central'] > 0]

    return final_df

def make_het_central_all(data, centrality: str):
    """
    get ehterozygosity and centrality data for all fragmentation types
    """
    for frag, data in data.items():
        final_df = make_het_central(data, centrality)
        # write csv
        final_df.to_csv(f'./csv/{centrality}_het_{frag}.csv', index=False)
    return final_df



def plot_node_centrality(df, step: int, centrality: str, frag: str, log=bool):


    sns.regplot(x='central', y='het', data=df, fit_reg=True, order=1,logistic=True,
                scatter_kws={'s': 50, 'alpha': 0.7, 'color': 'blue'})
    plt.ylabel("Heterozygosity", fontsize=18)

    plt.tick_params(axis='both', labelsize=16)
    plt.savefig(f'./figs/node_{centrality}_{step}_{frag}.svg', format="svg")
    plt.show()


def compute_correlation(df):
    """
    Calculate the Pearson correlation coefficient between 'het' and 'central'
    for each combination of 'step' and 'replica' in the DataFrame.

    Parameters:
        df (pd.DataFrame): Input DataFrame containing 'step', 'het', 'replica', and 'central' columns.

    Returns:
        pd.DataFrame: DataFrame with columns 'step', 'replica', and 'cor', containing the correlation values.
    """
    results = []
    grouped = df.groupby(['step', 'replica'])

    for (step, replica), group in grouped:
        if len(group) < 2: # need at least 2 data points to calculate correlation
            continue
        cor, pval = pearsonr(group['het'], group['central'])
        results.append({'step': step, 'replica': replica, 'cor': cor, 'pval': pval})

    results_df = pd.DataFrame(results)
    return results_df


##### plot heterozygisuty vs. node centrality

fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
# data = load_data(fragmentation_types)

### compute correlation between heterozygosity and centrality for all processes
### change centrality name
# make_het_central_all(data, 'bet')
# for frag in fragmentation_types:
#     df = pd.read_csv(f'./csv/bet_het_{frag}.csv')
#     results_df = compute_correlation(df)
#     results_df.to_csv(f'./csv/bet_het_{frag}_cor.csv', index=False)


##### plot single degree-het correlation
# df = pd.read_csv(f'./csv/bet_het_rand.csv')
# df = (df[(df['step'] == 0) & (df['replica'] == 0)])
# df = np.log10(df)
# ax = sns.regplot(x='central', y='het', data=df,
#             fit_reg=True, order=1)
# plt.ylabel('Heterozygosity (log)', fontsize=20)
# plt.xlabel('betweeness (log)', fontsize=20)
# plt.tick_params(axis='both', labelsize=16)
# # plt.legend(title='Fragmentation Type')
# plt.savefig(f'./figs/bet_het_singlecor.svg', format="svg")
# plt.show()

#### plot correlation het-degree for all replicates during fragmetiaon
# plt.figure(figsize=(10, 6))
# for frag in fragmentation_types:
#     df = pd.read_csv(f'./csv/degree_het_{frag}_cor.csv')
#     df = np.log10(df)
#     sns.lineplot(x='step', y='cor',data=df)
# plt.xlabel('Step', fontsize=20)
# plt.ylabel('Correlation (r)-betweenness', fontsize=20)
# # plt.ylim(0,1)
# plt.legend()
# plt.tick_params(axis='both', labelsize=16)
# # plt.savefig(f'./figs/bet_het_cor.svg', format="svg")
# plt.show()