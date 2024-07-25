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

from funcs import load_data


def prepare_het_df(dat, step: int):
    """
    Process the input DataFrame to add a 'node' column for each replica,
    considering up to the first 50 rows for each replica.

    Parameters:
    - df: pandas.DataFrame with columns ['replica', 'step', 'het']

    Returns:
    - het_df: pandas.DataFrame, the processed DataFrame with an additional 'node' column
    """
    all_het = []
    nodes = 50
    df = dat[2]

    # Iterate over each unique replica
    for replica in df['replica'].unique():
        # Select rows based on the 'step' value and up to the first 50 rows for the current replica
        replica_net = df[(df['replica'] == replica) & (df['step'] == step)].head(nodes)
        # Generate a sequence for the 'node' column within each replica slice
        replica_net['node'] = range(replica_net.shape[0])

        all_het.append(replica_net)

    # Concatenate all processed slices into a single DataFrame and reset the index
    het_df = pd.concat(all_het).reset_index(drop=True)

    return het_df




def merge_het_central(dat, step: int, centrality: str, frag: str, log=bool):
    het = prepare_het_df(dat, step)
    central = calculate_node_centrality(dat, step, centrality)
    # remove zero values
    central = central[central['central'] != 0]

    if log == True:
        central['central'] = np.log10(central['central'])

    final_df = pd.merge(het, central)
    return final_df


def plot_node_centrality(dat, step: int, centrality: str, frag: str, log=bool):

    final_df = merge_het_central(dat, step, centrality, frag, log)
    sns.regplot(x='central', y='het', data=final_df, fit_reg=False, order=2,
                scatter_kws={'s': 50, 'alpha': 0.1, 'color': 'grey'})
    plt.ylabel("Heterozygosity", fontsize=18)
    plt.xlabel(centrality.capitalize(), fontsize=18)
    plt.ylim(0, 1.8)
    plt.savefig(f'./figs/node_{centrality}_{step}_{frag}.jpg', format="jpg")
    plt.show()




def het_central_process_level(data, frag: str,centrality: str):
    """
    Create a DataFrame for all steps and replicas

    Parameters:
        data (pd.DataFrame): Input data.
        frag (str): Fragmentation strategy or other parameter used in merge_het_central.

    Returns:
        pd.DataFrame: Concatenated DataFrame for all steps and replicas.
    """
    # Initialize a list to store the DataFrames for each step
    df_list = []
    steps = range(0,150)
    # Iterate over the range of steps
    for step in steps:
        # Generate the DataFrame for the current step
        step_df = merge_het_central(data, step, centrality, frag, False)
        # Append the DataFrame to the list
        df_list.append(step_df)

    # Concatenate all the DataFrames in the list into a single DataFrame
    results_df = pd.concat(df_list, ignore_index=True)

    return results_df


def calculate_node_centrality(dat, step: int, centrality: str):
    """
    Calculates the betweenness centrality for the first network in each replica
    and organizes the results into a DataFrame, ensuring that the specified step
    index is available to avoid IndexError.

    Parameters:
    - dat: A nested list where dat[1] contains replicas, and each replica contains networks.
    - step: The step index to look for in each replica.

    Returns:
    - central_df: A DataFrame with columns 'replica', 'node', and 'central'.
    """
    central_dict = {}

    # Iterate over each replica in dat[1]
    for rep_index in range(len(dat[1])):
        # Ensure the step index is within the bounds of the list for this replica
        if step < len(dat[1][rep_index]):
            net = dat[1][rep_index][step]  # Safely get the network at the given step

            if centrality == 'betweenness':
                central = nx.betweenness_centrality(net)

            if centrality == 'degree':
                central = nx.degree_centrality(net)

            central_dict[rep_index] = central

        else:
            print(f"Skipping replica {rep_index}: step {step} is out of range.")
            continue

    # Convert the dictionary to a DataFrame
    central_df = pd.DataFrame([
        {'replica': rep, 'node': node, 'central': centrality}
        for rep, centrality_dict in central_dict.items()
        for node, centrality in centrality_dict.items()
    ])

    return central_df

def merge_het_central(dat, step: int, centrality: str, frag: str, log=bool):
    """
    This function merges heterozygosity and centrality data into a single DataFrame.

    Parameters:
    dat (list): A nested list where dat[1] contains replicas, and each replica contains networks.
    step (int): The step index to look for in each replica.
    centrality (str): The type of centrality to calculate. Can be 'betweenness' or 'degree'.
    frag (str): The fragmentation strategy.
    log (bool): If True, the centrality values are transformed using the log10 function.

    Returns:
    final_df (DataFrame): A DataFrame with merged heterozygosity and centrality data.
    """
    # Prepare the heterozygosity DataFrame
    het = prepare_het_df(dat, step)

    # Calculate the node centrality
    central = calculate_node_centrality(dat, step, centrality)

    # Remove zero values from the centrality data
    central = central[central['central'] != 0]

    # If log is True, transform the centrality values using the log10 function
    if log == True:
        central['central'] = np.log10(central['central'])
        het['het'] = np.log10(het['het'])

    # Merge the heterozygosity and centrality data into a single DataFrame
    final_df = pd.merge(het, central)

    return final_df


def plot_node_centrality(dat, step: int, centrality: str, frag: str, log=bool):

    final_df = merge_het_central(dat, step, centrality, frag, log)
    sns.regplot(x='central', y='het', data=final_df, fit_reg=True, order=2,line_kws={'color': 'red'},
                scatter_kws={'s': 50, 'alpha': 0.1, 'color': 'blue'})
    plt.ylabel("Heterozygosity (log10)", fontsize=18)
    plt.xlabel(f"{centrality.capitalize()} (log10)", fontsize=18)
    # plt.ylim(0, 1.8)
    plt.savefig(f'./figs/node_{centrality}_{step}_{frag}.jpg', format="jpg")
    plt.show()


def het_central_process_level(data, frag: str,centrality: str):
    """
    Create a DataFrame for all steps and replicas

    Parameters:
        data (pd.DataFrame): Input data.
        frag (str): Fragmentation strategy or other parameter used in merge_het_central.

    Returns:
        pd.DataFrame: Concatenated DataFrame for all steps and replicas.
    """
    # Initialize a list to store the DataFrames for each step
    df_list = []
    steps = range(0,300)
    # Iterate over the range of steps
    for step in steps:
        # Generate the DataFrame for the current step
        step_df = merge_het_central(data, step, centrality, frag, False)
        # Append the DataFrame to the list
        df_list.append(step_df)

    # Concatenate all the DataFrames in the list into a single DataFrame
    results_df = pd.concat(df_list, ignore_index=True)

    return results_df




def compute_correlation(data, frag: str,centrality:str):
    """
    Calculate the Pearson correlation coefficient between 'het' and 'central'
    for each combination of 'step' and 'replica' in the DataFrame.

    Parameters:
        df (pd.DataFrame): Input DataFrame containing 'step', 'het', 'replica', and 'central' columns.

    Returns:
        pd.DataFrame: DataFrame with columns 'step', 'replica', and 'cor', containing the correlation values.
    """
    # Initialize a list to store the results
    results = []

    df = het_central_process_level(data, centrality, frag)


    # Group by 'step' and 'replica'
    grouped = df.groupby(['step', 'replica'])

    # Iterate over each group
    for (step, replica), group in grouped:

        # Check if the group has at least 2 rows
        if len(group) <= 2:
            continue

        # Calculate the Pearson correlation coefficient
        correlation = pearsonr(x=group['central'],y=group['het'])[0]
        # Append the results as a dictionary
        results.append({'step': step, 'replica': replica, 'cor': correlation})

    # Convert the results list to a DataFrame
    results_df = pd.DataFrame(results)

    return results_df


def calculate_node_centrality(dat, step: int, centrality: str):
    """
    Calculates the betweenness centrality for the first network in each replica
    and organizes the results into a DataFrame, ensuring that the specified step
    index is available to avoid IndexError.

    Parameters:
    - dat: A nested list where dat[1] contains replicas, and each replica contains networks.
    - step: The step index to look for in each replica.

    Returns:
    - central_df: A DataFrame with columns 'replica', 'node', and 'central' for betweenness centrality.
    """
    central_dict = {}

    # Iterate over each replica in dat[1]
    for rep_index in range(len(dat[1])):
        # Ensure the step index is within the bounds of the list for this replica
        if step < len(dat[1][rep_index]):
            net = dat[1][rep_index][step]  # Safely get the network at the given step

            if centrality == 'betweenness':
                central = nx.betweenness_centrality(net)

            if centrality == 'degree':
                central = dict(nx.degree(net))

            central_dict[rep_index] = central

        else:
            print(f"Skipping replica {rep_index}: step {step} is out of range.")
            continue

    # Convert the dictionary to a DataFrame
    central_df = pd.DataFrame([
        {'replica': rep, 'node': node, 'central': centrality}
        for rep, centrality_dict in central_dict.items()
        for node, centrality in centrality_dict.items()
    ])

    return central_df


def calculate_statistics_cor(correlation_df,frag_type):
    """
    Calculate the mean and 95% confidence interval for each step across all replicas for a single frag_type.

    Parameters:
        correlation_df (pd.DataFrame): DataFrame containing the correlation results with columns 'step', 'replica', and 'cor'.

    Returns:
        pd.DataFrame: DataFrame with columns 'step', 'mean_cor', 'ci_lower', 'ci_upper', containing the mean correlation and confidence interval bounds for each step.
    """
    # Initialize a list to store the statistics
    statistics = []

    # Group by 'step'
    grouped = correlation_df.groupby('step')

    # Iterate over each group
    for step, group in grouped:
        # Calculate mean correlation
        mean_cor = group['cor'].mean()

        # Calculate standard error
        se = group['cor'].std() / np.sqrt(len(group['cor']))

        # Calculate the 95% confidence interval
        ci_lower, ci_upper = norm.interval(0.95, loc=mean_cor, scale=se)

        # Append the results
        statistics.append({
            'step': step,
            'mean_cor': mean_cor,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'frag': frag_type
        })

    # Convert the statistics list to a DataFrame
    statistics_df = pd.DataFrame(statistics)

    return statistics_df





def compute_correlation_all(data, centrality:str):

    all_correlations = []
    for frag_type, dataset in data.items():
        cor_frag = compute_correlation(dataset, centrality, frag_type)
        cor_frag = calculate_statistics_cor(cor_frag,frag_type)
        all_correlations.append(cor_frag)

    return all_correlations



def plot_mean_with_ci(data_list):
    """
    Plot the mean correlation with shaded 95% confidence intervals for each 'frag' type.

    Parameters:
        data_list (list): List of DataFrames, each containing 'step', 'mean_cor', 'ci_lower', 'ci_upper', and 'frag' columns.
    """

    # Iterate over the list of DataFrames
    for df in data_list:
        # Extract the frag type (assuming it's the same for all rows in the df)
        frag = df['frag'].iloc[0]
        plt.plot(df['step'], df['mean_cor'], label=frag)
        plt.fill_between(df['step'], df['ci_lower'], df['ci_upper'], alpha=0.2)

    plt.xlabel('Step', fontsize=20)
    plt.ylabel('Correlation (r)', fontsize=20)
    plt.legend(title='Fragmentation Type')
    plt.savefig(f'./figs/cor_central.jpg', format="jpg")
    plt.show()






##### plot heterozygisuty vs. node centrality


# fragmentation_types = ['rand']
# net = 'RGG'
# ignore = False
# data = load_data(fragmentation_types, net, ignore)
#
#
# #### plot snapshot of hetrozygosity- centrality correlation
# plot_node_centrality(data['rand'],step=0,centrality='degree',log=True,frag=fragmentation_types)

### plot correlation between centrality and heterozygosity for all processes
# results = compute_correlation_all(data,centrality='degree')
# plot_mean_with_ci(results)
