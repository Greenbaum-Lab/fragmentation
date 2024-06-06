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

from funcs_analysis import prepare_het_df


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
    sns.regplot(x='central', y='het', data=final_df, fit_reg=True, order=2,
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


frag = 'dist'
with open(f'RGG, {frag}_ignore_False.pickle', 'rb') as file:
    rand = pickle.load(file)
print('finish')
central='degree'
step=10

# x = calculate_node_centrality(rand,step=step,centrality=central)
# print(x)
x=het_central_process_level(rand,frag=frag,centrality=central)
print(x)
#
df_step_50 = x[(x['step'] == 50) & (x['replica'] == 5)]

# Create the scatter plot


plt.figure(figsize=(10, 6))
sns.scatterplot(x='central', y='het', data=df_step_50, s=100, alpha=0.3)

# Add labels and title
plt.xlabel('Centrality', fontsize=14)
plt.ylabel('Heterozygosity', fontsize=14)
plt.title('Centrality vs Heterozygosity for Step 50', fontsize=16)

# Display the plot
plt.show()
print( pearsonr(x=df_step_50['central'], y=df_step_50['het']))

df_step_50 = x[(x['step'] == 100) & (x['replica'] == 5)]

# Create the scatter plot
plt.figure(figsize=(10, 6))
sns.scatterplot(x='central', y='het', data=df_step_50, s=100, alpha=0.3)

# Add labels and title
plt.xlabel('Centrality', fontsize=14)
plt.ylabel('Heterozygosity', fontsize=14)
plt.title('Centrality vs Heterozygosity for Step 50', fontsize=16)

# Display the plot
plt.show()
print( pearsonr(x=df_step_50['central'], y=df_step_50['het']))


df_step_50 = x[(x['step'] == 144) & (x['replica'] == 5)]

# Create the scatter plot
plt.figure(figsize=(10, 6))
sns.scatterplot(x='central', y='het', data=df_step_50, s=100, alpha=0.3)

# Add labels and title
plt.xlabel('Centrality', fontsize=14)
plt.ylabel('Heterozygosity', fontsize=14)
plt.title('Centrality vs Heterozygosity for Step 50', fontsize=16)
for line in range(0, df_step_50.shape[0]):
    plt.text(df_step_50.central.iloc[line], df_step_50.het.iloc[line],
             df_step_50.node.iloc[line],  size='large')

# Display the plot
plt.show()
print( pearsonr(x=df_step_50['central'], y=df_step_50['het']))


net=rand[1][5][51]
nx.draw_networkx(net)
plt.show()

net=rand[1][5][101]
nx.draw_networkx(net)
plt.show()

net=rand[1][5][145]
pos = nx.spring_layout(net, k=0.15, iterations=20)
nx.draw_networkx(net,pos)
plt.show()
