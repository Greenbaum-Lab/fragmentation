import pickle

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from funcs import access_het_dist


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


def assign_node_numbers(df, nodes_per_step=50):
    """
    Assigns node numbers for each node in each replica.

    :param df: DataFrame containing the heterozygosity data (dat[2]).
    :param nodes_per_step: Number of nodes in each step.
    :return: DataFrame with the original data and an additional 'node_number' column.
    """
    dfs = []  # List to collect DataFrames
    df = access_het_dist(df)
    for replica in df['replica'].unique():
        df_replica = df[df['replica'] == replica]
        for index in range(nodes_per_step):
            node_rows = df_replica.iloc[index::nodes_per_step].copy()
            node_rows['node_number'] = index + 1  # Assign node number
            dfs.append(node_rows)
    df_with_node_numbers = pd.concat(dfs, ignore_index=True)  # Concatenate all DataFrames at once
    return df_with_node_numbers

pd.set_option('display.max_rows', None)

def get_max_het_node(df):
    """
    Gets the node in each replica that has the highest 'het' in the last step of each replica.

    :param df: DataFrame with 'replica', 'step', 'het', and 'node_number' columns.
    :return: Dictionary with replicas as keys and node numbers with the highest 'het' in the last step as values.
    """
    # Get the last step of each replica
    last_steps = df.groupby('replica')['step'].max()

    # Select only the rows that belong to the last step of each replica
    df_last_steps = df[df.apply(lambda row: row['step'] == last_steps[row['replica']], axis=1)]

    # Get the index of the row with the highest 'het' in each group
    max_het_indices = df_last_steps.groupby('replica')['het'].idxmax()

    # Select the rows with the highest 'het' in the last step of each replica
    df_max_het = df.loc[max_het_indices]

    # Create a dictionary with replicas as keys and node numbers as values
    max_het_nodes = df_max_het.set_index('replica')['node_number'].to_dict()

    return max_het_nodes


def extract_steps_for_nodes(df, max_het_nodes):
    """
    Extracts all steps for each node in its corresponding replica from the DataFrame.

    :param df: DataFrame with 'replica', 'step', 'het', and 'node_number' columns.
    :param max_het_nodes: Dictionary with replicas as keys and node numbers with the highest 'het' in the last step as values.
    :return: DataFrame with the extracted rows.
    """
    extracted_rows = pd.DataFrame()
    for replica, node_number in max_het_nodes.items():
        df_replica_node = df[(df['replica'] == replica) & (df['node_number'] == node_number)]
        extracted_rows = pd.concat([extracted_rows, df_replica_node], ignore_index=True)
    return extracted_rows


def export_het_csv(data, frag: str):
    """
    Export the heterozygosity data to a CSV file.
    :param data: full fragmentation data
    :param frag: fragmentation type
    :return: csv file
    """
    data = assign_node_numbers(rand)
    surviving_nodes = get_max_het_node(data)
    final_df = extract_steps_for_nodes(data, surviving_nodes)
    final_df.to_csv(f'{frag}_het.csv')
    print('File saved successfully')
    return final_df


frag = 'rand'
with open(f'RGG, {frag}_ignore_False_d06.pickle', 'rb') as file:
    rand = pickle.load(file)
print('finish')

#plot het for a specific node
df = assign_node_numbers(rand)
df = df[(df['replica'] == 69) & (df['node_number'] == 5)]
plt.plot(df['step'], df['het'])
plt.show()

export_het_csv(rand, frag)


