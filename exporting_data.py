import pickle

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from funcs import access_het_dist

# pd.set_option('display.max_rows', None)

# def select_nodes(df, num_nodes=1):
#     """
#     Selects a specified number of random node indices for each replica.
#
#     :param df: DataFrame containing the heterozygosity data (dat[2]).
#     :param num_nodes: Number of nodes to select per replica.
#     :return: A dictionary with replicas as keys and lists of selected node indices as values.
#     """
#     selection_dict = {}
#     for replica in df['replica'].unique():
#         df_replica = df[df['replica'] == replica]
#         nodes_per_replica = np.argmax(df_replica['step'].to_numpy()[1:] != df_replica['step'].to_numpy()[:-1]) + 1
#         random_indices = np.random.choice(nodes_per_replica, min(num_nodes, nodes_per_replica), replace=False)
#         selection_dict[replica] = random_indices
#     return selection_dict
#
#
# def extract_selected_nodes(df):
#     """
#     Extracts rows for the selected nodes across all steps for each replica.
#
#     :param df: DataFrame containing heterozygosity data.
#     :param selection_dict: A dictionary with replicas as keys and lists of selected node indices as values.
#     :return: DataFrame with the extracted rows, including a node_number column.
#     """
#
#     selection_dict = select_nodes(df, num_nodes=1)
#     selected_rows_across_replicas = pd.DataFrame()
#     for replica, indices in selection_dict.items():
#         df_replica = df[df['replica'] == replica]
#         nodes_per_replica = np.argmax(df_replica['step'].to_numpy()[1:] != df_replica['step'].to_numpy()[:-1]) + 1
#         for index in indices:
#             node_rows = df_replica.iloc[index::nodes_per_replica].copy()
#             node_rows['node_number'] = index + 1  # Assign node number
#             selected_rows_across_replicas = pd.concat([selected_rows_across_replicas, node_rows], ignore_index=True)
#     return selected_rows_across_replicas


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



def get_focal_step(df, num_nodes=10):
    """
    Gets the maximum step in each replica where there are at least 10 'het' values greater than 0.02.

    :param df: DataFrame with 'replica', 'step', and 'het' columns.
    :return: Dictionary with replicas as keys and maximum steps as values.
    """
    # Group by 'replica' and 'step', filter 'het' values greater than 0.02, and count them
    het_counts = df[df['het'] > 0.02].groupby(['replica', 'step']).size()

    # Filter the groups where the count is at least 10
    het_counts = het_counts[het_counts >= num_nodes]

    # For each replica, get the maximum step that satisfies the condition
    max_steps = het_counts.reset_index(level='step').groupby('replica')['step'].max()

    return max_steps.to_dict()



def get_max_het_nodes(df, num_nodes=10):
    """
    Gets the 10 nodes in each replica that have the highest 'het' in the 10th last step of each replica.

    :param df: DataFrame with 'replica', 'step', 'het', and 'node_number' columns.
    :return: Dictionary with replicas as keys and lists of node numbers with the highest 'het' in the 10th last step as values.
    """
    # Get the last step of each replica with n surviving nodes
    top_nodes = {}
    focal_step = get_focal_step(df, num_nodes=num_nodes)

    for replica, step in focal_step.items():
        # Select only the rows that belong to the 10th last step of the current replica
        df_step = df[(df['replica'] == replica) & (df['step'] == step)]

        # Sort the DataFrame by 'het' in descending order and get the 'node_number' of the top 10 rows
        top_nodes_replica = df_step.sort_values(by='het', ascending=False)['node_number'].head(num_nodes).tolist()
        top_nodes[replica] = top_nodes_replica

    return top_nodes

# this is functio nto get only a single node-the last survivor
# def extract_steps_for_nodes(df, max_het_nodes):
#     """
#     Extracts all steps for each node in its corresponding replica from the DataFrame.
#
#     :param df: DataFrame with 'replica', 'step', 'het', and 'node_number' columns.
#     :param max_het_nodes: Dictionary with replicas as keys and node numbers with the highest 'het' in the last step as values.
#     :return: DataFrame with the extracted rows.
#     """
#     extracted_rows = pd.DataFrame()
#     for replica, node_number in max_het_nodes.items():
#         df_replica_node = df[(df['replica'] == replica) & (df['node_number'] == node_number)]
#         extracted_rows = pd.concat([extracted_rows, df_replica_node], ignore_index=True)
#     return extracted_rows


def extract_steps_for_nodes(df, max_het_nodes):
    """
    Extracts all steps for each node in its corresponding replica from the DataFrame.

    :param df: DataFrame with 'replica', 'step', 'het', and 'node_number' columns.
    :param max_het_nodes: Dictionary with replicas as keys and node numbers with the highest 'het' in the last step as values.
    :return: DataFrame with the extracted rows.
    """

    extracted_rows = pd.DataFrame()
    for replica, all_nodes in max_het_nodes.items():

        for node_number in all_nodes:

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
    surviving_nodes = get_max_het_nodes(data, num_nodes=10)
    final_df = extract_steps_for_nodes(data, surviving_nodes)
    final_df.to_csv(f'{frag}_het.csv')
    print('File saved successfully')
    return final_df


frag = 'dist'
with open(f'{frag}_06.pickle', 'rb') as file:
    rand = pickle.load(file)
print('finish')

#plot het for a specific node
# df = assign_node_numbers(rand)
# df = df[(df['replica'] == 99) & (df['node_number'] == 40)]
# plt.plot(df['step'], df['het'])
# plt.show()

# export_het_csv(rand, frag)


