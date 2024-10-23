import pickle

import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from funcs import access_het_dist, access_networks

#show the first 1000 rows
# pd.set_option('display.max_rows', None)


def assign_node_numbers(df, nodes_per_step=50):
    """
    Assigns node numbers for each node in each replica.

    :param df: DataFrame containing the heterozygosity data (dat[2]).
    :param nodes_per_step: Number of nodes in each step. i.e. the number of nodes in the network.
    :return: DataFrame with the original data and an additional 'node_number' column.
    """
    dfs = []  # List to collect DataFrames
    df = access_het_dist(df)
    for replica in df['replica'].unique():
        df_replica = df[df['replica'] == replica]
        for index in range(nodes_per_step):
            node_rows = df_replica.iloc[index::nodes_per_step].copy()
            node_rows['node_number'] = index   # Assign node number
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
    data = assign_node_numbers(cor)
    surviving_nodes = get_max_het_nodes(data, num_nodes=10)
    final_df = extract_steps_for_nodes(data, surviving_nodes)
    final_df.to_csv(f'{frag}_het.csv')
    print('File saved successfully')
    return final_df










####################################################
######################### analysis #################

######## read the file with RGG (d-0.3) data

df = pd.read_csv('cor_het.csv')
time_series = df[(df['replica'] == 0) & (df['node_number'] == 49)]
print(df)


window_size = int(len(df) * 0.5)  # 50% of the data
step_size = 1  # Sliding step (you can adjust this)

def calc_autocorr_lag1(series):
    return series.autocorr(lag=1)

# Store autocorrelation results
autocorr_results = []

# Loop through the data with the sliding window
for start in range(0, len(df) - window_size + 1, step_size):
    window_data = df['het'][start:start + window_size]  # Extract the window of 'het'
    autocorr_value = calc_autocorr_lag1(window_data)  # Calculate lag-1 autocorrelation
    autocorr_results.append({
        'start_index': start,
        'end_index': start + window_size - 1,
        'autocorr_lag1': autocorr_value
    })

# Convert results to a DataFrame for better visualization
autocorr_df = pd.DataFrame(autocorr_results)


# het_series = pd.Series(df['het'].values)
# autocorr_lag1 = het_series.autocorr(lag=1)
print(autocorr_df)

plt.show(autocorr_df)
plt.show()


frag = 'cor'
with open(f'{frag}_06.pickle', 'rb') as file:
    cor = pickle.load(file)
print('finish')


# follow the largest component in the network for each replica across the fragmentation procss
def get_largest_component(nets):
    """
    get the largest component of each network in each replica and step.
    for early warning analysis.
    :param nets: list of lists of networks in the format of nets[replica][step]
    :return: a dataframe with the largest component of each network, the replica and the step
    """

    components = []
    for replica in range(len(nets)):
        for step in range(len(nets[0])):
            net = nets[replica][step]
            largest_component = max(nx.connected_components(net), key=len)
            for node in largest_component:
                components.append({'replica': replica, 'step': step, 'node_number': node})

    return pd.DataFrame(components)


nets = access_networks(cor)
het = assign_node_numbers(cor)
components = get_largest_component(nets)

# get the heterozygosity data for the corresponding largest component
full_data = pd.merge(het, components, on=['replica', 'step', 'node_number'])

# plot the heterozygosity data of all nodes of replica 0 with step on x axis and het as y axis
# df = full_data[(full_data['replica'] == 0) & (full_data['node_number'] == 43)]
# # print(df)
# plt.plot(df['step'], df['het'])
# plt.show()

# calculate the sd skewness and kurtosis for the heterozygosity data
def calculate_indicators(data):
    """
    Calculate the standard deviation, skewness, and kurtosis of the heterozygosity data.
    :param data: DataFrame with 'replica', 'step', and 'het' columns.
    :return: DataFrame with the calculated indicators for each step and replica.
    """
    grouped = data.groupby(['replica', 'step'])['het']
    indicators = grouped.agg(['std', 'skew']).reset_index()
    indicators['kurt'] = grouped.apply(pd.Series.kurtosis).values

    return indicators


indicators = calculate_indicators(full_data)
print(indicators)

df = indicators[(indicators['replica'] == 0)]
# print(df)
plt.plot(df['step'], df['std'])
plt.show()

df = indicators[(indicators['replica'] == 0)]
# print(df)
plt.plot(df['step'], df['skew'])
plt.show()

df = indicators[(indicators['replica'] == 0)]
# print(df)
plt.plot(df['step'], df['kurt'])
plt.show()



# Create a pandas series



#
########## plot het for a specific node
# df = assign_node_numbers(cor)
# df = df[(df['replica'] == 99) & (df['node_number'] == 40)]
# plt.plot(df['step'], df['het'])
# plt.show()
# #
# export_het_csv(rand, frag)


