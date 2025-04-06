import pickle
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from statsmodels.tsa.ar_model import AutoReg

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
    data = assign_node_numbers(data)
    surviving_nodes = get_max_het_nodes(data, num_nodes=10)
    final_df = extract_steps_for_nodes(data, surviving_nodes)
    final_df.to_csv(f'{frag}_het.csv')
    print('File saved successfully')
    return final_df










####################################################
######################### analysis #################

######## read the file with RGG (d-0.6) data

# with open(f'cor_d0.6_r1000.pickle', 'rb') as file:
#     cor = pickle.load(file)
# print('finish')
#
# frag = 'cor_d0.6_r1000'
########## plot het for a specific node
# df = assign_node_numbers(cor)
# df = df[(df['replica'] == 99) & (df['node_number'] == 40)]
# plt.plot(df['step'], df['het'])
# plt.show()
# #
# export_het_csv(cor, frag)



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


# load the data
nets = access_networks(cor)
het = assign_node_numbers(cor)
components = get_largest_component(nets)

# get the heterozygosity data for the corresponding largest component
component_data = pd.merge(het, components, on=['replica', 'step', 'node_number'])
component_data = component_data.sort_values(by=['replica', 'step', 'node_number'])
component_data = component_data[component_data['replica'].between(0, 2)]

# plot the heterozygosity data of all nodes of replica 0 with step on x axis and het as y axis
# df = component_data[(component_data['replica'] == 0) & (component_data['node_number'] == 43)]
# # print(df)
# plt.plot(df['step'], df['het'])
# plt.show()

# steps_to_plot = [0, 300, 600, 700]
# for step in steps_to_plot:
#     df = component_data[(component_data['replica'] == 0) & (component_data['step'] == step)]
#     plt.figure()
#     plt.hist(df['het'], bins=15, alpha=0.7)
#     plt.title(f'Histogram of het at step {step}')
#     plt.xlabel('het')
#     plt.ylabel('Frequency')
#     plt.show()
####### calculate the sd skewness and kurtosis for the heterozygosity data


def calculate_return_rate(df):
    """
    Calculate return rates (inverse of AR(1) coefficients) for each 'step'.
    :param df: DataFrame with 'step', 'het', and 'net' columns.
    """
    ar1_coefficients, return_rates, time_indices = [], [], []

    for step, group in df.groupby('step'):
        window_data = group['het']
        window_data_demeaned = window_data - np.mean(window_data)  # Demean the data

        if len(window_data) < 3:
            continue

        # Fit an AR(1) model with no intercept
        model = AutoReg(window_data_demeaned, lags=1, trend='n').fit()
        ar1 = model.params[0]
        return_rate = 1 / ar1

        ar1_coefficients.append(ar1)
        return_rates.append(return_rate)
        time_indices.append(group['step'].iloc[-1])

    # Create a DataFrame to store the results
    results_df = pd.DataFrame({
        'step': time_indices,
        'returnrate': return_rates
    })

    return results_df

#
# rr = calculate_return_rate(x)
# print(rr)
#
# plt.plot(rr['step'], rr['returnrate'])
# plt.ylim(-10, 10)
# plt.show()


def calculate_indicators(data):
    """
    Calculate the standard deviation, skewness, and kurtosis of the heterozygosity data.
    :param data: DataFrame with 'replica', 'step', and 'het' columns.
    :return: DataFrame with the calculated indicators for each step and replica.
    """

    grouped = data.groupby(['replica', 'step'])['het']
    indicators = grouped.agg(['std', 'skew']).reset_index()
    indicators['kurt'] = grouped.apply(pd.Series.kurtosis).values

    # Calculate return rates
    # rr_df = calculate_return_rate(data)
    # print(rr_df)
    # indicators = pd.merge(indicators, rr_df, on='step', how='left')

    return indicators


# indicators = calculate_indicators(component_data)

# save csv
# indicators.to_csv(f'indicators.csv', index=False)
# indicators = pd.read_csv(f'indicators.csv')
# df = indicators[(indicators['replica'] == 5)]
# print(df)
# plt.plot(df['step'], df['skew'])
# plt.show()



# summary = indicators.groupby('step').agg(
#     mean_std=('std', 'mean'),
#     mean_skew=('skew', 'mean'),
#     mean_kurt=('kurt', 'mean'),
#     ci95_std=('std', lambda x: 1.96 * x.std() / (len(x)**0.5)),
#     ci95_skew=('skew', lambda x: 1.96 * x.std() / (len(x)**0.5)),
#     ci95_kurt=('kurt', lambda x: 1.96 * x.std() / (len(x)**0.5))
# ).reset_index()
#
# print(summary)
#
#
# plt.plot(summary['step'], summary['mean_std'])
# plt.fill_between(summary['step'],
#                  summary['mean_std'] - summary['ci95_std'],
#                  summary['mean_std'] + summary['ci95_std'],
#                   alpha=0.2)
# plt.xlabel('Step')
# plt.ylabel('Standard Deviation')
# plt.savefig(f'./figs/ews_std.svg', format="svg")
# plt.show()
#
# plt.plot(summary['step'], summary['mean_skew'])
# plt.fill_between(summary['step'],
#                  summary['mean_skew'] - summary['ci95_skew'],
#                  summary['mean_skew'] + summary['ci95_skew'],
#                   alpha=0.2)
# plt.xlabel('Step')
# plt.ylabel('skew')
# plt.savefig(f'./figs/ews_skew.svg', format="svg")
# plt.show()
#
# plt.plot(summary['step'], summary['mean_kurt'])
# plt.fill_between(summary['step'],
#                  summary['mean_kurt'] - summary['ci95_kurt'],
#                  summary['mean_kurt'] + summary['ci95_kurt'],
#                   alpha=0.2)
# plt.xlabel('Step')
# plt.ylabel('kurt')
# plt.savefig(f'./figs/ews_kurt.svg', format="svg")
# plt.show()

