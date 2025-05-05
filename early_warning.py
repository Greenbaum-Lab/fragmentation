import pickle
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from statsmodels.tsa.ar_model import AutoReg
import random
from scipy.stats import skew, kurtosis

from funcs import access_het_dist, access_networks, normalize_steps


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
            node_rows['node_number'] = index  # Assign node number
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

######## read the pickle file with RGG (d-0.6) data

# with open(f'cor_d0.6_r1000.pickle', 'rb') as file:
#     cor = pickle.load(file)
# print('finish')
# frag = 'cor_d0.6_r1000'
#save the heterozygosity data to a csv file
# export_het_csv(cor, frag)


#plot nodes
# cor= pd.read_csv('cor_d0.6_r1000_het.csv')
# print(cor)
# plot_random_nodes(cor, num_nodes=10)

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
# nets = access_networks(cor)
# het = assign_node_numbers(cor)
# components = get_largest_component(nets)
#
# # get the heterozygosity data for the corresponding largest component
# component_data = pd.merge(het, components, on=['replica', 'step', 'node_number'])
# component_data = component_data.sort_values(by=['replica', 'step', 'node_number'])
# component_data = component_data[component_data['replica'].between(0, 2)]

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
    indicators.to_csv(f'indicators.csv', index=False)

    # Calculate return rates
    # rr_df = calculate_return_rate(data)
    # print(rr_df)
    # indicators = pd.merge(indicators, rr_df, on='step', how='left')

    return indicators


# indicators = calculate_indicators(component_data)


#### final plot


def plot_het_indicator(cor: pd.DataFrame,
                       indicators: pd.DataFrame,
                       indicator: str,
                       n_samples: int = 10) -> None:
    """
    Creates a combined plot with two y-axes:
      - Left y-axis: overall mean and standard deviation (SD) of heterozygosity (het) across all replicates,
        plus n_samples random individual replica curves for het.
      - Right y-axis: overall mean and SD of the specified indicator across all replicates,
        plus n_samples random individual replica curves for that indicator.

    The x-axis shows the 'step' converted to a percentage (% fragmentation).

    Parameters:
    -----------
    cor : pd.DataFrame
        DataFrame containing columns: ['step', 'het', 'replica', 'node_number'].
    indicators : pd.DataFrame
        DataFrame containing columns: ['replica', 'step', 'std', 'skew', 'kurt'].
    indicator : str
        The indicator to plot (e.g. 'skew', 'std', or 'kurt').
    n_samples : int, optional
        Number of random individual replica curves to plot for each metric. Default is 10.

    Returns:
    --------
    None
        Displays the combined plot.
    """
    # Use a common scaling for the x-axis: convert 'step' to % fragmentation.
    global_max = max(cor['step'].max(), indicators['step'].max())

    # ---------------------------
    # Overall Heterozygosity Stats
    # ---------------------------
    stats_het = cor.groupby('step')['het'].agg(mean='mean', std='std').reset_index()
    stats_het['step_pct'] = stats_het['step'] / global_max * 100

    # ---------------------------
    # Overall Indicator Stats
    # ---------------------------
    stats_ind = indicators.groupby('step')[indicator].agg(mean='mean', std='std').reset_index()
    stats_ind['step_pct'] = stats_ind['step'] / global_max * 100

    # ---------------------------
    # Set up the figure with two y-axes
    # ---------------------------
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax1 = plt.subplots(figsize=(12, 8))
    ax2 = ax1.twinx()
    color_het = 'darkorange'
    color_ind = 'forestgreen'

    # Plot overall het mean and SD on left axis (ax1)
    ax1.plot(stats_het['step_pct'], stats_het['mean'],
             color=color_het)
    ax1.fill_between(stats_het['step_pct'],
                     stats_het['mean'] - stats_het['std'],
                     stats_het['mean'] + stats_het['std'],
                     color=color_het, alpha=0.2)

    # Plot overall indicator mean and SD on right axis (ax2)
    ax2.plot(stats_ind['step_pct'], stats_ind['mean'],
             color=color_ind)
    ax2.fill_between(stats_ind['step_pct'],
                     stats_ind['mean'] - stats_ind['std'],
                     stats_ind['mean'] + stats_ind['std'],
                     color=color_ind, alpha=0.2)

    # ---------------------------
    # Plot individual sample curves for both het and indicator
    # ---------------------------
    # Identify replicas present in both DataFrames.
    common_replicas = list(set(cor['replica'].unique()).intersection(set(indicators['replica'].unique())))
    selected_replicas = random.sample(common_replicas, n_samples)

    for i, replica in enumerate(selected_replicas):
        # --- Heterozygosity for this replica ---
        rep_cor = cor[cor['replica'] == replica]
        rep_het = rep_cor.groupby('step', as_index=False)['het'].mean()
        rep_het['step_pct'] = rep_het['step'] / global_max * 100

        # --- Indicator for this replica ---
        rep_ind = indicators[indicators['replica'] == replica].sort_values('step')
        rep_ind['step_pct'] = rep_ind['step'] / global_max * 100

        # Plot sample curves with lower opacity.
        ax1.plot(rep_het['step_pct'], rep_het['het'], color=color_het, alpha=0.5)
        ax2.plot(rep_ind['step_pct'], rep_ind[indicator], color=color_ind, alpha=0.5)

    ax1.set_xlabel('% fragmentation', fontsize=32)
    ax1.set_ylabel('Heterozygosity', color=color_het, fontsize=32)
    ax2.set_ylabel("Kurtosis", color=color_ind, fontsize=32)

    ax1.tick_params(axis='y', labelsize=28, labelcolor=color_het)
    ax2.tick_params(axis='y', labelsize=28, labelcolor=color_ind)
    ax1.tick_params(axis='x', labelsize=28)

    plt.savefig(f'./figs/het_{indicator}_singlepop.svg', format='svg')
    plt.show()


# plot het+indicators of metapopilation data-change y label
# set random seed for reproducibility
# random.seed(1)
# cor = pd.read_csv('cor_d0.6_r1000_het.csv')
# indicators = pd.read_csv('indicators.csv')
# plot_het_indicator(cor, indicators, indicator='std', n_samples=10)

# plot het+indicators of single population data-change y label
# set random seed for reproducibility
# random.seed(5)
# cor = pd.read_csv('cor_d0.6_r1000_het.csv')
# indicators = pd.read_csv('indicators_singlepop.csv')

# plot_het_indicator(cor, indicators, indicator='kurt', n_samples=10)



with open(f'cor_d0.6_r1000.pickle', 'rb') as file:
    cor = pickle.load(file)
print('finish')
nets = access_networks(cor)
het = assign_node_numbers(cor)
components = get_largest_component(nets)
# get the heterozygosity data for the corresponding largest component
component_data = pd.merge(het, components, on=['replica', 'step', 'node_number'])
component_data = component_data.sort_values(by=['replica', 'step', 'node_number'])
print(component_data)
mean_het_per_step_replica = component_data.groupby(['replica', 'step'])['het'].mean()

# Save the resulting DataFrame to a CSV file
mean_het_per_step_replica.reset_index().to_csv('cor_d0.6_r1000_mean_het.csv', index=False)
print(mean_het_per_step_replica.reset_index())

######## truncate data before tipping point
# cor = pd.read_csv('cor_d0.6_r1000_het.csv')
# Initialize a list to store the results
# results = []
# for replica in cor['replica'].unique():
#     # Get the data for the current replica
#     replica_data = cor[cor['replica'] == replica]
#     # Loop over each unique node_number for the current replica
#     for node in replica_data['node_number'].unique():
#         # Get the data for the current node in the replica
#         node_data = replica_data[replica_data['node_number'] == node]
#         # Extract step and het values
#         time = node_data['step']
#         het = node_data['het']
#         # Compute the first derivative (slope) and second derivative (acceleration)
#         slope = np.gradient(het, time)
#         acceleration = np.gradient(slope, time)
#         # Find the step where acceleration is minimum
#         min_accel_index = np.argmin(acceleration)
#         min_accel_step = time.iloc[min_accel_index]  # Get the corresponding step
#         # Store the result in the results list
#         results.append({
#             'replica': replica,
#             'node_number': node,
#             'min_acceleration_step': min_accel_step
#         })
#
# results_df = pd.DataFrame(results)
# cor = pd.merge(cor, results_df, on=['replica', 'node_number'], how='left')
# cor_truncated = cor[cor['step'] < cor['min_acceleration_step']]
# cor_truncated = cor_truncated.drop(columns=['min_acceleration_step'])
# cor_truncated.to_csv('test.csv', index=False)




########## plot distribution of indicators of single pop
# indicators = pd.read_csv('kendall_singlepop_truncated.csv')
# indicator_columns = ['kurt.tau', 'kurt.p', 'sk.tau','sk.p', 'sd.tau', 'sd.p', 'returnrate.tau', 'returnrate.p']
#
# # Create a figure with 2 columns and 4 rows
# fig, axes = plt.subplots(4, 2, figsize=(15, 20))
# axes = axes.flatten()
#
# # Iterate over each indicator column and plot
# for i, col in enumerate(indicator_columns):
#     ax = axes[i]
#     ax.hist(indicators[col], bins=30, color='skyblue', edgecolor='black', alpha=0.7)
#     median_value = indicators[col].median()
#     ax.axvline(median_value, color='red', linestyle='dashed', linewidth=1.5)
#     # ax.set_xlabel('Kendall tau', fontsize=14)
#     ax.set_ylabel('Frequency', fontsize=14)
#     ax.set_title(f'{col}', fontsize=16)
#     ax.tick_params(axis='both', labelsize=12)
#
# # Adjust layout and show the plots
# plt.tight_layout()
# plt.show()
#

































# Define the logistic function with added noise
# def noisy_logistic_function(x, L=1, k=1, x0=0, noise_level=0.1):
#     logistic = L / (1 + np.exp(-k * (x - x0)))
#     noise = np.random.normal(0, noise_level, size=x.shape)
#     return logistic + noise
#
# x = np.linspace(-10, 10, 400)
# y = noisy_logistic_function(x, L=-1, k=3, x0=1, noise_level=0.1)
# plt.figure(figsize=(10, 6))
# plt.plot(x, y, label='Noisy Logistic Function', color='grey')
# plt.xlabel('Time', fontsize=30)
# plt.ylabel('Genetic diversity', fontsize=30)
# plt.tick_params(axis='both', labelsize=26)
# plt.savefig('./figs/logistic_func.svg', format='svg')
# plt.show()


# Create a Random Geometric Graph with 10 nodes
# num_nodes = 10
# radius = 0.5  # Radius within which nodes are connected
# rgg = nx.random_geometric_graph(num_nodes, radius)
# pos = nx.get_node_attributes(rgg, 'pos')
# nx.draw(rgg, pos, node_size=500)
# plt.savefig('./figs/ews_rgg.svg', format='svg')
# plt.show()




#need 0 skewness and 0 kurtosis and low variance
# Load the DataFrame
cor = pd.read_csv('cor_d0.6_r1000_het.csv')

# Filter the DataFrame for step 600 and the first 10 replicas
# filtered_cor = cor[(cor['step'] == 720) & (cor['replica'] < 100)]





# Plot the distribution of 'het' values
# plt.figure(figsize=(10, 6))
# plt.hist(filtered_cor['het'], bins=30, color='skyblue', edgecolor='black', alpha=0.7)
# plt.ylabel('Frequency', fontsize=30)
# plt.xlabel('Genetic diversity', fontsize=30)
# plt.tick_params(axis='both', labelsize=25)
# plt.savefig('./figs/ews_dist2.svg', format='svg')
# plt.show()
# #calculate skewness and kurtosis
# skewness = filtered_cor['het'].skew()
# kurtosis = filtered_cor['het'].kurtosis()
# sd= filtered_cor['het'].std()
# print("Skewness:", skewness)
# print("Kurtosis:", kurtosis)
# print("Standard Deviation:", sd)




# Parameters for the normal distribution
# mean = 0
# std_dev = 0.5  # Low standard deviation for low variance
# size = 1000  # Number of samples
#
# # Generate the data
# data = np.random.normal(loc=mean, scale=std_dev, size=size)
#
# # Calculate skewness, kurtosis, and standard deviation
# data_skewness = skew(data)
# data_kurtosis = kurtosis(data)  # Excess kurtosis
# data_std = np.std(data)
#
# # Plot the distribution
# plt.figure(figsize=(10, 6))
# plt.hist(data, bins=30, color='skyblue', edgecolor='black', alpha=0.7)
# plt.title(f'Distribution with Skewness: {data_skewness:.2f}, Kurtosis: {data_kurtosis:.2f}, Std Dev: {data_std:.2f}')
# plt.xlabel('Value')
# plt.ylabel('Frequency')
# plt.savefig('./figs/ews_dist1.svg', format='svg')
# plt.show()

