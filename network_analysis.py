import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from funcs import load_data, normalize_steps, calculate_statistics, compute_modularity, calculate_centrality, measure_giant_component


def plot_het_central(data: dict, measure: str):
    fragmentation_types = list(data.keys())
    plt.figure()
    color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette

    for i, frag_type in enumerate(fragmentation_types):
        het = data[frag_type][3]
        central = calculate_centrality(data[frag_type][1], measure=measure)
        merged = pd.merge(het, central, how='outer')
        merged = merged[merged[measure] != 0]

        if measure == 'component':
            sns.regplot(x='component', y='avg', data=merged, fit_reg=True, order=2,
                        truncate=True, scatter_kws={'s': 50, 'alpha': 0.01, 'color': color_palette(i)},
                        line_kws={'lw': 2, 'label': frag_type})

            # add a diagonal line
            plt.plot([0.05, 1], [0.05, 1], linestyle='--', color='black',linewidth=1)
            plt.xlabel('Fraction of nodes in the largest component', fontsize=16)


        if measure == 'modularity':
            sns.regplot(x='modularity', y='avg', data=merged, fit_reg=True, order=3,
                        truncate=True, scatter_kws={'s': 50, 'alpha': 0.01, 'color': color_palette(i)},
                        line_kws={'lw': 2, 'label': frag_type})

            plt.gca().invert_xaxis()
            plt.ylim(-0.05, 1.05)
            plt.xlabel('Modularity', fontsize=16)

    plt.ylabel('Heterozygosity', fontsize=16)
    plt.legend()

    plt.savefig(f'./figs/het_{measure}.jpg', format="jpg")
    plt.show()


################################
################################ stack plot
def measure_isolated_nodes(network: nx.Graph) -> int:
    """
    Measure the number of isolated nodes in the network.
    :param network: NetworkX graph
    :return: Number of isolated nodes
    """
    isolated_nodes = list(nx.isolates(network))
    return len(isolated_nodes) / len(network)

def measure_components(network: nx.Graph, min_size: int = 4) -> int:
    """
    Measure the number of components with a size greater than or equal to a given threshold,
    excluding the giant component.
    :param network: NetworkX graph
    :param min_size: Minimum size of components to be counted
    :return: Number of nodes in large components excluding the giant component
    """
    largest_component = max(nx.connected_components(network), key=len)

    components = [
        comp for comp in nx.connected_components(network)
        if (comp != largest_component or len(comp) == min_size) and len(comp) >= min_size
    ]

    return sum(len(comp) for comp in components) / len(network)


def measure_waste(network: nx.Graph, max_size: int = 3, min_size: int = 2) -> int:
    """
    Measure the number of components with a size greater than or equal to a given threshold,
    excluding the giant component.
    :param network: NetworkX graph
    :param min_size: Minimum size of components to be counted
    :return: Number of nodes in large components excluding the giant component
    """
    components = [comp for comp in nx.connected_components(network) if min_size <= len(comp) <= max_size]
    num_nodes_in_medium_components = sum(len(comp) for comp in components)
    return num_nodes_in_medium_components / len(network)


def measure_network_metrics(networks: list) -> pd.DataFrame:
    """
    Measure various metrics of the networks and return them as a DataFrame:
    - Size of the giant component
    - Number of isolated nodes
    - Number of components with 4 or more nodes excluding the giant component
    :param networks: List of NetworkX graphs
    :return: DataFrame with metrics for each network
    """
    metrics = []

    for step, network in enumerate(networks):
        giant_component = measure_giant_component(network)
        isolated_nodes = measure_isolated_nodes(network)
        components = measure_components(network)
        waste = measure_waste(network)

        total = giant_component + isolated_nodes + components + waste

        # Round the first three metrics
        giant = round(giant_component / total, 2)
        isolated = round(isolated_nodes / total, 2)
        components = round(components / total, 2)

        # Adjust the last metric so the total sums up to 1
        waste = 1 - giant - isolated - components

        scaled_metrics = {
            "step": step,
            "giant": giant,
            "isolated": isolated,
            "components": components,
            "waste": waste,
        }

        metrics.append(scaled_metrics)

    return pd.DataFrame(metrics)


def measure_network_metrics_replicas(replicas: list) -> pd.DataFrame:
    """
    Measure metrics for a list of lists of networks (replicas) and return a DataFrame
    including a column for the replica index.
    :param replicas: List of lists of NetworkX graphs
    :return: DataFrame with metrics for each network and replica
    """
    all_metrics = []

    for replica_index, networks in enumerate(replicas):
        replica_metrics = measure_network_metrics(networks)
        replica_metrics['replica'] = replica_index
        all_metrics.append(replica_metrics)

    return pd.concat(all_metrics, ignore_index=True)



def calculate_statistics(df):
    """Calculate mean and 95% confidence interval for all columns in the dataframe."""
    result = []

    # Select all columns except 'step' and 'replica'
    columns_to_analyze = df.columns.difference(['step', 'replica'])

    for column in columns_to_analyze:
        mean_values = round(df.groupby('step')[column].mean(),3)

        # Create a DataFrame for this column's statistics
        column_stats = pd.DataFrame({
            'step': mean_values.index,
            f'{column}': mean_values.values,
            # f'{column}_ci': confidence_interval.values
        })

        result.append(column_stats)

    # Concatenate all column statistics DataFrames along the 'step' index
    result_df = pd.concat(result, axis=1)

    # Remove duplicate 'step' columns
    result_df = result_df.loc[:, ~result_df.columns.duplicated()]
    # result_df['waste'] = 1 - result_df['giant'] - result_df['isolated'] - result_df['components']

    return result_df


def plot_network_stacked_area(df: pd.DataFrame, frag: str):
    """
    Plot the metrics as stacked area charts.
    :param df: DataFrame containing the metrics to plot
    :param frag: Fragmentation type
    """

    # Ensure the DataFrame is sorted by 'step'
    df = df.sort_values(by='step')
    df = normalize_steps(df)

    # Create a new figure and axes with a specific size
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define the columns to plot and the colors to use
    columns = ['waste', 'isolated', 'components', 'giant']
    colors = plt.cm.Dark2.colors[:len(columns)]

    # Prepare the data for the stackplot
    x_values = df['step'].values
    y_values = [df[col].values for col in columns]

    # Create the stackplot
    ax.stackplot(x_values, y_values, labels=columns, colors=colors, alpha=0.8)

    # Set parameters
    ax.set_xlabel('Fragmentation (%)', fontsize=20)
    ax.set_ylabel('Proportion of the network', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=18)  # Increase the size of the tick labels
    ax.set_title(frag)
    plt.ylim(0, 1)
    # plt.legend(loc='upper left')

    plt.savefig(f'./figs/stack_{frag}.jpg')
    plt.show()
#
#
# def calculate_mantel_correlation(data):
#     # Initialize a list to store the results
#     mantel_results = []
#
#     # Get the number of replicas
#     num_replicas = len(data[7])
#
#     # Loop over all replicas
#     for rep in range(num_replicas):
#         # Get the number of steps for the current replica
#         num_steps = len(data[7][rep])
#         num_steps = 10
#
#         # Loop over all steps
#         for step in range(num_steps):
#             # Get the FST matrix and network for the current step
#             matrix = data[7][rep][step]
#             net = data[1][rep][step]
#
#             distance_matrix = get_euclidean_matrix(net)
#
#             correlation = perform_mantel_test(matrix, distance_matrix, perms=999,
#                                               method='pearson',print=False)[0]
#
#             mantel_results.append({
#                 'replica': rep,
#                 'step': step,
#                 'correlation': correlation
#             })
#
#     return pd.DataFrame(mantel_results)
#
#
# def plot_step_vs_correlation(df):
#     plt.figure(figsize=(10, 6))
#     plt.plot(df['step'], df['correlation'])
#     plt.xlabel('Step')
#     plt.ylabel('Correlation')
#     plt.show()
#
#
# def calculate_statistics(df):
#     # Group the data by 'step'
#     grouped = df.groupby('step')
#
#     # Calculate the mean correlation for each step
#     mean_correlation = grouped['correlation'].mean()
#
#     # Calculate the standard error of the mean for each step
#     sem_correlation = grouped['correlation'].sem()
#
#     ci_correlation = 1.96 * sem_correlation  # 95% confidence interval
#
#     # Create a new DataFrame for the results
#     results = pd.DataFrame({
#         'step': mean_correlation.index,
#         'correlation_mean': mean_correlation.values,
#         'correlation_ci': ci_correlation.values
#     })
#
#     return results
#
# def plot_correlation_with_ci(df):
#     plt.figure(figsize=(10, 6))
#     sns.lineplot(x='step', y='correlation_mean', data=df)
#     plt.fill_between(df['step'], df['correlation_mean'] - df['correlation_ci'], df['correlation_mean'] + df['correlation_ci'], color='b', alpha=0.1)
#     plt.xlabel('Step')
#     plt.ylabel('Correlation')
#     plt.show()
#
#
# def plot_correlation_with_ci(df, fragmentation_types):
#     plt.figure(figsize=(10, 6))
#
#     # If a single fragmentation type is provided, convert it to a list
#     if isinstance(fragmentation_types, str):
#         fragmentation_types = [fragmentation_types]
#
#     # Loop over all fragmentation types
#     for frag_type in fragmentation_types:
#         # Filter the data for the current fragmentation type
#         df_filtered = df[df['fragmentation_type'] == frag_type]
#
#         # Plot the data for the current fragmentation type
#         sns.lineplot(x='step', y='correlation_mean', data=df_filtered, label=frag_type)
#         plt.fill_between(df_filtered['step'], df_filtered['correlation_mean'] - df_filtered['correlation_ci'], df_filtered['correlation_mean'] + df_filtered['correlation_ci'], alpha=0.1)
#
#     plt.xlabel('Step')
#     plt.ylabel('Correlation')
#     plt.legend()
#     plt.show()
#
#
# def process_multiple_fragmentation_types(data, fragmentation_types):
#     for frag_type in fragmentation_types:
#         frag_data = data[frag_type]
#         correlations = calculate_mantel_correlation(frag_data)
#         stats = calculate_statistics(correlations)
#         plot_correlation_with_ci(stats,fragmentation_types=frag_type)




###########################################@@#####################
###############################  analysis  #######################

# fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt']
# fragmentation_types = ['opt']
# net = 'RGG'
# ignore = False
# data = load_data(fragmentation_types, net, ignore)

#########################plot centrality vs heterozygosity
# plot_het_central(data, measure='modularity')


##############################plot stacks
# networks = data[fragmentation_types[0]][1]
# matrices = measure_network_metrics_replicas(networks)
# stats = calculate_statistics(matrices)
# plot_network_stacked_area(stats,frag=fragmentation_types[0])

##############################plot centrality vs fragmnetation
# plot_centrality(data,centrality='connectivity')

# x=calculate_mantel_correlation(data)
# print(x)
# stats= calculate_statistics(x)
# print(stats)
# plot_correlation_with_ci(stats, fragmentation_types)

##plot fst-distance relationship
# data = data[frag]
# matrix = data[7][0][20]
# net = data[1][0][20]
# distance_matrix = get_euclidean_matrix(net)
# perform_mantel_test(matrix,distance_matrix)

# distance_matrix = get_distance_matrix(net)
# plot_matrix_relationship(distance_matrix=distance_matrix,fst_matrix=matrix)


