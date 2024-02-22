from statistics import mean
import pickle
import pandas as pd
from joypy import joyplot
from matplotlib import pyplot as plt

from funcs_analysis import load_data, plot_data, filter_intervals

# from funcs3 import make_networks, make_replicates_new

# pd.set_option('display.max_rows',None)


########### run pipeline
# print("here i start!")
# n = 50  # no. of nodes
# n_rep = 100
# net = "RGG"
# ignore = False
#
# # # create list off nets
# nets = make_networks(n_nets=n_rep, n_nodes=n, net_type=net)
#
# # run the pipeline for all fragmentation types
# rand = make_replicates_new(nets=nets, frag_type='rand', ignore=ignore)
# print("1")
# cor = make_replicates_new(nets=nets, frag_type='cor', ignore=ignore)
# print("2")
# int = make_replicates_new(nets=nets, frag_type='int', ignore=ignore)
# print("3")
# reg = make_replicates_new(nets=nets, frag_type='reg', ignore=ignore)
# print("4")
# div = make_replicates_new(nets=nets, frag_type='div', ignore=ignore)
# print("5")
# dist = make_replicates_new(nets=nets, frag_type='dist', ignore=ignore)
# print("6")
# dist = make_replicates_new(nets=nets, frag_type='opt', ignore=ignore)
# print("7")
#
## save files as tuple
# pickle_filename = f'{net}, rand_ignore_{ignore}.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(rand, file)
#
# pickle_filename = f'{net}, cor_ignore_{ignore}.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(cor, file)
#
# pickle_filename = f'{net}, int_ignore_{ignore}.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(int, file)
#
# pickle_filename = f'{net}, reg_ignore_{ignore}.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(reg, file)
#
# pickle_filename = f'{net}, div_ignore_{ignore}.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(div, file)
#
# pickle_filename = f'{net}, dist_ignore_{ignore}.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(dist, file)
#
# pickle_filename = f'{net}, opt_ignore_{ignore}.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(opt, file)
########### finish pipeline



fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)

# Plot fst and het along fragmentation
plot_data(data, 5, 'Pairwise Fst',measure='fst', save=True)
plot_data(data, 3, 'Heterozygosity',measure='heterozygosity', save=True)

##############plot distributions


def plot_fragmentation_types(data):
    """
    Plots density plots for each fragmentation type separately in a single figure.

    :param data: Dictionary containing data for all fragmentation types.
    """
    # Determine the number of fragmentation types to decide the number of subplots
    num_types = len(data)
    fig, axes = plt.subplots(num_types, 1, figsize=(8, 4 * num_types))  # Adjust size as needed

    for i, (frag_type, frag_data) in enumerate(data.items()):
        df = filter_intervals(frag_data[4])  # Assuming this returns a DataFrame

        # Plotting each density plot separately
        joyplot(
            data=df[['fst', 'step']],
            by='step', overlap=1,
            colormap=plt.cm.viridis, fade=True, range_style='all',
            linecolor="black", linewidth=0.1,
            ax=axes[i]  # This specifies which subplot to use
        )

        axes[i].set_title(f'{frag_type.capitalize()} Fragmentation', fontsize=18)  # Customize title as needed
        axes[i].tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.show()


plot_fragmentation_types(data)



# with open('RGG, opt_ignore_True.pickle', 'rb') as file:
#     rand = pickle.load(file)
# # First plot for 'fst'
# df = filter_intervals(rand[4])
#
# plt.figure(figsize=(8, 8))
# fig, axes = joyplot(
#     data=df[['fst', 'step']],
#     by='step', overlap=1,
#     colormap=plt.cm.viridis, fade=True, range_style='all',
#     linecolor="black", linewidth=0.1
# )
#
# fig.suptitle('Fst - distance fragmentation', fontsize=18)
# for ax in axes:
#     ax.tick_params(axis='both', which='major', labelsize=16)
#
# plt.show()
#
#
# # Second plot for 'het'
# df = filter_intervals(rand[2])
#
# plt.figure(figsize=(8, 8))
# fig, axes = joyplot(
#     data=df[['het', 'step']],
#     by='step', overlap=2,
#     colormap=plt.cm.viridis, fade=True, range_style='all',
#     linecolor="black", linewidth=0.1)
#
# fig.suptitle('Heterozygosity - distance fragmentation', fontsize=18)
# for ax in axes:
#     ax.tick_params(axis='both', which='major', labelsize=16)
#     ax.set_xlim(-0.1, 1.3)
# plt.show()




#
# def ignore_isolated(df:pd.DataFrame):
#     """
#     Filter the DataFrame to include only non-isolated populations
#     in fst remove fst=1
#     in heterozygosity remove het=0.02
#     :param df: dataframe with fst and heterozygosity distributions
#     :return: filtered dataframe
#     """
#     if 'het' in df.columns:
#         df_filtered = df[df['het'] != 0.02]
#     if 'fst' in df.columns:
#         df_filtered = df[df['fst'] != 1]
#
#     return df_filtered
#
#
#
# def make_het_stat(df: pd.DataFrame) -> pd.DataFrame:
#     """
#     Calculate the mean heterozygosity of each step for each replica.
#
#     :param df: DataFrame of heterozygosity distribution in each step and replica
#     :return: DataFrame of average heterozygosity for each step and replica
#     """
#
#     if 'het' in df.columns:
#         group_means = df.groupby(['step', 'replica'])['het'].mean().reset_index()
#         group_means = group_means.rename(columns={'het': 'avg'})
#
#     if 'fst' in df.columns:
#         group_means = df.groupby(['step', 'replica'])['fst'].mean().reset_index()
#         group_means = group_means.rename(columns={'fst': 'avg'})
#     print(group_means)
#     return group_means
#
#
# def plot_data_from_list(dataframes, column, labels, ylabel, title):
#     """
#     Plot data from a list of dataframes.
#
#     :param dataframes: List of dataframes to be plotted.
#     :param column: The column name to be used for plotting.
#     :param labels: List of labels for each dataframe.
#     :param ylabel: Label for the Y-axis.
#     :param title: Title of the plot.
#     :param filename: Filename for saving the plot.
#     """
#     plt.figure()
#
#     # Iterate over each dataframe and its corresponding label
#     for df, label in zip(dataframes, labels):
#         mean_values = df.groupby('step')['avg'].mean()
#         confidence = df.groupby('step')['avg'].std()
#
#         plt.plot(mean_values, label=label)
#         plt.fill_between(mean_values.index, mean_values - confidence, mean_values + confidence, alpha=0.2)
#
#     plt.xlabel('Fragmentation step')
#     plt.ylabel(ylabel)
#     plt.title(title)
#     plt.legend()
#     # plt.savefig(f"{filename}.svg", format="svg")
#     # plt.close()
#     plt.show()  # Uncomment if you want to display the plot as well
#
#
# df_items = list(data.values())
# df_items = [item[4] for item in df_items]
# filtered_dfs = [ignore_isolated(df) for df in df_items]
#
# x = [make_het_stat(df) for df in filtered_dfs]
#
# plot_data_from_list(x, 'avg',
#                     fragmentation_types,
#                     'Y-axis Label', 'Plot Title')

