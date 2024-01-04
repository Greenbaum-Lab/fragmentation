from statistics import mean
import pickle
from joypy import joyplot
from matplotlib import pyplot as plt

# from funcs3 import make_networks, make_replicates_new
from processes import find_breaking_point, find_breakink_point_list



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





color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette



def load_data(fragmentation_types, net, ignore):
    data = {}
    for frag_type in fragmentation_types:
        filename = f'RGG, {frag_type}_ignore_{ignore}.pickle'
        with open(filename, 'rb') as file:
            data[frag_type] = pickle.load(file)
    print("finish load")
    return data

def plot_data(data, index, ylabel, title, net, ignore, filename):
    plt.figure()
    for frag_type, datasets in data.items():
        mean_values = datasets[index].groupby('step')['avg'].mean()
        confidence = datasets[index].groupby('step')['avg'].std()
        plt.plot(mean_values, label=frag_type.capitalize())
        plt.fill_between(mean_values.index, mean_values - confidence, mean_values + confidence, alpha=0.2)

    # Add breaking points and other plot details
    for i, (frag_type, datasets) in enumerate(data.items()):
        breaking_point = mean(find_breakink_point_list(datasets[1]))
        plt.axvline(x=breaking_point, color=color_palette(i), ymax=0.1)

    plt.xlabel('Fragmentation step')
    plt.ylabel(ylabel)
    plt.title(f'{net}, ignoring isolated={ignore}')
    plt.legend()
    # plt.savefig(f"{filename} {net} ignore={ignore}.svg", format="svg")
    # plt.close()
    plt.show()

# Main script
fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div']
net = 'RGG'  # Example network type
ignore = False

data = load_data(fragmentation_types, net, ignore)

# Plotting
plot_data(data, 5, 'Pairwise Fst', 'Title here', net, ignore, 'fst')
plot_data(data, 3, 'Unscaled heterozygosity', 'Title here', net, ignore, 'het')


# pd.set_option('display.max_rows',None)

###############plot distributions
# plot distributions interval

def filter_intervals(df, interval_percentage=10):
    """
    Filter the DataFrame to include only specific intervals of steps.

    Args:
    df (pd.DataFrame): The original DataFrame with 'step' and 'replica' columns.
    interval_percentage (int): The percentage interval for filtering steps.

    Returns:
    pd.DataFrame: Filtered DataFrame.
    """
    # Determine the maximum step value
    max_step = df['step'].max()

    # Calculate interval step based on the percentage
    interval_step = max_step * interval_percentage // 100

    # Create a list of steps to include
    steps_to_include = list(range(0, max_step, interval_step))

    # Filter the DataFrame to include only these steps
    filtered_df = df[df['step'].isin(steps_to_include)]

    return filtered_df


with open('RGG, div_ignore_False.pickle', 'rb') as file:
    rand = pickle.load(file)


# First plot for 'fst'
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
# fig.suptitle('Fst - divisive fragmentation', fontsize=18)
# for ax in axes:
#     ax.tick_params(axis='both', which='major', labelsize=16)
#
# plt.show()
#
#
# # Second plot for 'het'
# df = filter_intervals(rand[2])

# plt.figure(figsize=(8, 8))
# fig, axes = joyplot(
#     data=df[['het', 'step']],
#     by='step', overlap=2,
#     colormap=plt.cm.viridis, fade=True, range_style='all',
#     linecolor="black", linewidth=0.1)
#
# fig.suptitle('Heterozygosity - divisive fragmentation', fontsize=18)
# for ax in axes:
#     ax.tick_params(axis='both', which='major', labelsize=16)
#     ax.set_xlim(-0.1, 1.3)
# plt.show()
#
#
#
#
#
#
#
#


#
# ############### plot networks along steps of fragmentation
############### to demonstrate the process

# with open('RGG, dist_ignore_False.pickle', 'rb') as file:
#     rand = pickle.load(file)
#
# with open('RGG, cor_ignore_False.pickle', 'rb') as file:
#     cor = pickle.load(file)
#
# with open('RGG, dist_ignore_False.pickle', 'rb') as file:
#     dist = pickle.load(file)
#
# with open('RGG, int_ignore_False.pickle', 'rb') as file:
#     int = pickle.load(file)
#
# with open('RGG, reg_ignore_False.pickle', 'rb') as file:
#     reg = pickle.load(file)
#
# with open('RGG, div_ignore_False.pickle', 'rb') as file:
#     div = pickle.load(file)
#
# print("finish load !!!")
#
#
# ##### plot fragmentation process
# net50=rand[1][10][50]
# net100=rand[1][10][100]
# net150=rand[1][10][150]
# net200=rand[1][10][200]
# net250=rand[1][10][250]
#
# pos = nx.spring_layout(net50, k=0.2, iterations=20,seed=50)
#
# fig, axes = plt.subplots(1, 5, figsize=(20, 4))
#
# # Draw each network
# nx.draw_networkx(net50, pos=pos, ax=axes[0], node_size=20, with_labels=False)
# axes[0].set_title("50",fontsize=22)
#
# nx.draw_networkx(net100, pos=pos, ax=axes[1], node_size=20, with_labels=False)
# axes[1].set_title("100",fontsize=22)
#
# nx.draw_networkx(net150, pos=pos, ax=axes[2], node_size=20, with_labels=False)
# axes[2].set_title("150",fontsize=22)
#
# nx.draw_networkx(net200, pos=pos, ax=axes[3], node_size=20, with_labels=False)
# axes[3].set_title("200",fontsize=22)
#
# nx.draw_networkx(net250, pos=pos, ax=axes[4], node_size=20, with_labels=False)
# axes[4].set_title("250",fontsize=22)
#
# # Adjust layout and display the figure
# plt.tight_layout()
# plt.show()