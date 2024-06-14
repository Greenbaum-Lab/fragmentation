from statistics import mean
import pickle
import pandas as pd
from joypy import joyplot
from matplotlib import pyplot as plt

from funcs_analysis import load_data, plot_data, filter_intervals, plot_all_distributions, plot_distribution

# from funcs_initial_data import make_networks, make_replicates_new

# pd.set_option('display.max_rows',None)


########## run pipeline
# print("here i start!")
# n = 50  # no. of nodes
# n_rep = 100  # no. of replicates
# net = "RGG"
# ignore = False
# #
# # # # create list off nets
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
# opt = make_replicates_new(nets=nets, frag_type='opt', ignore=ignore)
# print("7")
# opt2 = make_replicates_new(nets=nets, frag_type='opt2', ignore=ignore)
# print("8")
# wrst = make_replicates_new(nets=nets, frag_type='wrst', ignore=ignore)
# print("9")
#
#
#
# # save files as tuple
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
#
# pickle_filename = f'{net}, opt2_ignore_{ignore}.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(opt2, file)
#
# pickle_filename = f'{net}, wrst_ignore_{ignore}.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(wrst, file)
# ########## finish pipeline



fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']
# fragmentation_types = ['int']
net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)


# Plot fst and het along fragmentation
plot_data(data, 5, 'Pairwise Fst',measure='fst')
plot_data(data, 3, 'Heterozygosity',measure='heterozygosity')


##############plot distributions



def filter_intervals(df, interval_percentage=25):
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
    steps_to_include = steps_to_include[:4]

    # Filter the DataFrame to include only these steps
    filtered_df = df[df['step'].isin(steps_to_include)]
    return filtered_df




def plot_distribution(df, measure='het' or 'fst',type=str):
    """ Plot the distribution of a given measure across different steps.
    """
    # Create a figure and axes
    fig, ax = plt.subplots()
    # Get unique steps
    unique_steps = df['step'].unique()

    # Generate reversed color gradient
    colors = plt.cm.YlOrRd(np.linspace(0.3, 1, len(unique_steps)))[::-1]

    # Plot histogram for each step with increasing alpha
    for i, step in enumerate(unique_steps):
        if measure == 'fst':
            values = df[df['step'] == step]['fst']
        if measure == 'het':
            values = df[df['step'] == step]['het']
        ax.hist(values, bins=40, alpha=0.7, label=f'Step {step}', density=True,
                color=colors[i], edgecolor='black')

    # Set titles and labels
    ax.set_xlabel('Fst' if measure == 'fst' else 'Heterozygosity')
    ax.set_ylabel('Density (%) ')
    ax.legend()

    # Optional: set x and y limits
    # ax.set_xlim(0, 1.4)
    ax.set_ylim(0, 20)

    # Show the plot
    plt.savefig(f'./figs/dist_{measure}_{type}.jpg', format="jpg")
    plt.show()


fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt']

frag='div'
fragmentation_types = [frag]
net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)


df = filter_intervals(data[frag][2])

plot_distribution(df,measure='het',type=frag)