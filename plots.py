import random
import statistics
from statistics import mean
import time
import pickle

import networkx as nx
from joypy import joyplot
from matplotlib import pyplot as plt
import pandas as pd

from processes import find_breaking_point, find_breakink_point_list

# from funcs3 import make_replicates, make_replicates_new, make_networks

#
# print("here i start!")
n = 50  # no. of nodes
n_rep = 100
net = "RGG"
ignore = False
#
color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette
#
# # Record the starting time
# start_time = time.time()
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
#
# # Path and filename for the saved file using tuple
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
# #
# with open('ER, rand_ignore_False.pickle', 'rb') as file:
#     rand = pickle.load(file)
#
# with open('ER, cor_ignore_False.pickle', 'rb') as file:
#     cor = pickle.load(file)
#
# with open('ER, int_ignore_False.pickle', 'rb') as file:
#     int = pickle.load(file)
#
#
# breaking_point_rand = mean(find_breakink_point_list(rand[1]))
# breaking_point_cor = mean(find_breakink_point_list(cor[1]))
# breaking_point_int = mean(find_breakink_point_list(int[1]))
#
# #########plot fst
# # Calculate the mean and median values over the 'step' column
# mean_rand = rand[5].groupby('step')['avg'].mean()
# mean_cor = cor[5].groupby('step')['avg'].mean()
# mean_int = int[5].groupby('step')['avg'].mean()
#
#
# # Calculate the confidence interval
# confidence_rand = rand[5].groupby('step')['avg'].std()
# confidence_cor = cor[5].groupby('step')['avg'].std()
# confidence_int = int[5].groupby('step')['avg'].std()
#
# plt.figure()
# # Plotting the line graph with mean and median values
# plt.plot(mean_rand, label='Random')
# plt.plot(mean_cor, label='Correlated')
# plt.plot(mean_int, label='Intrusive')
#
# # Adding the confidence interval
# plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
# plt.fill_between(mean_cor.index, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
# plt.fill_between(mean_int.index, mean_int - confidence_int, mean_int + confidence_int, alpha=0.2)
#
# # add breaking points
# plt.axvline(x=breaking_point_rand, color=color_palette(0), ymax=0.1)
# plt.axvline(x=breaking_point_cor, color=color_palette(1), ymax=0.1)
# plt.axvline(x=breaking_point_int, color=color_palette(2), ymax=0.1)
#
# # Add labels and legend
# plt.xlabel('Fragmentation step')
# plt.ylabel('Pairwise Fst')
# plt.title(f'{net}, ignoring isolated={ignore} ')
# plt.legend()
#
# # Display the plot
# plt.savefig(f"fst {net} ignore={ignore}.svg", format="svg")
# plt.close()
# ####heterozygosity
#
# # Calculate the mean and median values over the 'step' column
# mean_rand = rand[3].groupby('step')['avg'].mean()
# mean_cor = cor[3].groupby('step')['avg'].mean()
# mean_int = int[3].groupby('step')['avg'].mean()
#
# # Calculate the confidence interval
# confidence_rand = rand[3].groupby('step')['avg'].std()
# confidence_cor = cor[3].groupby('step')['avg'].std()
# confidence_int = int[3].groupby('step')['avg'].std()
#
# plt.figure()
# # Plotting the line graph with mean and median values
# plt.plot(mean_rand, label='Random')
# plt.plot(mean_cor, label='Correlated')
# plt.plot(mean_int, label='Intrusive')
#
# # Adding the confidence interval
# plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
# plt.fill_between(mean_cor.index, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
# plt.fill_between(mean_int.index, mean_int - confidence_int, mean_int + confidence_int, alpha=0.2)
#
# # add breaking points
# plt.axvline(x=breaking_point_rand, color=color_palette(0), ymax=0.1)
# plt.axvline(x=breaking_point_cor, color=color_palette(1), ymax=0.1)
# plt.axvline(x=breaking_point_int, color=color_palette(2), ymax=0.1)
#
# # Add labels and legend
# plt.xlabel('Fragmentation step')
# plt.ylabel('Unscaled heterozygosity')
# plt.title(f'{net}, ignoring isolated={ignore}')
# plt.legend()
#
# # Display the plot
# plt.savefig(f"het {net} ignore={ignore}.svg", format="svg")
# plt.close()



#
# color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette

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
### save files as tuple
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


# with open('RGG, rand_ignore_False.pickle', 'rb') as file:
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
# breaking_point_rand = mean(find_breakink_point_list(rand[1]))
# breaking_point_cor = mean(find_breakink_point_list(cor[1]))
# breaking_point_int = mean(find_breakink_point_list(int[1]))
# breaking_point_reg = mean(find_breakink_point_list(reg[1]))
# breaking_point_div = mean(find_breakink_point_list(div[1]))
# breaking_point_dist = mean(find_breakink_point_list(dist[1]))
#
# #########plot fst
# # Calculate the mean and median values over the 'step' column
# mean_rand = rand[5].groupby('step')['avg'].mean()
# mean_cor = cor[5].groupby('step')['avg'].mean()
# mean_int = int[5].groupby('step')['avg'].mean()
# mean_reg = reg[5].groupby('step')['avg'].mean()
# mean_div = div[5].groupby('step')['avg'].mean()
# mean_dist = dist[5].groupby('step')['avg'].mean()
#
#
# # Calculate the confidence interval
# confidence_rand = rand[5].groupby('step')['avg'].std()
# confidence_cor = cor[5].groupby('step')['avg'].std()
# confidence_int = int[5].groupby('step')['avg'].std()
# confidence_reg = reg[5].groupby('step')['avg'].std()
# confidence_div = div[5].groupby('step')['avg'].std()
# confidence_dist = dist[5].groupby('step')['avg'].std()
#
# plt.figure()
# # Plotting the line graph with mean and median values
# plt.plot(mean_rand, label='Random')
# plt.plot(mean_cor, label='Correlated')
# plt.plot(mean_int, label='Intrusive')
# plt.plot(mean_reg, label='Regressive')
# plt.plot(mean_div, label='Divisive')
# plt.plot(mean_dist, label='Distance')
#
# # Adding the confidence interval
# plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
# plt.fill_between(mean_cor.index, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
# plt.fill_between(mean_int.index, mean_int - confidence_int, mean_int + confidence_int, alpha=0.2)
# plt.fill_between(mean_reg.index, mean_reg - confidence_reg, mean_reg + confidence_reg, alpha=0.2)
# plt.fill_between(mean_div.index, mean_div - confidence_div, mean_div + confidence_div, alpha=0.2)
# plt.fill_between(mean_dist.index, mean_dist - confidence_dist, mean_dist + confidence_dist, alpha=0.2)
#
# # add breaking points
# plt.axvline(x=breaking_point_rand, color=color_palette(0), ymax=0.1)
# plt.axvline(x=breaking_point_cor, color=color_palette(1), ymax=0.1)
# plt.axvline(x=breaking_point_int, color=color_palette(2), ymax=0.1)
# plt.axvline(x=breaking_point_reg, color=color_palette(3), ymax=0.1)
# plt.axvline(x=breaking_point_div, color=color_palette(4), ymax=0.1)
# plt.axvline(x=breaking_point_dist, color=color_palette(5), ymax=0.1)
#
# # Add labels and legend
# plt.xlabel('Fragmentation step')
# plt.ylabel('Pairwise Fst')
# plt.title(f'{net}, ignoring isolated={ignore} ')
# plt.legend()
#
# # Display the plot
# plt.savefig(f"fst {net} ignore={ignore}.svg", format="svg")
# plt.close()
#
# ####heterozygosity
#
# # Calculate the mean and median values over the 'step' column
# mean_rand = rand[3].groupby('step')['avg'].mean()
# mean_cor = cor[3].groupby('step')['avg'].mean()
# mean_int = int[3].groupby('step')['avg'].mean()
# mean_reg = reg[3].groupby('step')['avg'].mean()
# mean_div = div[3].groupby('step')['avg'].mean()
# mean_dist = dist[3].groupby('step')['avg'].mean()
#
# # Calculate the confidence interval
# confidence_rand = rand[3].groupby('step')['avg'].std()
# confidence_cor = cor[3].groupby('step')['avg'].std()
# confidence_int = int[3].groupby('step')['avg'].std()
# confidence_reg = reg[3].groupby('step')['avg'].std()
# confidence_div = div[3].groupby('step')['avg'].std()
# confidence_dist = dist[3].groupby('step')['avg'].std()
#
# plt.figure()
# # Plotting the line graph with mean and median values
# plt.plot(mean_rand, label='Random')
# plt.plot(mean_cor, label='Correlated')
# plt.plot(mean_int, label='Intrusive')
# plt.plot(mean_reg, label='Regressive')
# plt.plot(mean_div, label='Divisive')
# plt.plot(mean_dist, label='Distance')
#
# # Adding the confidence interval
# plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
# plt.fill_between(mean_cor.index, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
# plt.fill_between(mean_int.index, mean_int - confidence_int, mean_int + confidence_int, alpha=0.2)
# plt.fill_between(mean_reg.index, mean_reg - confidence_reg, mean_reg + confidence_reg, alpha=0.2)
# plt.fill_between(mean_div.index, mean_div - confidence_div, mean_div + confidence_div, alpha=0.2)
# plt.fill_between(mean_dist.index, mean_dist - confidence_dist, mean_dist + confidence_dist, alpha=0.2)
#
# # add breaking points
# plt.axvline(x=breaking_point_rand, color=color_palette(0), ymax=0.1)
# plt.axvline(x=breaking_point_cor, color=color_palette(1), ymax=0.1)
# plt.axvline(x=breaking_point_int, color=color_palette(2), ymax=0.1)
# plt.axvline(x=breaking_point_reg, color=color_palette(3), ymax=0.1)
# plt.axvline(x=breaking_point_div, color=color_palette(4), ymax=0.1)
# plt.axvline(x=breaking_point_dist, color=color_palette(5), ymax=0.1)
#
# # Add labels and legend
# plt.xlabel('Fragmentation step')
# plt.ylabel('Unscaled heterozygosity')
# plt.title(f'{net}, ignoring isolated={ignore}')
# plt.legend()
#
# # Display the plot
# plt.savefig(f"het {net} ignore={ignore}.svg", format="svg")
# plt.close()
#
#
#
#
#
#



###############plot distributions

# distribution fst
# plt.figure()
# joyplot(
#     data=rand[4][['fst', 'step']],
#     by='step', ylim=0, overlap=1,
#     colormap=plt.cm.viridis_r, fade=True,range_style='all',
#     figsize=(12, 8),  linecolor = "white", linewidth=0.1
# )
# plt.title('pairwise Fst in intance fragmentation',
#           fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
# plt.savefig("fst_distribution_ignore.png", format="png")
# plt.show()
# plt.figure()
#
#
# #distribution heterozygosity
# plt.figure()
# joyplot(
#     data=rand[2][['het', 'step']],
#     by='step', ylim=0, overlap=1,
#     colormap=plt.cm.viridis_r, fade=True,range_style='all',
#     figsize=(12, 8),  linecolor = "white", linewidth=0.1
# )
# plt.title('Heterozygosity in intance fragmentation',
#           fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
# plt.savefig("het_distribution_ignore.png", format="png")
# plt.show()
# plt.figure()

# with open('RGG, div_ignore_False.pickle', 'rb') as file:
#     rand = pickle.load(file)
#
# # plot distributions interval
# def plot_snapshot(df):
#
#     # Find the total number of rows
#     num_rows = df.shape[0]
#     # Calculate the 5% interval
#     interval = round(num_rows * 0.1)
#     # Create an auxiliary 'group' column
#     df = df.sort_values('step').reset_index(drop=True)  # ensure 'step' is sorted
#     df['group'] = df.index // interval
#
#     return df
#
# df = plot_snapshot(rand[4])
#
# #distribution heterozygosity
# plt.figure()
# joyplot(
#     data=df[['fst', 'group']],
#     by='group', ylim=0, overlap=1,
#     colormap=plt.cm.viridis_r, fade=True,range_style='all',
#     figsize=(12, 8),  linecolor = "white", linewidth=0.1
# )
# plt.title('fst in divisive fragmentation',
#           fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
# # plt.savefig("fst_distribution_include_int.svg", format="svg")
#
# plt.show()
# plt.figure()
#
# df = plot_snapshot(rand[2])
#
# #distribution heterozygosity
# plt.figure()
# joyplot(
#     data=df[['het', 'group']],
#     by='group', ylim=0, overlap=1,
#     colormap=plt.cm.viridis_r, fade=True,range_style='all',
#     figsize=(12, 8),  linecolor = "white", linewidth=0.1
# )
# plt.title('Heterozygosity in divisive fragmentation',
#           fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
# # plt.savefig("het_distribution_include_int.svg", format="svg")
# plt.show()
# plt.figure()






############### plot networks along steps of fragmentation
############### to demonstrate the process

with open('RGG, dist_ignore_False.pickle', 'rb') as file:
    rand = pickle.load(file)
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

net50=rand[1][10][50]
net100=rand[1][10][100]
net150=rand[1][10][150]
net200=rand[1][10][200]
net250=rand[1][10][250]

pos = nx.spring_layout(net50, k=0.2, iterations=20,seed=50)

fig, axes = plt.subplots(1, 5, figsize=(20, 4))

# Draw each network
nx.draw_networkx(net50, pos=pos, ax=axes[0], node_size=20, with_labels=False)
axes[0].set_title("50",fontsize=22)

nx.draw_networkx(net100, pos=pos, ax=axes[1], node_size=20, with_labels=False)
axes[1].set_title("100",fontsize=22)

nx.draw_networkx(net150, pos=pos, ax=axes[2], node_size=20, with_labels=False)
axes[2].set_title("150",fontsize=22)

nx.draw_networkx(net200, pos=pos, ax=axes[3], node_size=20, with_labels=False)
axes[3].set_title("200",fontsize=22)

nx.draw_networkx(net250, pos=pos, ax=axes[4], node_size=20, with_labels=False)
axes[4].set_title("250",fontsize=22)

# Adjust layout and display the figure
plt.tight_layout()
plt.show()