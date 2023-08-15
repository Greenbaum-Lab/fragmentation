import random
import statistics
from statistics import mean
import time
import pickle


from matplotlib import pyplot as plt
import pandas as pd

from processes import find_breaking_point, find_breakink_point_list

from funcs3 import make_replicates, make_replicates_new, make_networks


print("here i start!")
n = 50  # no. of nodes
n_rep = 100
net = "AB"
ignore = False

color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette

# Record the starting time
start_time = time.time()

# # create list off nets
nets = make_networks(n_nets=n_rep, n_nodes=n, net_type=net)

# run the pipeline for all fragmentation types
rand = make_replicates_new(nets=nets, frag_type='rand', ignore=ignore)
print("1")
cor = make_replicates_new(nets=nets, frag_type='cor', ignore=ignore)
print("2")
int = make_replicates_new(nets=nets, frag_type='int', ignore=ignore)
print("3")

# Path and filename for the saved file using tuple
pickle_filename = f'{net}, rand_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(rand, file)

pickle_filename = f'{net}, cor_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(cor, file)

pickle_filename = f'{net}, int_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(int, file)

# 
# # # Load the tuple using pickle
# with open('rand_include.pickle', 'rb') as file:
#     rand = pickle.load(file)
# 
# # # Load the tuple using pickle
# with open('cor_include.pickle', 'rb') as file:
#     cor = pickle.load(file)
# 
# # # Load the tuple using pickle
# with open('int_include.pickle', 'rb') as file:
#     int = pickle.load(file)



breaking_point_rand = mean(find_breakink_point_list(rand[1]))
breaking_point_cor = mean(find_breakink_point_list(cor[1]))
breaking_point_int = mean(find_breakink_point_list(int[1]))

#########plot fst
# Calculate the mean and median values over the 'step' column
mean_rand = rand[5].groupby('step')['avg'].mean()
mean_cor = cor[5].groupby('step')['avg'].mean()
mean_int = int[5].groupby('step')['avg'].mean()


# Calculate the confidence interval
confidence_rand = rand[5].groupby('step')['avg'].std()
confidence_cor = cor[5].groupby('step')['avg'].std()
confidence_int = int[5].groupby('step')['avg'].std()

plt.figure()
# Plotting the line graph with mean and median values
plt.plot(mean_rand, label='Random')
plt.plot(mean_cor, label='Correlated')
plt.plot(mean_int, label='Intrusive')

# Adding the confidence interval
plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
plt.fill_between(mean_cor.index, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
plt.fill_between(mean_int.index, mean_int - confidence_int, mean_int + confidence_int, alpha=0.2)

# add breaking points
plt.axvline(x=breaking_point_rand, color=color_palette(0), ymax=0.1)
plt.axvline(x=breaking_point_cor, color=color_palette(1), ymax=0.1)
plt.axvline(x=breaking_point_int, color=color_palette(2), ymax=0.1)

# Add labels and legend
plt.xlabel('Fragmentation step')
plt.ylabel('Pairwise Fst')
plt.title(f'{net}, ignoring isolated={ignore} ')
plt.legend()

# Display the plot
plt.savefig(f"fst {net} ignore={ignore}.png", format="png")
plt.close()
####heterozygosity

# Calculate the mean and median values over the 'step' column
mean_rand = rand[3].groupby('step')['avg'].mean()
mean_cor = cor[3].groupby('step')['avg'].mean()
mean_int = int[3].groupby('step')['avg'].mean()

# Calculate the confidence interval
confidence_rand = rand[3].groupby('step')['avg'].std()
confidence_cor = cor[3].groupby('step')['avg'].std()
confidence_int = int[3].groupby('step')['avg'].std()

plt.figure()
# Plotting the line graph with mean and median values
plt.plot(mean_rand, label='Random')
plt.plot(mean_cor, label='Correlated')
plt.plot(mean_int, label='Intrusive')

# Adding the confidence interval
plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
plt.fill_between(mean_cor.index, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
plt.fill_between(mean_int.index, mean_int - confidence_int, mean_int + confidence_int, alpha=0.2)

# add breaking points
plt.axvline(x=breaking_point_rand, color=color_palette(0), ymax=0.1)
plt.axvline(x=breaking_point_cor, color=color_palette(1), ymax=0.1)
plt.axvline(x=breaking_point_int, color=color_palette(2), ymax=0.1)

# Add labels and legend
plt.xlabel('Fragmentation step')
plt.ylabel('Unscaled heterozygosity')
plt.title(f'{net}, ignoring isolated={ignore}')
plt.legend()

# Display the plot
plt.savefig(f"het {net} ignore={ignore}.png", format="png")
plt.close()

running_time = time.time() - start_time
print("Running time:", running_time, "seconds")






print("here i start!")
n = 50  # no. of nodes
n_rep = 100
net = "SW"
ignore = False

color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette

# Record the starting time
start_time = time.time()

# # create list off nets
nets = make_networks(n_nets=n_rep, n_nodes=n, net_type=net)

# run the pipeline for all fragmentation types
rand = make_replicates_new(nets=nets, frag_type='rand', ignore=ignore)
print("1")
cor = make_replicates_new(nets=nets, frag_type='cor', ignore=ignore)
print("2")
int = make_replicates_new(nets=nets, frag_type='int', ignore=ignore)
print("3")

# Path and filename for the saved file using tuple
pickle_filename = f'{net}, rand_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(rand, file)

pickle_filename = f'{net}, cor_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(cor, file)

pickle_filename = f'{net}, int_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(int, file)

#
# # # Load the tuple using pickle
# with open('rand_include.pickle', 'rb') as file:
#     rand = pickle.load(file)
#
# # # Load the tuple using pickle
# with open('cor_include.pickle', 'rb') as file:
#     cor = pickle.load(file)
#
# # # Load the tuple using pickle
# with open('int_include.pickle', 'rb') as file:
#     int = pickle.load(file)



breaking_point_rand = mean(find_breakink_point_list(rand[1]))
breaking_point_cor = mean(find_breakink_point_list(cor[1]))
breaking_point_int = mean(find_breakink_point_list(int[1]))

#########plot fst
# Calculate the mean and median values over the 'step' column
mean_rand = rand[5].groupby('step')['avg'].mean()
mean_cor = cor[5].groupby('step')['avg'].mean()
mean_int = int[5].groupby('step')['avg'].mean()


# Calculate the confidence interval
confidence_rand = rand[5].groupby('step')['avg'].std()
confidence_cor = cor[5].groupby('step')['avg'].std()
confidence_int = int[5].groupby('step')['avg'].std()

plt.figure()
# Plotting the line graph with mean and median values
plt.plot(mean_rand, label='Random')
plt.plot(mean_cor, label='Correlated')
plt.plot(mean_int, label='Intrusive')

# Adding the confidence interval
plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
plt.fill_between(mean_cor.index, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
plt.fill_between(mean_int.index, mean_int - confidence_int, mean_int + confidence_int, alpha=0.2)

# add breaking points
plt.axvline(x=breaking_point_rand, color=color_palette(0), ymax=0.1)
plt.axvline(x=breaking_point_cor, color=color_palette(1), ymax=0.1)
plt.axvline(x=breaking_point_int, color=color_palette(2), ymax=0.1)

# Add labels and legend
plt.xlabel('Fragmentation step')
plt.ylabel('Pairwise Fst')
plt.title(f'{net}, ignoring isolated={ignore} ')
plt.legend()

# Display the plot
plt.savefig(f"fst {net} ignore={ignore}.png", format="png")
plt.close()
####heterozygosity

# Calculate the mean and median values over the 'step' column
mean_rand = rand[3].groupby('step')['avg'].mean()
mean_cor = cor[3].groupby('step')['avg'].mean()
mean_int = int[3].groupby('step')['avg'].mean()

# Calculate the confidence interval
confidence_rand = rand[3].groupby('step')['avg'].std()
confidence_cor = cor[3].groupby('step')['avg'].std()
confidence_int = int[3].groupby('step')['avg'].std()

plt.figure()
# Plotting the line graph with mean and median values
plt.plot(mean_rand, label='Random')
plt.plot(mean_cor, label='Correlated')
plt.plot(mean_int, label='Intrusive')

# Adding the confidence interval
plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
plt.fill_between(mean_cor.index, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
plt.fill_between(mean_int.index, mean_int - confidence_int, mean_int + confidence_int, alpha=0.2)

# add breaking points
plt.axvline(x=breaking_point_rand, color=color_palette(0), ymax=0.1)
plt.axvline(x=breaking_point_cor, color=color_palette(1), ymax=0.1)
plt.axvline(x=breaking_point_int, color=color_palette(2), ymax=0.1)

# Add labels and legend
plt.xlabel('Fragmentation step')
plt.ylabel('Unscaled heterozygosity')
plt.title(f'{net}, ignoring isolated={ignore}')
plt.legend()

# Display the plot
plt.savefig(f"het {net} ignore={ignore}.png", format="png")
plt.close()

running_time = time.time() - start_time
print("Running time:", running_time, "seconds")








# ###############plot distributions
#
# # distribution fst
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
#

## plot distributions interval
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
# df = plot_snapshot(int[4])
#
# #distribution heterozygosity
# plt.figure()
# joyplot(
#     data=df[['fst', 'group']],
#     by='group', ylim=0, overlap=1,
#     colormap=plt.cm.viridis_r, fade=True,range_style='all',
#     figsize=(12, 8),  linecolor = "white", linewidth=0.1
# )
# plt.title('fst in intance fragmentation',
#           fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
# plt.savefig("fst_distribution_include_int.svg", format="svg")
#
# plt.show()
# plt.figure()
#
# df = plot_snapshot(int[2])
# #distribution heterozygosity
# plt.figure()
# joyplot(
#     data=df[['het', 'group']],
#     by='group', ylim=0, overlap=1,
#     colormap=plt.cm.viridis_r, fade=True,range_style='all',
#     figsize=(12, 8),  linecolor = "white", linewidth=0.1
# )
# plt.title('Heterozygosity in intance fragmentation',
#           fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
# plt.savefig("het_distribution_include_int.svg", format="svg")
# plt.show()
# plt.figure()


#
#
#
# # run the pipeline for all fragmentation types
# rand = make_replicates_new(nets=nets, frag_type='rand', ignore=False)
# print("1")
# cor = make_replicates_new(nets=nets, frag_type='cor', ignore=False)
# print("2")
# int = make_replicates_new(nets=nets, frag_type='int', ignore=False)
# print("3")
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
# # Calculate the confidence interval
# confidence_rand = rand[5].groupby('step')['avg'].std()
# confidence_cor = cor[5].groupby('step')['avg'].std()
# confidence_int = int[5].groupby('step')['avg'].std()
#
# # Plotting the line graph with mean and median values
# plt.plot(mean_rand, label='Random')
# plt.plot(mean_cor, label='Correlated')
# plt.plot(mean_int, label='intance')
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
# plt.title(f'RGG({n},{p}), {n_rep} reps, include ones ')
# plt.legend()
#
# # Display the plot
# plt.savefig("fst_include.png", format="png")
# plt.show()
#
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
# # Plotting the line graph with mean and median values
# plt.plot(mean_rand, label='Random')
# plt.plot(mean_cor, label='Correlated')
# plt.plot(mean_int, label='intance')
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
# plt.title(f'RGG({n},{p}), {n_rep} reps, include ones ')
# plt.legend()
#
# # Display the plot
# plt.savefig("het_include.png", format="png")
# plt.show()
# running_time = time.time() - start_time
# print("Running time:", running_time, "seconds")
#
#
#
#
# # Path and filename for the saved file using tuple
# pickle_filename = 'rand_include.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(rand, file)
#
# pickle_filename = 'cor_include.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(cor, file)
#
# pickle_filename = 'int_include.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(int, file)
#
#
# ###############plot distributions
#
# # distribution fst
# plt.figure()
# joyplot(
#     data=rand[4][['fst', 'step']],
#     by='step', ylim=0, overlap=1,
#     colormap=plt.cm.viridis_r, fade=True,range_style='all',
#     figsize=(12, 8),  linecolor = "white", linewidth=0.1
# )
# plt.title('pairwise Fst in intance fragmentation',
#           fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
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
# plt.show()
# plt.figure()
#

### plot distributions interval
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
# #distribution heterozygosity
# plt.figure()
# joyplot(
#     data=df[['fst', 'group']],
#     by='group', ylim=0, overlap=1,
#     colormap=plt.cm.viridis_r, fade=True,range_style='all',
#     figsize=(12, 8),  linecolor = "white", linewidth=0.1
# )
# plt.title('Heterozygosity in intance fragmentation',
#           fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
# plt.show()
# plt.figure()
