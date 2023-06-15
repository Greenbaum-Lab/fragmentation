import random

import networkx as nx
from joypy import joyplot
from matplotlib import pyplot as plt
import pandas as pd
from funcs import make_networks, make_iterations, frag_random_giant_comp, frag_cor_giant_comp, frag_dist_giant_comp, \
    calculate_centrality, het_rand, het_cor, het_dist
from funcs import make_networks
from processes import find_breaking_point

n = 50  # no. of nodes
p = 0.5  # probability to connect nodes
seed = 98

random.seed(65)  # set random seed



# nets=make_networks(n_nets=3, n_nodes=n, connectivity=p, net_type='ER')
#
# rand = make_iterations(nets, fragmentation='rand')
# cor = make_iterations(nets, fragmentation='cor')
# dist = make_iterations(nets, fragmentation='dist')




# # Calculate the mean and median values over the 'step' column
# mean_rand = rand[0].groupby('step')['avg'].mean()
# mean_cor = cor[0].groupby('step')['avg'].mean()
# mean_dist = dist[0].groupby('step')['avg'].mean()
#
# # Calculate the confidence interval
# confidence_rand = rand[0].groupby('step')['avg'].std()
# confidence_cor = cor[0].groupby('step')['avg'].std()
# confidence_dist = dist[0].groupby('step')['avg'].std()
#
# # Plotting the line graph with mean and median values
# plt.plot(mean_rand, label='Rand')
# plt.plot(mean_cor, label='Cor')
# plt.plot(mean_dist, label='Dist')
#
# # Adding the confidence interval
# plt.fill_between(mean_rand.index, mean_rand - confidence_rand, mean_rand + confidence_rand, alpha=0.2)
# plt.fill_between(mean_cor.index, mean_cor - confidence_cor, mean_cor + confidence_cor, alpha=0.2)
# plt.fill_between(mean_dist.index, mean_dist - confidence_dist, mean_dist + confidence_dist, alpha=0.2)
#
# # Add labels and legend
# plt.xlabel('Fragmentation step')
# plt.ylabel('Pairwise Fst')
# plt.legend()
#
# # Display the plot
# plt.show()




















##############single net plots

# net = nx.erdos_renyi_graph(n=n, p=p)  # create network
# net = nx.erdos_renyi_graph(n=n,p=0.8)
net = nx.barabasi_albert_graph(n=50, m=2,seed=8)


rand = frag_random_giant_comp(net=net)
cor = frag_cor_giant_comp(net=net)
dist = frag_dist_giant_comp(net=net)

central_rand = calculate_centrality(rand[2])
central_cor = calculate_centrality(cor[2])
central_dist = calculate_centrality(dist[2])

brk_rand = find_breaking_point(rand[2])
brk_cor = find_breaking_point(cor[2])
brk_dist = find_breaking_point(dist[2])


# plotting avg and median

color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette

plt.plot(rand[0]['step'], rand[0]['avg'], label="avg rand", color=color_palette(0))
plt.plot(rand[0]['step'], rand[0]['median'], label="med rand", color=color_palette(0), linestyle='dashed')
plt.plot(cor[0]['step'], cor[0]['avg'], label="avg cor", color=color_palette(1))
plt.plot(cor[0]['step'], cor[0]['median'], label="med cor", color=color_palette(1), linestyle='dashed')
plt.plot(dist[0]['step'], dist[0]['avg'], label="avg dist", color=color_palette(2))
plt.plot(dist[0]['step'], dist[0]['median'], label="med dist", color=color_palette(2), linestyle='dashed')

plt.axvline(x=brk_rand, color=color_palette(0),ymax=0.1)
plt.axvline(x=brk_cor, color=color_palette(1),ymax=0.1)
plt.axvline(x=brk_dist, color=color_palette(2),ymax=0.1)
plt.xlabel('Fragmentation process')
plt.ylabel('Pairwise fst')
plt.legend()
plt.show()


# # # merge centrality measures and fst
merged_rand = pd.merge(central_rand, rand[0], on='step')
merged_cor = pd.merge(central_cor, cor[0], on='step')
merged_dist = pd.merge(central_dist, dist[0], on='step')


#find breaking point in clustering. if the axis is starting at 1 to 0 i take the complemetary (-1)
line1=max(central_rand['clustering'])-((brk_rand/len(rand[2]))*max(central_rand['clustering']))
line2=max(central_cor['clustering'])-((brk_cor/len(cor[2]))*max(central_cor['clustering']))
line3=max(central_dist['clustering'])-((brk_dist / len(dist[2])) * max(central_dist['clustering']))

plt.plot(merged_rand['clustering'], merged_rand['avg'], label="Random", color=color_palette(0))
plt.plot(merged_cor['clustering'], merged_cor['avg'], label="Correlated", color=color_palette(1))
plt.plot(merged_dist['clustering'], merged_dist['avg'], label="Distance", color=color_palette(2))

plt.axvline(x=line1, color=color_palette(0),ymax=0.1)
plt.axvline(x=line2, color=color_palette(1),ymax=0.1)
plt.axvline(x=line3, color=color_palette(2),ymax=0.1)
plt.xlim(1, 0)
plt.xlabel("Clustering")
plt.ylabel("Average Fst")
plt.legend()
plt.show()


# find breaking point in betweenness
line1=(brk_rand/len(rand[2]))*max(central_rand['betweenness'])
line2=(brk_cor/len(cor[2]))*max(central_cor['betweenness'])
line3=(brk_dist/len(dist[2]))*max(central_dist['betweenness'])

plt.plot(merged_rand['betweenness'], merged_rand['avg'], label="Random", color=color_palette(0))
plt.plot(merged_cor['betweenness'], merged_cor['avg'], label="Correlated", color=color_palette(1))
plt.plot(merged_dist['betweenness'], merged_dist['avg'], label="Distance", color=color_palette(2))

plt.axvline(x=line1, color=color_palette(0),ymax=0.1)
plt.axvline(x=line2, color=color_palette(1),ymax=0.1)
plt.axvline(x=line3, color=color_palette(2),ymax=0.1)
plt.xlabel("Betweenness")
plt.ylabel("Average Fst")
plt.legend()
plt.show()

desired_steps = list(range(0, 20, 2))
filtered_data = rand[1][rand[1]['step'].isin(desired_steps)]
# # plot distribution of Fst values - ridge-lines
plt.figure()
joyplot(
    data=filtered_data[['fst', 'step']],
    by='step', ylim=0, overlap=0.5,
    colormap=plt.cm.autumn, fade=True,
    figsize=(12, 8)
)

plt.title('pairwise Fst along random fragmentation', fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
plt.savefig('dens rand.svg', format='svg')
plt.show()


filtered_data = cor[1][cor[1]['step'].isin(desired_steps)]

plt.figure()
joyplot(
    data=filtered_data[['fst', 'step']],
    by='step', ylim=0, overlap=0.5,
    colormap=plt.cm.autumn, fade=True,
    figsize=(12, 8)
)
plt.title('pairwise Fst along correlated fragmentation', fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
plt.show()

filtered_data = dist[1][dist[1]['step'].isin(desired_steps)]

plt.figure()
joyplot(
    data=filtered_data[['fst', 'step']],
    by='step', ylim=0, overlap=0.5,
    colormap=plt.cm.autumn, fade=True,
    figsize=(12, 8)
)
plt.title('pairwise Fst along distance fragmentation', fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
plt.show()








#######################
####################### heterozygosity

rand_het = het_rand(net)
cor_het = het_cor(net)
dist_het = het_dist(net)


plt.plot(rand_het[0]['step'], rand_het[0]['avg'], label="avg rand", color=color_palette(0))
plt.plot(rand_het[0]['step'], rand_het[0]['median'], label="med rand", color=color_palette(0), linestyle='dashed')
plt.plot(cor_het[0]['step'], cor_het[0]['avg'], label="avg rand", color=color_palette(1))
plt.plot(cor_het[0]['step'], cor_het[0]['median'], label="med rand", color=color_palette(1), linestyle='dashed')
plt.plot(dist_het[0]['step'], dist_het[0]['avg'], label="avg dist", color=color_palette(2))
plt.plot(dist_het[0]['step'], dist_het[0]['median'], label="med dist", color=color_palette(2), linestyle='dashed')

plt.axvline(x=brk_rand, color=color_palette(0), ymax=0.1)
plt.axvline(x=brk_cor, color=color_palette(1), ymax=0.1)
plt.axvline(x=brk_dist, color=color_palette(2), ymax=0.1)
plt.xlabel('Fragmentation process')
plt.ylabel('Heterozygosity')
plt.legend()
plt.show()
