from funcs import make_networks, make_iterations
from funcs import make_networks

n = 15  # no. of nodes
p = 0.8  # probability to connect nodes
seed = 987

# pos = nx.spring_layout(net, seed=98)  # set the fixed position for plotting the network
# random.seed(65)  # set random seed



nets=make_networks(n_nets=3, n_nodes=n, connectivity=p, net_type='ER')

rand = make_iterations(nets, fragmentation='rand')
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
# net = nx.random_geometric_graph(n=n, radius=0.8)

# rand = frag_random_giant_comp(net=net)
# cor = frag_cor_giant_comp(net=net)
# dist = frag_dist_giant_comp(net=net)
#
# central_rand = calculate_centrality(rand[2])
# central_cor = calculate_centrality(cor[2])
# central_dist = calculate_centrality(dist[2])
#
# brk_rand = find_breaking_point(rand[2])
# brk_cor = find_breaking_point(cor[2])
# brk_dist = find_breaking_point(dist[2])
#
#
# # plotting avg and median
#
# color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette
#
# plt.plot(rand[0]['step'], rand[0]['avg'], label="avg rand", color=color_palette(0))
# plt.plot(rand[0]['step'], rand[0]['median'], label="med rand", color=color_palette(0), linestyle='dashed')
# plt.plot(cor[0]['step'], cor[0]['avg'], label="avg cor", color=color_palette(1))
# plt.plot(cor[0]['step'], cor[0]['median'], label="med cor", color=color_palette(1), linestyle='dashed')
# plt.plot(dist[0]['step'], dist[0]['avg'], label="avg dist", color=color_palette(2))
# plt.plot(dist[0]['step'], dist[0]['median'], label="med dist", color=color_palette(2), linestyle='dashed')
#
# plt.axvline(x=brk_rand, color=color_palette(0),ymax=0.1)
# plt.axvline(x=brk_cor, color=color_palette(1),ymax=0.1)
# plt.axvline(x=brk_dist, color=color_palette(2),ymax=0.1)
# plt.xlabel('Fragmentation process')
# plt.ylabel('Pairwise fst')
# plt.legend()
# plt.show()
#
# #
# #
# #
# # # # merge centrality measures and fst
# merged_rand = pd.merge(central_rand, rand[0], on='step')
# merged_cor = pd.merge(central_cor, cor[0], on='step')
# merged_dist = pd.merge(central_dist, dist[0], on='step')
#
#
# #find breaking point in clustering
# line1=(brk_rand/len(rand[2]))*max(central_rand['clustering'])
# line2=(brk_cor/len(cor[2]))*max(central_cor['clustering'])
# line3=(brk_dist/len(dist[2]))*max(central_dist['clustering'])
#
# plt.plot(merged_rand['clustering'], merged_rand['avg'], label="Random", color=color_palette(0))
# plt.plot(merged_cor['clustering'], merged_cor['avg'], label="Correlated", color=color_palette(1))
# plt.plot(merged_dist['clustering'], merged_dist['avg'], label="Distance", color=color_palette(2))
#
# plt.axvline(x=line1, color=color_palette(0),ymax=0.1)
# plt.axvline(x=line2, color=color_palette(1),ymax=0.1)
# plt.axvline(x=line3, color=color_palette(2),ymax=0.1)
# plt.xlim(1, 0)
# plt.xlabel("Clustering")
# plt.ylabel("Average Fst")
# plt.legend()
# plt.show()
#
#
#
# #find breaking point in betweenness
# line1=(brk_rand/len(rand[2]))*max(central_rand['betweenness'])
# line2=(brk_cor/len(cor[2]))*max(central_cor['betweenness'])
# line3=(brk_dist/len(dist[2]))*max(central_dist['betweenness'])
#
# plt.plot(merged_rand['betweenness'], merged_rand['avg'], label="Random", color=color_palette(0))
# plt.plot(merged_cor['betweenness'], merged_cor['avg'], label="Correlated", color=color_palette(1))
# plt.plot(merged_dist['betweenness'], merged_dist['avg'], label="Distance", color=color_palette(2))
#
# plt.axvline(x=line1, color=color_palette(0),ymax=0.1)
# plt.axvline(x=line2, color=color_palette(1),ymax=0.1)
# plt.axvline(x=line3, color=color_palette(2),ymax=0.1)
# plt.xlabel("Betweenness")
# plt.ylabel("Average Fst")
# plt.legend()
# plt.show()
#

# # # plot distribution of Fst values - ridge-lines
# plt.figure()
# joyplot(
#     data=rand[1][['fst', 'step']],
#     by='step', ylim=0, overlap=0.5,
#     colormap=plt.cm.autumn, fade=True,
#     figsize=(12, 8)
# )
#
# plt.title('pairwise Fst along random fragmentation', fontsize=16)
# plt.show()
#

#
# plt.figure()
# joyplot(
#     data=dist[1][['fst', 'step']],
#     by='step', ylim=0, overlap=0.5,
#     colormap=plt.cm.autumn, fade=True,
#     figsize=(12, 8)
# )
# plt.title('pairwise Fst along distance fragmentation', fontsize=16)
# plt.show()
#
#
# plt.figure()
# joyplot(
#     data=cor[1][['fst', 'step']],
#     by='step', ylim=0, overlap=0.5,
#     colormap=plt.cm.autumn, fade=True,
#     figsize=(12, 8)
# )
# plt.title('pairwise Fst along correlated fragmentation', fontsize=16)
# plt.show()
#
#
#
#
#
#
#
# #######################
# ####################### heterozygosity
#
# het_rand = het_rand(net)
# het_cor = het_cor(net)
# het_dist = het_dist(net)
#
#
# plt.plot(het_rand[0]['step'], het_rand[0]['avg'], label="avg rand", color=color_palette(0))
# plt.plot(het_rand[0]['step'], het_rand[0]['median'], label="med rand", color=color_palette(0), linestyle='dashed')
# plt.plot(het_cor[0]['step'], het_cor[0]['avg'], label="avg rand", color=color_palette(1))
# plt.plot(het_cor[0]['step'], het_cor[0]['median'], label="med rand", color=color_palette(1), linestyle='dashed')
# plt.plot(het_dist[0]['step'], het_dist[0]['avg'], label="avg dist", color=color_palette(2))
# plt.plot(het_dist[0]['step'], het_dist[0]['median'], label="med dist", color=color_palette(2), linestyle='dashed')
#
# plt.axvline(x=brk_rand, color=color_palette(0), ymax=0.1)
# plt.axvline(x=brk_cor, color=color_palette(1), ymax=0.1)
# plt.axvline(x=brk_dist, color=color_palette(2), ymax=0.1)
# plt.xlabel('Fragmentation process')
# plt.ylabel('Heterozygosity')
# plt.legend()
# plt.show()
