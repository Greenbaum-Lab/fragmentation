import random
import networkx as nx
import matplotlib.pyplot as plt
from joypy import joyplot
from fragmentation import remove_edge_random
from fragmentation import remove_edge_correlated
from fragmentation import remove_edges_distance
from funcs import make_fst_list
from funcs import make_fst_dist
from funcs import make_fst_stat
from funcs import make_het_list
from funcs import make_het_dist
from funcs import make_het_stat


n = 20  # no. of nodes
p = 0.8  # probability to connect nodes
# seed = 666
net = nx.erdos_renyi_graph(n, p)  # create network
n_frag = 500  # no. of fragmentation steps
pos = nx.spring_layout(net, seed=55)  # set the fixed position for plotting the network
# random.seed(12)  # set random seed


migration_list = remove_edge_random(migration=net, n=n_frag)
fst = make_fst_list(migration_list=migration_list)
fst_dens = make_fst_dist(fst)
fst_stat = make_fst_stat(fst_dens)
print(fst_stat)

# # # plot distribution of Fst values - ridge-lines
plt.figure()
joyplot(
    data=fst_dens[['fst', 'step']],
    by='step',
    colormap=plt.cm.autumn, fade=True,
    figsize=(12, 8)
)
plt.title('pairwise Fst along fragmentation', fontsize=16)
plt.show()




# heterozygosity = calculate_het(net, remove_edge_random)
# heterozygosity_dens = make_het_dist(heterozygosity)
# het_stat = make_het_stat(heterozygosity_dens)


# plotting avg and median
# plt.plot(fst_stat['step'], fst_stat['avg'], label="average")
# plt.plot(fst_stat['step'], fst_stat['median'], label="median")
#
# plt.xlabel('fragmentation process')
# plt.ylabel('pairwise fst')
# plt.legend()
# plt.show()
#

# # plotting avg and median
# plt.plot(het_stat['step'], het_stat['avg'], label="average")
# plt.plot(het_stat['step'], het_stat['median'], label="median")
#
# plt.xlabel('fragmentation process')
# plt.ylabel('unscaled heterozygosity')
# plt.legend()
# plt.show()
#
