import random
import networkx as nx
import matplotlib.pyplot as plt

n = 5  # no. of nodes
p = 0.8  # probability to connect nodes
net = nx.erdos_renyi_graph(n, p)  # create ER network
n_frag = 20  # no. of fragmentation steps

# plotting avg and median
plt.plot(fst_data['step'], fst_data['avg'], label="average")
plt.plot(fst_data['step'], fst_data['median'], label="median")

plt.xlabel('fragmentation process')
plt.ylabel('Fst')
plt.legend()
plt.show()

# plot distribution of Fst values - ridge-lines
from joypy import joyplot

plt.figure()
joyplot(
    data=fst_dens[['fst', 'step']],
    by='step',
    colormap=plt.cm.autumn, fade=True,
    figsize=(12, 8)
)
plt.title('pairwise Fst along fragmentation', fontsize=20)
plt.show()

#
#
#
plt.figure()
joyplot(
    data=x[['het', 'step']],
    by='step',
    colormap=plt.cm.autumn, fade=True,
    figsize=(12, 8)
)
plt.title('pairwise Fst along fragmentation', fontsize=20)
plt.show()

#


# how to allow isolated populations in the analysis
# should  I put all in a class?
# func is slow in big networks
