import networkx
import random
import statistics
from statistics import mean
import time
import seaborn as sns

import networkx as nx
from joypy import joyplot
from matplotlib import pyplot as plt
import pandas as pd

from processes import find_breaking_point, find_breakink_point_list, remove_edge_random, remove_edge_correlated, \
    remove_edge_distance
from funcs2 import frag_rand, frag_cor, frag_dist, het_rand, het_cor, het_dist, make_networks

n = 50  # no. of nodes
p = 0.4  # probability to connect nodes
n_rep = 100

# Record the starting time
start_time = time.time()

# create list off nets
nets = make_networks(n_nets=n_rep, n_nodes=n, connectivity=p, net_type='ER')
print("finish nets")
from multiprocessing import Pool


def parallelize_list_comprehension(nets, function):
    with Pool() as pool:
        return pool.map(function, nets)


breaking_point_rand = parallelize_list_comprehension(nets, remove_edge_random)
breaking_point_rand = find_breakink_point_list(breaking_point_rand)

breaking_point_cor = parallelize_list_comprehension(nets, remove_edge_correlated)
breaking_point_cor = find_breakink_point_list(breaking_point_cor)

breaking_point_dist = parallelize_list_comprehension(nets, remove_edge_distance)
breaking_point_dist = find_breakink_point_list(breaking_point_dist)


running_time = time.time() - start_time
print("Running time:", running_time, "seconds")
bins = 100
sns.histplot(data=breaking_point_rand, bins=bins, kde=True, color='blue', label='rand')
sns.histplot(data=breaking_point_cor, bins=bins, kde=True, color='red', label='cor')
sns.histplot(data=breaking_point_dist, bins=bins, kde=True, color='green', label='dist')

plt.xlabel('Value')
plt.ylabel('Count')
plt.title('Histogram')
plt.legend()
plt.savefig("breaking.png", format="png")

plt.show()

#######add parameters of percolation
def calculate_centrality(net: list) -> pd.DataFrame:
    """
    Calculate the degree of network
    :param net: list of migration networks
    :return: dataframe of degree and clustering for each step
    """
    m = net.copy()
    clustering = list(map(lambda x: nx.average_clustering(x), m))
    betweenness = list(map(lambda x: sum(nx.betweenness_centrality(x).values()) / len(x), m))
    step = range(len(m))
    d = {'step': step, 'clustering': clustering, 'betweenness': betweenness}
    df = pd.DataFrame(data=d)
    print(df)
    return df

