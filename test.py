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
from funcs2 import frag_rand, frag_cor, frag_dist, het_rand, het_cor, het_dist, make_networks, make_iterations_fst, \
    calculate_centrality, make_iterations_het, make_iterations_new, make_iterations_het_new

from funcs3 import make_replicates, make_replicates_new

# # # Load the tuple using pickle
# with open(pickle_filename, 'rb') as file:
#     loaded_tuple = pickle.load(file)


n = 50  # no. of nodes
p = 0.3  # probability to connect nodes
n_rep = 100
# seed = 98

color_palette = plt.get_cmap('tab10')  # You can change 'tab10' to any other available palette

# Record the starting time
start_time = time.time()

# # create list off nets
nets = make_networks(n_nets=n_rep, n_nodes=n, connectivity=p, net_type='RGG')

# run the pipeline for all fragmentation types
rand = make_replicates_new(nets=nets, frag_type='rand', ignore=True)
print("1")
cor = make_replicates_new(nets=nets, frag_type='cor', ignore=True)
print("2")
dist = make_replicates_new(nets=nets, frag_type='dist', ignore=True)
print("3")

#
# # Path and filename for the saved file using tuple
# pickle_filename = 'rand_ignore.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(rand, file)
#
# pickle_filename = 'cor_ignore.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(cor, file)
#
# pickle_filename = 'dist_ignore.pickle'
# with open(pickle_filename, 'wb') as file:
#     pickle.dump(dist, file)
#
#
#
