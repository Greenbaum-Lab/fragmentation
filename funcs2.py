import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from Transformation import m_to_f
import random
net = nx.random_geometric_graph(8, 0.2, seed=12)
random.seed(56)

# function to extract all the connected components from a network
def get_components(net):

    # Get the connected components
    components = list(nx.connected_components(net))
    # Create subgraphs for each component
    subgraphs = []

    for component in components:
        if len(component) > 1:
            subgraph = net.subgraph(component)
            subgraphs.append(subgraph)
    return subgraphs


# function to combine arrays
def combine_arrays(arrays):

    # Calculate the total number of rows and columns
    total_rows = sum(array.shape[0] for array in arrays)
    total_cols = sum(array.shape[1] for array in arrays)

    # Initialize the combined array with zeros
    combined_array = np.zeros((total_rows, total_cols))

    # Copy the arrays into the combined array
    row_start = 0
    col_start = 0
    for array in arrays:
        row_end = row_start + array.shape[0]
        col_end = col_start + array.shape[1]
        combined_array[row_start:row_end, col_start:col_end] = array
        row_start = row_end
        col_start = col_end

    return combined_array


test1 = get_components(net)

for i in range(len(test1)):
    test1[i]=nx.attr_matrix(test1[i])[0]

combo=combine_arrays(test1)
print(combo)

#plot network
nx.draw_networkx(nx.Graph(combo))
plt.show()

from funcs import normalize
combo = normalize(combo)
print(combo)
#it works!
print(m_to_f(combo))


# #when i do the same thing with similar network it doesn't work
# net2 = nx.random_geometric_graph(8, 0.2, seed=123)
# test2 = get_components(net2)
#
# for i in range(len(test2)):
#     test2[i]=nx.attr_matrix(test2[i])[0]
#
# combo2=combine_arrays(test2)
# print(combo2)
#
# #plot network
# nx.draw_networkx(nx.Graph(combo2))
# plt.show()
#
# print(m_to_f(combo2))
