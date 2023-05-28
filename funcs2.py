import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Transformation import m_to_f
import random
net = nx.random_geometric_graph(10, 0.2,seed=337)
nx.draw_networkx(net, with_labels=False)
plt.show()

random.seed(56)



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


test=get_components(net)
print(test)


for i in range(len(test)):
    test[i]=nx.attr_matrix(test[i])[0]

print(test)


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

combo=combine_arrays(test)
print(combo)

nx.draw_networkx(nx.Graph(combo))
plt.show()

array1 = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
array2 = np.array([[0, 1], [1, 0]])

x=combine_arrays([array1,array2])
print(x)
print(m_to_f(x))

