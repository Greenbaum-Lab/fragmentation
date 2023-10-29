import numpy as np
import networkx
import pickle
import pandas as pd

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from community import community_louvain

from network_analysis import plot_network_realization


def add_replica(df, step_column='step'):
    """
    add a column with the number of replica for each row based on 'step'
    """
    # Find where the step resets
    df.loc[:, 'replica'] = (df[step_column] < df[step_column].shift(1)).astype(int)

    # Compute the cumulative sum to get the replica number
    df.loc[:, 'replica'] = df['replica'].cumsum()

    return df

def extract_nodes(df):
    """
    tracking each node separately and extract its heterozygosity
    :param df: df of heterozygous of all nodes for all replicates (het dens)
    :return: df for each node along the fragmentation for each replica
    """
    # count how many nodes are in a network
    nodes = np.argmax(df['step'] != 0)


    # Randomly select a number between 0 and 49
    random_index = np.random.randint(0, nodes)

    # Select the rows corresponding to the randomly selected index for each step
    selected_rows = df.iloc[100::nodes].copy()

    return selected_rows



with open('RGG, dist_ignore_False.pickle', 'rb') as file:
    rand = pickle.load(file)


het = rand[2]
het.reset_index(drop=True, inplace=True)
het_values = extract_nodes(het)
het_values = add_replica(het_values)

# Create a pivot DataFrame
pivot_df = het_values.pivot(columns='replica', index='step', values='het')

# Plot using the pivot DataFrame
plt.figure(figsize=(10, 6))
plt.plot(pivot_df.index, pivot_df, color='grey', alpha=0.5)

plt.xlabel('Step')
plt.ylabel('Heterozygosity')
plt.title('Plot of Het against Step-distance')
plt.grid(True)
plt.show()




