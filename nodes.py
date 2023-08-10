import numpy as np
import networkx
import pickle
import pandas as pd

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def add_replica(df, step_column='step'):
    """
    add a column with the number of replica for each row based on 'step'
    """
    # Find where the step resets
    df['replica'] = (df[step_column] < df[step_column].shift(1)).astype(int)

    # Compute the cumulative sum to get the replica number
    df['replica'] = df['replica'].cumsum()

    return df


def extract_nodes(df):
    """
    tracking each node sepreatly and extract its heterozygosity
    :param df: df of heyterozygosity of all nodes for all replicates (het dens)
    :return: df for each node along the fragmentation for each replica
    """
    # count how many nodes are in a network
    nodes = np.argmax(df['step'] != 0)
    print(nodes)
    nodes=2450
    # Randomly select a number between 0 and 49
    random_index = np.random.randint(0, nodes)

    # Select the rows corresponding to the randomly selected index for each step
    selected_rows = df.iloc[100::nodes]

    return selected_rows


with open('dist_include.pickle', 'rb') as file:
    rand = pickle.load(file)

het = rand[4]
het = het.iloc[:650000]
het.reset_index(drop=True, inplace=True)
print(het)
het_values = extract_nodes(het)
het_values = add_replica(het_values)
print(het_values)

# plot all lines\populations together
# Creating a pivot DataFrame to organize the data by replica
pivot_df = het_values.pivot(columns='replica',index='step', values='fst')
print(pivot_df)
plt.plot(het_values['step'].unique(), pivot_df, color='lightgrey',alpha=0.5)

plt.xlabel('Step')
plt.ylabel('Heterozygosity')
plt.title('Plot of Het against Step')
plt.show()




