import numpy as np
import networkx
import pickle
import pandas as pd

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def add_replica(df, step_column='step'):
    # Find where the step resets
    df['replica'] = (df[step_column] < df[step_column].shift(1)).astype(int)

    # Compute the cumulative sum to get the replica number
    df['replica'] = df['replica'].cumsum()

    return df

nodes=50

def extract_nodes(df):

    # count how many nodes are in a network
    nodes = df['step'].value_counts().get(0, 0)
    print(nodes)
    # Randomly select a number between 0 and 49
    random_index = np.random.randint(0, nodes)

    # Select the rows corresponding to the randomly selected index for each step
    selected_rows = df.iloc[random_index::nodes]

    return selected_rows


with open('rand_include.pickle', 'rb') as file:
    rand = pickle.load(file)

het = rand[2]

het = add_replica(het)

het_values = extract_nodes(het)




# Creating a pivot DataFrame to organize the data by replica
pivot_df = het_values.pivot(columns='replica', values='het')

# Plotting all the lines in grey
plt.plot(het_values['step'].unique(), pivot_df, color='lightgrey',alpha=0.2)

plt.xlabel('Step')
plt.ylabel('Het')
plt.title('Plot of Het against Step')
plt.show()

# plt.plot(het_values['step'], het_values['het'])
# plt.xlabel('Step')
# plt.ylabel('Het')
# plt.title('Plot of Het against Step')
# plt.show()





