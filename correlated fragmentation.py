import random
import networkx as nx
import matplotlib.pyplot as plt


n = 5  # no. of nodes
p = 0.8  # probability to connect nodes
net = nx.erdos_renyi_graph(n, p)  # create ER network
n_frag = 20  # no. of fragmentation steps


def remove_edge_correlated(m, iterations: int, prob_increase: float):
    """
    remove edges from a migration network in a correlated sequence where the probability
    of removing an edge from a specific node increases if an edge was removed from the
    node in the previous iteration.
    :param m: migration network
    :param iterations: number of edges to remove
    :param prob_increase: the probability to get a correlated edge removal compared to random
    :return: migration network after edge removal
    """
    # create a list of all edges in the network
    edges = list(m.edges())

    # iterate for the specified number of iterations
    for i in range(iterations):
        # choose a random edge from the remaining edges
        edge = random.choice(edges)

        # remove the chosen edge from the network
        m.remove_edge(*edge)

        # update the list of remaining edges
        edges.remove(edge)

        # increase the probability of removing another edge from the same node
        if random.random() < prob_increase:
            node = random.choice(list(m.nodes()))
            edges = [e for e in edges if e[0] != node and e[1] != node]

    # return the modified network
    return m



import random
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas
import pandas as pd
from Transformation import m_to_f
from Transformation import m_to_t
from statistics import mean
from statistics import median


n = 10  # no. of nodes
p = 0.8  # probability to connect nodes
net = nx.erdos_renyi_graph(n, p, seed=55)  # create ER network
n_frag = 5  # no. of fragmentation steps


def calculate_fst(m) -> list:
    """
    Calculate Fst of M matrix after each edge removal
    :param m: initial migration network M of networkx
    :return: list of lists of Fst matrices for each step of fragmentation
    """
    fst_list = []
    for i in range(n_frag):
        M = nx.attr_matrix(m)[0]  # take the matrix of the net
        F = m_to_f(M)  # migration to fst function
        remove_edge_correlated(m, 1)  # use the remove edge function
        fst_list.append(F)  # add another item (fragmentation step) to the list

    return fst_list


array = calculate_fst(net)


def make_fst_dens(f: list) -> pd.DataFrame:
    """
    take a list of F metrics and return a dataframe without diagonal values (zero)
    :param f: list of fst metrics
    :return: dataframe with column represent each matrix
    """
    fst_dens = []
    for i in range(len(f)):
        F_no_diag = f[i][~np.eye(len(f[i]), dtype=bool)]  # remove diagonals of zero and concatenate array
        F_no_diag = np.ndarray.tolist(F_no_diag)  # transform to list
        fst_dens.append(F_no_diag)  # add another item (fragmentation step) to the list
    df = pd.DataFrame(fst_dens)
    df = df.transpose()
    df = df.stack().rename_axis(('delete', 'step')).reset_index(name='fst')
    df = df.drop(columns=['delete'])
    return df


fst_dens = make_fst_dens(array)
print(fst_dens)
# make dataframe with one column of all values-need to remove "a"


print(fst_dens)

def calculate_fst_data(f: pd.DataFrame) -> pd.DataFrame:
    """
     calculate the mean and median fst of each step
    :param f: dataframe of fst distribution in each step
    :return: dataframe of average and median for each step
    """
    avg = []
    med = []
    for i in range(n_frag):
        fst_avg = f[f['step'] == i]['fst']
        avg.append(mean(fst_avg))
        fst_med = f[f['step'] == i]['fst']
        med.append(median(fst_med))
        step = range(n_frag)
    d = {'step': step, 'avg': avg, 'median': med}
    df = pd.DataFrame(data=d)
    return df


fst_data = calculate_fst_data(fst_dens)

def calculate_het(m) -> pd.DataFrame:
    """
    calculate the unscaled heterozygosity from a diagonal of coalescence matrix
    :return: list of heterozygosity values for each population
    """
    h = []
    het = []
    df = pd.DataFrame()
    for i in range(n_frag):
        h = m_to_t(nx.attr_matrix(m)[0])
        h = np.diag(h)
        h = np.ndarray.tolist(h)
        het.append(h)
        remove_edge_correlated(m, 1)
    df = pd.DataFrame(het)
    # df = df.stack()
    # df = pd.DataFrame(df)
    df = df.stack().rename_axis(('step', 'delete')).reset_index(name='het')
    df=df.drop(columns=['delete'])
    return df

x = calculate_het(net)

print(x)



################ plotting avg and median
# plt.plot(fst_data['step'], fst_data['avg'], label="average")
# plt.plot(fst_data['step'], fst_data['median'], label="median")
#
# plt.xlabel('fragmentation process')
# plt.ylabel('Fst')
# plt.legend()
# plt.show()
#
# # plot distribution of Fst values - ridge-lines
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

