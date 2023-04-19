import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from Transformation import m_to_f
from Transformation import m_to_t
from statistics import mean
from statistics import median
from fragmentation import remove_edge_random
from fragmentation import remove_edge_correlated

# create a network
n = 10  # no. of nodes
p = 0.8  # probability to connect nodes
net = nx.erdos_renyi_graph(n, p, seed=55)
n_frag = 5  # no. of fragmentation steps

# pos = nx.spring_layout(net, seed=55)
nx.draw(net, with_labels=True)
plt.show()


def calculate_fst(m: np.ndarray, frag_process) -> list:
    """
    Calculate Fst of M matrix after each random edge removal
    :param frag_process: type of fragmentation (random, correlated)
    :param m: initial migration network M of networkx
    :return: list of lists of Fst matrices for each step of fragmentation
    """
    fst_list = []
    for i in range(n_frag):
        M = nx.attr_matrix(m)[0]  # take the matrix of the net
        F = m_to_f(M)  # migration to fst function
        frag_process(m=m, n=1)  # use the remove edge function
        fst_list.append(F)  # add another item (fragmentation step) to the list

    return fst_list


def calculate_het(m: np.ndarray, frag_process) -> list:
    """
    Calculate heterozygosity based on coalescence matrix diagonal of M matrix after each random edge removal
    :param frag_process: type of fragmentation (random, correlated)
    :param m: initial migration network M of networkx
    :return: list of lists of coalescence matrices for each step of fragmentation
    """
    het_list = []
    h = []
    for i in range(n_frag):
        M = nx.attr_matrix(m)[0]  # take the matrix of the net
        T = m_to_t(M)  # migration to fst function
        h = np.diag(T)
        h = np.ndarray.tolist(h)
        frag_process(m=m, n=1)  # use the remove edge function
        het_list.append(h)  # add another item (fragmentation step) to the list

    return het_list



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

x=calculate_fst(net,remove_edge_random)
print(make_fst_dens(x).head())

def calculate_het(m: np.ndarray) -> pd.DataFrame:
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
        remove_edge_random(m, 1)
    df = pd.DataFrame(het)
    df = df.stack().rename_axis(('step', 'delete')).reset_index(name='het')
    df = df.drop(columns=['delete'])
    return df



# make dataframe with one column of all values-need to remove "a"




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

