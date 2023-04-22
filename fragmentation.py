import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import random


def remove_edge_random(m, n: int):
    """
    Remove a random edge from net m of type networkx
    :param m: initial migration net m
    :param n: no. of fragmentation steps
    :return: net after edge removal
    """

    for i in range(n):
        edges = list(nx.edges(m))
        edges_to_remove = (random.sample(edges, k=1))  # choose a random edge
        m.remove_edge(*(edges_to_remove[0]))

    return m


def remove_edge_correlated(m, n: int):
    """
    remove edges from a migration network in a correlated sequence where the probability
    of removing an edge from a specific node increases if an edge was removed from the
    node in the previous iteration.
    :param m: migration network
    :param n: number of edges to remove
    :return: migration network after edge removal
    """
    prob_increase = 0.8
    # create a list of all edges in the network
    edges = list(m.edges())

    # iterate for the specified number of iterations
    for i in range(n):
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
