import copy

import networkx as nx
import random
import math

import numpy as np
from matplotlib import pyplot as plt


def remove_edge_random(net: nx.Graph, n: int) -> list:
    """
    Remove a random edge from migration network of type networkx
    :param net:  initial migration network
    :param n: no. of fragmentation steps
    :return: list of networks after n edge removal
    """
    migration_list = [net.copy()]
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    pos = nx.spring_layout(net.copy(), seed=55)

    for i in range(n):
        migration = migration_list[-1].copy()  # takes the last item in list

        #  plot initial network
        if i == 0:
            nx.draw_networkx(migration, pos=pos, ax=axs[i])
            axs[i].set_title(f"Original network")

        edges = list(nx.edges(migration))
        edges_to_remove = (random.sample(edges, k=1))  # choose a random edge
        migration.remove_edge(*(edges_to_remove[0]))

        # stop when network breaks
        if nx.is_connected(migration):
            migration_list.append(migration)

        else:

            #  plot last network before breaking
            nx.draw_networkx(migration_list[i], pos=pos, ax=axs[1])
            axs[1].set_title(f"{i} edges removed randomly")

            break

    plt.show()
    return migration_list


def remove_edge_correlated(net: nx.Graph, n: int) -> list:
    """
    Remove a correlated edge from migration network of type networkx
    :param net: initial migration network
    :param n: no. of fragmentation steps
    :return: list of networks after n edge removal
    """
    migration_list = [net.copy()]
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    pos = nx.spring_layout(net.copy(), seed=55)

    # takes the last item in list
    migration = migration_list[-1].copy()

    # choose a random edge from the remaining edges
    edge = random.choice(list(migration.edges))

    # remove the chosen edge from the network
    migration.remove_edge(*edge)

    # add the network to the list
    migration_list.append(migration)

    for i in range(n):

        #  plot initial network
        if i == 0:
            nx.draw_networkx(migration, pos=pos, ax=axs[i])
            axs[i].set_title(f"Original network")

        # takes the last item in list
        migration = migration_list[-1].copy()

        # takes the edges of the nodes of the removed edge
        edges_a = list(migration.edges(edge[1]))
        edges_b = list(migration.edges(edge[0]))

        # make a list of the edges of the two nodes
        edges = edges_a + edges_b

        # sample an edge from the edges
        edge = random.choice(edges)

        # remove the chosen edge from the network
        migration.remove_edge(*edge)

        if nx.is_connected(migration):
            migration_list.append(migration)
        else:

            #  plot last network before breaking
            nx.draw_networkx(migration_list[i], pos=pos, ax=axs[1])
            axs[1].set_title(f"{i} edges removed by correlation")
            break

    plt.show()
    return migration_list


def remove_edge_distance(net: nx.Graph, n: int) -> list:
    """
    Remove edge from migration network of type networkx
    :param net:  initial migration network
    :param n: no. of fragmentation steps
    :return: list of networks after n edge removal
    """

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    pos = nx.spring_layout(net.copy(), seed=55)
    migration_list = [net.copy()]
    # Get all the edges of the graph
    edges = net.copy().edges()
    # Calculate the euclidean distance between all nodes
    distances = {edge: round(
        ((pos[edge[1]][1] - pos[edge[0]][1]) ** 2 + (pos[edge[1]][0] - pos[edge[0]][0]) ** 2) ** 0.5, 2)
        for edge in edges}
    # Sort the edges by their distances in descending order
    edges = sorted(edges, key=distances.get, reverse=True)

    for i in range(n):
        # takes the last item in list
        migration = migration_list[-1].copy()

        if i == 0:
            #  plot initial network
            nx.draw_networkx(migration, pos=pos, ax=axs[i])
            axs[i].set_title(f"Original network")

        # Remove the longest edge from network
        migration.remove_edge(*edges[0])
        # Remove the longest edge from the edges list
        edges.pop(0)

        #  add new network to list if it is connected
        if nx.is_connected(migration):
            migration_list.append(migration)

        else:
            #  plot last network before breaking
            nx.draw_networkx(migration_list[i], pos=pos, ax=axs[1])
            axs[1].set_title(f"{i} edges removed by distance")
            break

    plt.show()
    return migration_list


# def remove_edge_random(net: nx.Graph, n: int) -> list:
#     """
#     Remove a random edge from migration network of type networkx
#     :param net:  initial migration network
#     :param n: no. of fragmentation steps
#     :return: list of networks after n edge removal
#     """
#     migration_list = [net.copy()]
#
#     for i in range(n):
#         migration = migration_list[-1].copy()  # takes the last item in list
#
#         edges = list(nx.edges(migration))
#         edges_to_remove = (random.sample(edges, k=1))  # choose a random edge
#         migration.remove_edge(*(edges_to_remove[0]))
#         edges.remove(edges_to_remove[0])
#         migration_list.append(migration)
#
#     return migration_list


def remove_edge_random_giant_comp(net: nx.Graph) -> list:
    """
    Remove a random edge from the input network in each iteration.
    If the resulting graph is disconnected, it only keeps the giant component.
    Keep track of all generated graphs in a list until only two nodes remain.

    :param net: initial networkx object
    :return: list of networkx objects
    """
    migration = net.copy()
    migration_list = []  # initialize list with the original network
    while len(migration) > 2:  # stop when the network includes only two connected nodes
        # choose a random edge and remove it
        edges = list(migration.edges())
        edge = random.choice(edges)
        migration.remove_edge(*edge)

        # check if the resulting graph is connected
        if not nx.is_connected(migration):
            # if not, keep only the giant component
            largest_cc = max(nx.connected_components(migration), key=len)
            migration = migration.subgraph(largest_cc).copy()

        # add the resulting graph to the list
        migration_list.append(migration.copy())

    return migration_list


def remove_edge_correlated_giant_comp(net: nx.Graph) -> list:
    """
    Remove a correlated edge from migration network of type networkx
    each iteration we choose an edge to remove from the list of edges connected to it
    if the node is no longer part of the network we choose a random node
    :param net: initial migration network
    :param n: no. of fragmentation steps
    :return: list of networks after n edge removal
    """

    migration = net.copy()
    migration_list = [net.copy()]  # initialize list with the original network

    # choose a random edge for start
    edge = random.choice(list(migration.edges))

    # take the initial nodes to remove edges from
    edge_a = edge[0]
    edge_b = edge[1]

    while len(migration) > 2:  # stop when the network includes only two connected nodes
        # print(nx.attr_matrix(migration))
        # takes the edges of the nodes of the removed edge
        edges_a = list(migration.edges(edge_a))
        edges_b = list(migration.edges(edge_b))

        # make a list of the edges of the two nodes
        edges = edges_a + edges_b

        # sample an edge from the edges
        edge = random.choice(list(edges))

        # remove the chosen edge from the network
        migration.remove_edge(*edge)

        # check if the resulting graph is connected
        if not nx.is_connected(migration):
            # if not, keep only the giant component
            largest_cc = max(nx.connected_components(migration), key=len)
            migration = migration.subgraph(largest_cc).copy()

        # choose the nodes to remove edges from
        # if node doesn't exist in new network choose a random node
        if edge[0] in migration.nodes:
            edge_a = edge[0]
        else:
            edge_a = random.choice(list(migration.nodes))

        if edge[1] in migration.nodes:
            edge_b = edge[1]
        else:
            edge_b = random.choice(list(migration.nodes))

        # add the resulting graph to the list
        migration_list.append(migration.copy())

    return migration_list


def remove_edge_distance_giant_comp(net: nx.Graph) -> list:
    """
    Remove edge based on distance from migration network of type networkx
    :param net:  initial migration network
    :return: list of networks after n edge removal
    """

    pos = nx.spring_layout(net.copy(), seed=55)

    migration = net.copy()
    migration_list = []  # initialize list of networks

    while len(migration) > 2:  # stop when the network includes only two connected nodes

        edges = migration.edges()

        # Calculate the euclidean distance between all nodes and create a dict
        distances = {edge: round(
            ((pos[edge[1]][1] - pos[edge[0]][1]) ** 2 + (pos[edge[1]][0] - pos[edge[0]][0]) ** 2) ** 0.5, 2)
            for edge in edges}

        edges = sorted(edges, key=distances.get, reverse=True)

        # Remove the longest edge from network
        migration.remove_edge(*edges[0])

        # Remove the longest edge from the edges list
        edges.pop(0)

        # check if the resulting graph is connected
        if not nx.is_connected(migration):
            # if not, keep only the giant component
            largest_cc = max(nx.connected_components(migration), key=len)
            migration = migration.subgraph(largest_cc).copy()

        # add the resulting graph to the list
        migration_list.append(migration.copy())

    return migration_list

