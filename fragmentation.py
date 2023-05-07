import copy

import networkx as nx
import random
import math
from matplotlib import pyplot as plt


# def remove_edge_random(net: nx.Graph, n: int) -> list:
#     """
#     Remove a random edge from migration network of type networkx
#     :param net:  initial migration network
#     :param n: no. of fragmentation steps
#     :return: list of networks after n edge removal
#     """
#     migration_list = [net.copy()]
#     fig, axs = plt.subplots(1, 2, figsize=(10, 5))
#     pos = nx.spring_layout(net.copy(), seed=55)
#
#     for i in range(n):
#         migration = migration_list[-1].copy()  # takes the last item in list
#
#         #  plot initial network
#         if i == 0:
#             nx.draw_networkx(migration, pos=pos, ax=axs[i])
#             axs[i].set_title(f"Original network")
#
#         edges = list(nx.edges(migration))
#         edges_to_remove = (random.sample(edges, k=1))  # choose a random edge
#         migration.remove_edge(*(edges_to_remove[0]))
#
#         # stop when network breaks
#         if nx.is_connected(migration):
#             migration_list.append(migration)
#
#         else:
#
#             #  plot last network before breaking
#             nx.draw_networkx(migration_list[i], pos=pos, ax=axs[1])
#             axs[1].set_title(f"{i} edges removed randomly")
#
#             break
#
#     plt.show()
#     return migration_list


# def remove_edge_correlated(net: nx.Graph, n: int) -> list:
#     """
#     Remove a correlated edge from migration network of type networkx
#     :param net: initial migration network
#     :param n: no. of fragmentation steps
#     :return: list of networks after n edge removal
#     """
#     migration_list = [net.copy()]
#     prob_increase = 0.99
#     fig, axs = plt.subplots(1, 2, figsize=(10, 5))
#     pos = nx.spring_layout(net.copy(), seed=55)
#
#     for i in range(n):
#
#         # takes the last item in list
#         migration = migration_list[-1].copy()
#
#         #  plot initial network
#         if i == 0:
#             nx.draw_networkx(migration, pos=pos, ax=axs[i])
#             axs[i].set_title(f"Original network")
#         # make a list of all edges
#         edges = list(migration.edges())  # create a list of all edges in the network
#
#         # choose a random edge from the remaining edges
#         edge = random.choice(edges)
#
#         # remove the chosen edge from the network
#         migration.remove_edge(*edge)
#
#         # increase the probability of removing another edge from the same node
#         if random.random() < prob_increase:
#             node = random.choice(list(migration.nodes()))
#             edges = [e for e in edges if e[0] != node and e[1] != node]
#
#         if nx.is_connected(migration):
#             migration_list.append(migration)
#         else:
#
#             #  plot last network before breaking
#             nx.draw_networkx(migration_list[i], pos=pos, ax=axs[1])
#             axs[1].set_title(f"{i} edges removed by correlation")
#             break
#
#     plt.show()
#     return migration_list
#

#
# net = nx.erdos_renyi_graph(n=5, p=0.8, seed=5)
# print(net.edges)
# edge = random.choice(list(net.edges))
# print(edge[0])
# print(edge[1])
# print(edge)
# print(net.edges(0))
# print(net.edges(2))
# edges_b = list(net.edges(2))
# edges_a = list(net.edges(0))
#
# edges = edges_a + edges_b
#
# print(edges)


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


def remove_edge_random(net: nx.Graph, n: int) -> list:
    """
    Remove a random edge from migration network of type networkx
    :param net:  initial migration network
    :param n: no. of fragmentation steps
    :return: list of networks after n edge removal
    """
    migration_list = [net.copy()]

    for i in range(n):
        migration = migration_list[-1].copy()  # takes the last item in list

        edges = list(nx.edges(migration))
        edges_to_remove = (random.sample(edges, k=1))  # choose a random edge
        migration.remove_edge(*(edges_to_remove[0]))

        migration_list.append(migration)

    return migration_list


# def remove_edge_random_gc(net: nx.Graph, n: int) -> list:
#     """
#     Remove a random edge from migration network of type networkx
#     :param net:  initial migration network
#     :param n: no. of fragmentation steps
#     :return: list of networks after n edge removal
#     """
#     migration_list = [net.copy()]
#
#     for i in range(n):
#         giant_net = migration_list[-1]
#         edges = list(nx.edges(giant_net))
#         edges_to_remove = random.sample(edges, k=1)
#         giant_net.remove_edge(*edges_to_remove[0])
#         edges.remove(edges_to_remove[0])
#         if nx.is_connected(giant_net):
#             migration_list.append(giant_net.copy())
#
#         else:
#             # if len(giant_net) > 2:
#             compon = max(nx.connected_components(giant_net), key=len)
#             giant_component = giant_net.subgraph(compon)
#             if nx.is_connected(giant_component):
#                 migration_list.append(giant_component.copy())
#
#     return migration_list


import networkx as nx
import random


def remove_edge_random_gc(net: nx.Graph) -> list:
    """
    Remove a random edge from a networkx graph
    :param net: initial graph
    :param n: number of edge removals
    :return: list of connected graphs after edge removal
    """
    migration = net.copy()
    migration_list = []
    while nx.is_connected(migration):
        migration_list.append(migration.copy())

        edges = list(nx.edges(migration))
        edges_to_remove = random.sample(edges, k=1)
        migration.remove_edge(*edges_to_remove[0])
        print(migration_list[-1].edges)
        print(nx.is_connected(migration_list[-1]))

        if not nx.is_connected(migration):

            components = sorted(nx.connected_components(migration), key=len, reverse=True)
            giant_component = migration.subgraph(components[0])
            migration_list.append(giant_component.copy())
            print(migration_list[-1].edges)
            print(nx.is_connected(migration_list[-1]))
            # if len(giant_component) < 2:
            #     break

    return migration_list


net = nx.erdos_renyi_graph(10, 0.8,seed=54)
nx.draw_networkx(net)
plt.show()

print(len(net))
#
# nx.draw_networkx(x[70])
# plt.show()

y = remove_edge_random_gc(net)
print(len(y[-1]))
print(len(y))