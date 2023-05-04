import networkx as nx
import random
import math
from matplotlib import pyplot as plt


def remove_edge_random(migration, n: int) -> list:
    """
    Remove a random edge from migration network of type networkx
    :param migration:  initial migration network
    :param n: no. of fragmentation steps
    :return: list of networks after n edge removal
    """
    migration_list = [migration]
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    pos = nx.spring_layout(migration, seed=55)

    for i in range(n):
        if i == 0:
            nx.draw_networkx(migration, pos=pos, ax=axs[i])
            axs[i].set_title(f"Original network")
        if not nx.is_connected(migration_list[i]):  # stop when network breaks and plot last network
            nx.draw_networkx(migration_list[i - 2], pos=pos, ax=axs[1])
            axs[1].set_title(f"Network after {i} edges removed randomly")
            break
        edges = list(nx.edges(migration))
        edges_to_remove = (random.sample(edges, k=1))  # choose a random edge
        migration.remove_edge(*(edges_to_remove[0]))
        print(nx.is_connected(migration))
        migration_list.append(migration.copy())
    plt.show()
    return migration_list



def remove_edge_correlated(migration, n: int) -> list:
    """
    Remove a random edge from migration network of type networkx
    :param migration:  initial migration network
    :param n: no. of fragmentation steps
    :return: list of networks after n edge removal
    """
    migration_list = []
    prob_increase = 0.8
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    pos = nx.spring_layout(migration, seed=55)

    for i in range(n):
        if i == 0:
            nx.draw_networkx(migration, pos=pos, ax=axs[i])
            axs[i].set_title(f"Original network")
        if not nx.is_connected(migration):  # stop when network breaks and plot last network
            nx.draw_networkx(migration_list[i - 2], pos=pos, ax=axs[1])
            axs[1].set_title(f"{i} edges removed by correlataion ")
            break
        # make a list of all edges
        edges = list(migration.edges())  # create a list of all edges in the network

        # choose a random edge from the remaining edges
        edge = random.choice(edges)

        # remove the chosen edge from the network
        migration.remove_edge(*edge)

        # increase the probability of removing another edge from the same node
        if random.random() < prob_increase:
            node = random.choice(list(migration.nodes()))
            edges = [e for e in edges if e[0] != node and e[1] != node]

        migration_list.append(migration.copy())

    plt.show()
    return migration_list


def remove_edges_distance(migration, n: int) -> list:
    """
    Remove edge from migration network of type networkx
    :param migration:  initial migration network
    :param n: no. of fragmentation steps
    :return: list of networks after n edge removal
    """
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    migration_list = []
    pos = nx.spring_layout(migration, seed=55)
    edges = migration.edges()  # Get all the edges of the graph
    distances = {edge: round(
        ((pos[edge[1]][1] - pos[edge[0]][1]) ** 2 + (pos[edge[1]][0] - pos[edge[0]][0]) ** 2) ** 0.5, 2)
                 for edge in edges}  # Calculate the euclidean distance between all nodes

    # Sort the edges by their distances in descending order
    edges = sorted(edges, key=distances.get, reverse=True)

    for i in range(n):
        if i == 0:
            nx.draw_networkx(migration, pos=pos, ax=axs[i])
            axs[i].set_title(f"Original network")
        if not nx.is_connected(migration):  # stop when network breaks and plot last network
            nx.draw_networkx(migration_list[i - 2], pos=pos, ax=axs[1])
            axs[1].set_title(f"{i} edges removed by distance ")
            break

        else:
            # Remove the longest edge
            migration.remove_edge(*edges[0])
            print(nx.is_connected(migration))
            migration_list.append(migration.copy())

            # Remove the longest edge from the edges list
            edges.pop(0)

    plt.show()
    return migration_list
