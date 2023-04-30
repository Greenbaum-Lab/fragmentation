import networkx as nx
import random
import math
from matplotlib import pyplot as plt




def remove_edge_random(net: nx.Graph, n: int) -> list:
    """
    Remove a random edge from migration network of type networkx
    :param migration:  initial migration network
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
            axs[1].set_title(f"Network after {i} edges removed randomly")

            break

    plt.show()
    return migration_list


def remove_edge_correlated(net: nx.Graph, n: int) -> list:
    """
    Remove a random edge from migration network of type networkx
    :param net: initial migration network
    :param n: no. of fragmentation steps
    :return: list of networks after n edge removal
    """
    migration_list = [net.copy()]
    prob_increase = 0.8
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    pos = nx.spring_layout(net.copy(), seed=55)

    for i in range(n):

        # takes the last item in list
        migration = migration_list[-1].copy()

        #  plot initial network
        if i == 0:
            nx.draw_networkx(migration, pos=pos, ax=axs[i])
            axs[i].set_title(f"Original network")
        # if not nx.is_connected(migration):  # stop when network breaks and plot last network
        #     nx.draw_networkx(migration_list[i - 2], pos=pos, ax=axs[1])
        #     axs[1].set_title(f"{i} edges removed by correlataion ")
        #     break
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
