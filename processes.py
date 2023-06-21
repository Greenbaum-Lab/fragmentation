import random
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt


def remove_edge_random(net: nx.Graph) -> list:
    """
    Remove a random edge from the input network in each iteration.
    Keep track of all generated graphs in a list until only two nodes remain.

    :param net: initial networkx object
    :return: list of networkx objects
    """
    migration = net.copy()
    migration_list = []  # initialize list with the original network

    while nx.number_of_edges(migration) > 2:  # stop when the network includes only two connected nodes
        # choose a random edge and remove it
        edges = list(migration.edges())
        edge = random.choice(edges)
        migration.remove_edge(*edge)

        # add the resulting graph to the list
        migration_list.append(migration.copy())

    return migration_list


def remove_edge_correlated(net: nx.Graph) -> list:
    """
    Remove a correlated edge from migration network of type networkx
    each iteration we choose an edge to remove from the list of edges connected to it
    :param net: initial migration network
    :return: list of networks after n edge removal
    """

    migration = net.copy()
    migration_list = [net.copy()]  # initialize list with the original network

    # choose a random edge for start
    edge = random.choice(list(migration.edges))

    # take the initial nodes to remove edges from
    node_a = edge[0]
    node_b = edge[1]

    while nx.number_of_edges(migration) > 2:  # stop when the network includes only two connected nodes
        # takes the edges of the nodes of the removed edge
        edges_a = list(migration.edges(node_a))
        edges_b = list(migration.edges(node_b))

        # make a list of the edges of the two nodes
        edges = edges_a + edges_b

        # if nodes doesn't have more connected edges choose a random edge
        if len(edges) < 1:
            connected_nodes = get_connected_nodes(migration)
            connected_edges = get_connected_edges(migration, connected_nodes)
            edge = random.choice(list(connected_edges))

        else:
            # sample an edge from the edges
            edge = random.choice(list(edges))

        # remove the chosen edge from the network
        migration.remove_edge(*edge)

        # choose the nodes to remove edges from
        # if node doesn't exist in new network choose a random node
        node_a = edge[0]
        node_b = edge[1]

        # add the resulting graph to the list
        migration_list.append(migration.copy())

    return migration_list


def remove_edge_distance(net: nx.Graph) -> list:
    """
    Remove edge based on distance from migration network of type networkx
    :param net:  initial migration network
    :return: list of networks after n edge removal
    """

    pos = nx.spring_layout(net.copy(), seed=55)

    migration = net.copy()
    migration_list = []  # initialize list of networks

    while nx.number_of_edges(migration) > 2:  # stop when the network includes only two connected nodes

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

        # add the resulting graph to the list
        migration_list.append(migration.copy())

    return migration_list


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


def get_connected_nodes(net: nx.Graph) -> set:
    """
    Get all the nodes in a network that are connected (exclude isolated nodes).
    used in  correlated fragmentation

    :param net: networkx graph object
    :return: set of connected nodes
    """
    connected_nodes = set()
    components = nx.connected_components(net)

    for component in components:
        if len(component) > 1:  # Exclude isolated nodes
            connected_nodes.update(component)

    return connected_nodes


def get_connected_edges(net: nx.Graph, connected_nodes: set) -> list:
    """
    Get all the edges that are connected to the specified connected nodes.
    used in  correlated fragmentation

    :param net: networkx graph object
    :param connected_nodes: set of connected nodes
    :return: list of connected edges
    """
    connected_edges = []

    for u, v in net.edges():
        if u in connected_nodes or v in connected_nodes:
            connected_edges.append((u, v))

    return connected_edges


def intervals(lst):
    if len(lst) <= 50:
        return lst

    interval = max((len(lst) - 1) // 49, 1)
    return lst[:49 * interval:interval] + [lst[-1]]


# def find_breaking_point(lst):
#     """
#     find the index of the list where the network is no longer connected
#     """
#     first_element = lst[0]
#
#     for i, element in enumerate(lst):
#         if len(element) < len(first_element):
#             return i
#
#     return -1  # Return -1 if no element is found


def find_breaking_point(networks):
    """
    find the index of the list where the network is no longer connected
    """
    for index, network in enumerate(networks):
        if not nx.is_connected(network):
            return index
    return None

def find_breakink_point_list(networks: list):
    breaking_point = []
    for net in networks:
        x = find_breaking_point(net)
        breaking_point.append(x)
    return breaking_point

# # random.seed(6)
# # net = nx.barabasi_albert_graph(100,2,seed=2)
# # Create a new figure and axis
# fig, ax = plt.subplots(figsize=(20, 20))
# pos = nx.spring_layout(net, seed=98)  # set the fixed position for plotting the network
# # Draw the networkx graph on the axis
# nx.draw_networkx(net, with_labels=False, ax=ax, width=0.5, pos=pos)
# # plt.show()
# plt.savefig("main_AB.svg", format="svg")
#
# rand = remove_edge_random(net)
# rand = intervals(rand)
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
# nx.draw_networkx(rand[10], with_labels=False, ax=ax2, width=0.5, pos=pos)
# nx.draw_networkx(rand[18], with_labels=False, ax=ax1, width=0.5, pos=pos)
# # plt.show()
# plt.savefig("rand_AB.svg", format="svg")
#
# cor = remove_edge_correlated(net)
# cor = intervals(cor)
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
# nx.draw_networkx(cor[10], with_labels=False, ax=ax2, width=0.5, pos=pos)
# nx.draw_networkx(cor[18], with_labels=False, ax=ax1, width=0.5, pos=pos)
# # plt.show()
# plt.savefig("cor_AB.svg", format="svg")
#
# dist = remove_edge_distance(net)
# dist = intervals(dist)
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
# nx.draw_networkx(dist[10], with_labels=False, ax=ax2, width=0.5, pos=pos)
# nx.draw_networkx(dist[18], with_labels=False, ax=ax1, width=0.5, pos=pos)
# # plt.show()
# plt.savefig("dist_AB.svg", format="svg")
