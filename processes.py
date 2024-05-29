import random
import statistics

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from collections import OrderedDict

from networkx.algorithms.community import girvan_newman

import networkx as nx
import random

def remove_edge_random(net: nx.Graph) -> list:
    """
    Iteratively removes a random edge from the input NetworkX graph until only two edges remain.
    Keeps track of all intermediate graphs and returns them in a list.

    Parameters:
    net (nx.Graph): The initial NetworkX graph object.

    Returns:
    list: A list of NetworkX graph objects representing the state of the graph after each edge removal.
    """
    # Copy the input graph to avoid modifying the original
    migration = net.copy()
    # Initialize a list to store the intermediate graphs
    migration_list = [migration.copy()]
    # Keep track of the number of edges to avoid recalculating in each iteration
    edge_count = migration.number_of_edges()

    # Continue removing edges until only two remain
    while edge_count > 2:
        # Get the list of current edges in the graph
        edges = list(migration.edges())
        # Select a random edge to remove
        edge = random.choice(edges)
        # Remove the selected edge from the graph
        migration.remove_edge(*edge)
        # Append a copy of the current state of the graph to the list
        migration_list.append(migration.copy())
        # Decrement the edge count
        edge_count -= 1

    return migration_list


def remove_edge_intrusive(net: nx.Graph) -> list:
    """
    Iteratively removes edges connected to a randomly chosen node until that node is isolated.
    Continues the process by selecting new random nodes until only two edges remain.
    Keeps track of all intermediate graphs and returns them in a list.

    Parameters:
    net (nx.Graph): The initial NetworkX graph object.

    Returns:
    list: A list of NetworkX graph objects representing the state of the graph after each edge removal.
    """
    # Copy the input graph to avoid modifying the original
    migration = net.copy()
    # Initialize a list to store the intermediate graphs
    migration_list = [migration.copy()]
    # Get a list of all nodes in the graph
    nodes = list(migration.nodes)
    # Keep track of the number of edges to avoid recalculating in each iteration
    edge_count = migration.number_of_edges()

    # Continue removing edges until only two remain
    while edge_count > 2:
        # Choose a random node
        node = random.choice(nodes)
        # Get all edges connected to the chosen node
        edges = list(migration.edges(node))

        # Remove the node from the nodes list to avoid reselecting it
        nodes.remove(node)

        # Remove edges connected to the node until it's isolated
        while edges and edge_count > 2:
            edge = random.choice(edges)
            migration.remove_edge(*edge)
            edges.remove(edge)
            migration_list.append(migration.copy())
            edge_count -= 1

    return migration_list


net = nx.random_geometric_graph(10,0.5)
pos=nx.spring_layout(net,seed=4)
nx.draw_networkx(net,pos)
plt.show()

x=remove_edge_intrusive(net)

nx.draw_networkx(x[5],pos)
plt.show()

nx.draw_networkx(x[10],pos)
plt.show()

nx.draw_networkx(x[15],pos)
plt.show()

nx.draw_networkx(x[18],pos)
plt.show()
nx.draw_networkx(x[22],pos)
plt.show()


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




def get_edges_connected_to_nodes(G, nodes_list):
    """
    extract the edges connected to a list of nodes and returen an ordered list.
    used in remove edge regressive.
    :param G: network
    :param nodes_list: list of nodes ordered by distance from the side
    :return:
    """

    # Use an OrderedDict to store the edges without duplicates
    edges_ordered_dict = OrderedDict()

    # Iterate through the nodes and add the edges to the OrderedDict
    for node in nodes_list:
        for edge in G.edges(node):
            # The order of the nodes in the edge might not be consistent, so we use a tuple with sorted nodes
            edges_ordered_dict[tuple(sorted(edge))] = None

    # Convert the keys of the OrderedDict to a list
    edges_list = list(edges_ordered_dict.keys())

    return edges_list




def remove_edge_regressive(net: nx.Graph) -> list:
    """
    Remove edge from the side\edge of the migration network
    :param net:  initial migration networkx network
    :return: list of networks with edges removed
    """

    pos = nx.spring_layout(net.copy(), seed=55)

    migration = net.copy()
    migration_list = []  # initialize list of networks

    # Get the positions of the nodes
    positions = nx.get_node_attributes(migration, 'pos')

    # Calculate the distance to the left side (x = 1) for each node
    distances = {node: 1 - pos[0] for node, pos in positions.items()}

    # create a list of nodes ordered by distance
    nodes = sorted(distances, key=distances.get, reverse=True)

    edges = get_edges_connected_to_nodes(migration,nodes)

    while nx.number_of_edges(migration) > 2:  # stop when the network includes only two connected nodes

        # Remove the longest edge from network
        migration.remove_edge(*edges[0])

        # Remove the longest edge from the edges list
        edges.pop(0)

        # add the resulting graph to the list
        migration_list.append(migration.copy())

    return migration_list


def find_edges_crossed_by_line(net):
    """
    create a path that crosses the network and get all the edges
    that will be fragmented.
    :param net: network
    :return: ordered list of the edges to remove
    """

    # Step 2: Define random endpoints of the straight line path
    start_point = np.random.rand(2)
    end_point = np.random.rand(2)

    # Step 3: Calculate equation of the line
    delta = end_point - start_point
    slope = delta[1] / delta[0]
    intercept = start_point[1] - slope * start_point[0]

    # Step 4: Check intersection with edges
    edges_on_line = []
    edge_positions = []
    for edge in net.edges():
        p1 = np.array(net.nodes[edge[0]]['pos'])
        p2 = np.array(net.nodes[edge[1]]['pos'])

        # Determine the sides of the line for both points
        side_p1 = slope * p1[0] - p1[1] + intercept
        side_p2 = slope * p2[0] - p2[1] + intercept

        # Check if points are on opposite sides of the line to know the order (left-right; up-down)
        if side_p1 * side_p2 < 0:
            edges_on_line.append(edge)
            # Store average position for ordering
            if -1 < slope < 1:
                edge_positions.append((p1[0] + p2[0]) / 2)
            else:
                edge_positions.append((p1[1] + p2[1]) / 2)

    # Order the edges along the line based on average x or y positions
    ordered_edges_on_line = [edges_on_line[i] for i in np.argsort(edge_positions)]

    return ordered_edges_on_line


def remove_edge_divisive(net: nx.Graph) -> list:
    """
    Remove edges from a network based on a path of fragmentation
    :param net:  initial migration networkx network
    :return: list of networks with edges removed
    """

    migration = net.copy()
    migration_list = []  # initialize list of networks

    while nx.number_of_edges(migration) > 2:  # stop when the network includes only two connected nodes

        # get a list of edges that lie on the fragmented path
        edges = find_edges_crossed_by_line(migration)

        for edge in range(len(edges)): # for each edge from the current fragmented path

            # Remove edge from network
            migration.remove_edge(*edges[0])

            # Remove the longest edge from the edges list
            edges.pop(0)

            # add the resulting graph to the list
            migration_list.append(migration.copy())

    return migration_list


def remove_edge_optimal(net: nx.Graph) -> list:
    """
    Remove edges from the network to maximize connectivity between nodes for as long as possible.
    Track all intermediate states of the network in a list until the stopping condition is met.

    :param net: initial networkx object
    :return: list of networkx objects showing the network's evolution
    """
    migration = net.copy()
    migration_list = [migration.copy()]  # start with the original network

    while nx.number_of_edges(migration) > 2:  # Adjust the condition as needed
        # Compute edge betweenness centrality to determine the importance of edges
        centrality = nx.edge_betweenness_centrality(migration)

        # Find the edge with the lowest centrality (minimal impact on connectivity)
        edge_to_remove = min(centrality, key=centrality.get)

        # Remove the selected edge
        migration.remove_edge(*edge_to_remove)

        # Add the resulting graph to the list
        migration_list.append(migration.copy())
    return migration_list


def remove_edge_optimal_no_update(net: nx.Graph) -> list:
    """
    Remove edges from the network to maximize connectivity between nodes.
    unlike optimal, this process calculates the edges importance at the begining once, without updating.
    Track all intermediate states of the network in a list until the stopping condition is met.

    :param net: initial networkx object
    :return: list of networkx objects showing the network's evolution
    """
    migration = net.copy()
    migration_list = [migration.copy()]  # start with the original network

    # Compute edge betweenness centrality to determine the importance of edges
    centrality = nx.edge_betweenness_centrality(migration)
    edges = sorted(centrality, key=centrality.get)

    while nx.number_of_edges(migration) > 2:
        # Find the edge with the lowest centrality (minimal impact on connectivity)
        edge_to_remove = edges[0]

        # Remove the selected edge
        migration.remove_edge(*edge_to_remove)

        # Add the resulting graph to the list
        migration_list.append(migration.copy())

        # Remove the edge with the highest betweenness from the edges list
        edges.pop(0)

    return migration_list


def remove_edge_worst(net: nx.Graph) -> list:
    """
    Remove edges from the network to minimze connectivity between nodes.
    Track all intermediate states of the network in a list until the stopping condition is met.

    :param net: initial networkx object
    :return: list of networkx objects showing the network's evolution
    """
    migration = net.copy()
    migration_list = [migration.copy()]  # start with the original network

    while nx.number_of_edges(migration) > 2:  # Adjust the condition as needed
        # Compute edge betweenness centrality to determine the importance of edges
        centrality = nx.edge_betweenness_centrality(migration)

        # Find the edge with the lowest centrality (minimal impact on connectivity)
        edge_to_remove = max(centrality, key=centrality.get)

        # Remove the selected edge
        migration.remove_edge(*edge_to_remove)

        # Add the resulting graph to the list
        migration_list.append(migration.copy())
    return migration_list

net1 = nx.random_geometric_graph(10,0.6,seed=7)

remove_edge_optimal_no_update(net1)
# def remove_edge_random_giant_comp(net: nx.Graph) -> list:
#     """
#     Remove a random edge from the input network in each iteration.
#     If the resulting graph is disconnected, it only keeps the giant component.
#     Keep track of all generated graphs in a list until only two nodes remain.
#
#     :param net: initial networkx object
#     :return: list of networkx objects
#     """
#     migration = net.copy()
#     migration_list = []  # initialize list with the original network
#     while len(migration) > 2:  # stop when the network includes only two connected nodes
#         # choose a random edge and remove it
#         edges = list(migration.edges())
#         edge = random.choice(edges)
#         migration.remove_edge(*edge)
#
#         # check if the resulting graph is connected
#         if not nx.is_connected(migration):
#             # if not, keep only the giant component
#             largest_cc = max(nx.connected_components(migration), key=len)
#             migration = migration.subgraph(largest_cc).copy()
#
#         # add the resulting graph to the list
#         migration_list.append(migration.copy())
#
#     return migration_list
#
#
# def remove_edge_correlated_giant_comp(net: nx.Graph) -> list:
#     """
#     Remove a correlated edge from migration network of type networkx
#     each iteration we choose an edge to remove from the list of edges connected to it
#     if the node is no longer part of the network we choose a random node
#     :param net: initial migration network
#     :param n: no. of fragmentation steps
#     :return: list of networks after n edge removal
#     """
#
#     migration = net.copy()
#     migration_list = [net.copy()]  # initialize list with the original network
#
#     # choose a random edge for start
#     edge = random.choice(list(migration.edges))
#
#     # take the initial nodes to remove edges from
#     edge_a = edge[0]
#     edge_b = edge[1]
#
#     while len(migration) > 2:  # stop when the network includes only two connected nodes
#         # print(nx.attr_matrix(migration))
#         # takes the edges of the nodes of the removed edge
#         edges_a = list(migration.edges(edge_a))
#         edges_b = list(migration.edges(edge_b))
#
#         # make a list of the edges of the two nodes
#         edges = edges_a + edges_b
#
#         # sample an edge from the edges
#         edge = random.choice(list(edges))
#
#         # remove the chosen edge from the network
#         migration.remove_edge(*edge)
#
#         # check if the resulting graph is connected
#         if not nx.is_connected(migration):
#             # if not, keep only the giant component
#             largest_cc = max(nx.connected_components(migration), key=len)
#             migration = migration.subgraph(largest_cc).copy()
#
#         # choose the nodes to remove edges from
#         # if node doesn't exist in new network choose a random node
#         if edge[0] in migration.nodes:
#             edge_a = edge[0]
#         else:
#             edge_a = random.choice(list(migration.nodes))
#
#         if edge[1] in migration.nodes:
#             edge_b = edge[1]
#         else:
#             edge_b = random.choice(list(migration.nodes))
#
#         # add the resulting graph to the list
#         migration_list.append(migration.copy())
#
#     return migration_list
#
#
# def remove_edge_distance_giant_comp(net: nx.Graph) -> list:
#     """
#     Remove edge based on distance from migration network of type networkx
#     :param net:  initial migration network
#     :return: list of networks after n edge removal
#     """
#
#     pos = nx.spring_layout(net.copy(), seed=55)
#
#     migration = net.copy()
#     migration_list = []  # initialize list of networks
#
#     while len(migration) > 2:  # stop when the network includes only two connected nodes
#
#         edges = migration.edges()
#
#         # Calculate the euclidean distance between all nodes and create a dict
#         distances = {edge: round(
#             ((pos[edge[1]][1] - pos[edge[0]][1]) ** 2 + (pos[edge[1]][0] - pos[edge[0]][0]) ** 2) ** 0.5, 2)
#             for edge in edges}
#
#         edges = sorted(edges, key=distances.get, reverse=True)
#
#         # Remove the longest edge from network
#         migration.remove_edge(*edges[0])
#
#         # Remove the longest edge from the edges list
#         edges.pop(0)
#
#         # check if the resulting graph is connected
#         if not nx.is_connected(migration):
#             # if not, keep only the giant component
#             largest_cc = max(nx.connected_components(migration), key=len)
#             migration = migration.subgraph(largest_cc).copy()
#
#         # add the resulting graph to the list
#         migration_list.append(migration.copy())
#
#     return migration_list
#

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

