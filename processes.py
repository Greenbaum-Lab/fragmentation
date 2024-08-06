
import statistics
import numpy as np
from matplotlib import pyplot as plt
from collections import OrderedDict
import networkx as nx
import random
from math import sqrt

def remove_edge_random(net: nx.Graph) -> list:
    """
    Iteratively removes a random edge from the input NetworkX network until no edge remains.
    Keeps track of all intermediate networks and returns them in a list.

    Parameters:
    net (nx.Graph): The initial NetworkX network object.

    Returns:
    list: A list of networks after each edge removal - from initial network to isolated nodes.
    """
    # Copy the input graph to avoid modifying the original
    migration = net.copy()
    # Initialize a list to store all networks
    migration_list = [migration.copy()]
    # Keep track of the number of edges
    edge_count = migration.number_of_edges()

    # Continue removing edges until only two remain
    while edge_count:
        # Get the list of current edges in the graph
        edges = list(migration.edges())
        # Select a random edge to remove
        edge = random.choice(edges)
        # Remove the selected edge from the graph
        migration.remove_edge(*edge)
        # Append a copy of the current network to the list
        migration_list.append(migration.copy())
        # Decrement the edge count
        edge_count -= 1

    return migration_list


def remove_edge_intrusive(net: nx.Graph) -> list:
    """
    Iteratively removes edges connected to a randomly chosen node until that node is isolated.
    Continues the process by selecting new random nodes until no edge remains.
    Keeps track of all intermediate networks and returns them in a list.

    Parameters:
    net (nx.Graph): The initial NetworkX network object.

    Returns:
    list: A list of networks after each edge removal - from initial network to isolated nodes.
    """
    # Copy the input graph to avoid modifying the original
    migration = net.copy()
    # Initialize a list to store the intermediate graphs
    migration_list = [migration.copy()]
    # Get a list of all nodes in the graph
    nodes = list(migration.nodes)
    # Keep track of the number of edges
    edge_count = migration.number_of_edges()

    # Continue removing edges until only two remain
    while edge_count:
        # Choose a random node
        node = random.choice(nodes)
        # Get all edges connected to the chosen node
        edges = list(migration.edges(node))
        # Remove the node from the nodes list to avoid reselecting it
        nodes.remove(node)
        # Remove edges connected to the node until it's isolated
        while edges and edge_count:
            edge = random.choice(edges)
            migration.remove_edge(*edge)
            edges.remove(edge)
            migration_list.append(migration.copy())
            edge_count -= 1

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



def remove_edge_correlated(net: nx.Graph) -> list:
    """
    Iteratively removes an edge from migration network.
    We begin by randomly choosing an edge and removing it.
    Then, the nodes that are connected by this edge are chosen and their edges are stored.
    Among the edges stored, we randomly choose another edge and remove it.
    We repeat this process until the netowrk is completely broken
    each iteration we choose an edge to remove from the list of edges connected to it

    Parameters:
    net (nx.Graph): The initial NetworkX network object.

    Returns:
    list: A list of networks after each edge removal - from initial network to isolated nodes.
    """

    migration = net.copy()
    migration_list = [net.copy()]  # initialize list with the original network
    # choose a random edge for start
    edge = random.choice(list(migration.edges))
    # take the initial nodes to remove edges from
    node_a = edge[0]
    node_b = edge[1]
    while nx.number_of_edges(migration):
        # takes the edges of the nodes of the removed edge

        edges_a = list(migration.edges(node_a))
        edges_b = list(migration.edges(node_b))
        # make a list of the edges of the two nodes
        edges = edges_a + edges_b
        # if nodes doesn't have more connected edges choose a random edge
        if not edges:
            connected_nodes = get_connected_nodes(migration)
            connected_edges = get_connected_edges(migration, connected_nodes)
            edge = random.choice(list(connected_edges))

        else:
            edge = random.choice(edges)

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
    Remove edge based on distance.
    Calculate the position of nodes and their pairwise distance.
    Remove edges from the largest distance to smaller.
    Parameters:
    net (nx.Graph): The initial NetworkX network object.
    Returns:
    list: A list of networks after each edge removal - from initial network to isolated nodes.
    """

    migration = net.copy()
    migration_list = []

    pos = nx.get_node_attributes(migration, 'pos')
    edges = migration.edges()

    # Calculate the Euclidean distance between all nodes and create a dict
    distances = {}
    for node in edges:
        x_dist = ((pos[node[1]][1] - pos[node[0]][1]) ** 2)
        y_dist = ((pos[node[1]][0] - pos[node[0]][0]) ** 2)
        eucl_dist = sqrt(x_dist + y_dist)
        distances.update({node: eucl_dist})

    # sort edges by distance
    edges = sorted(edges, key=distances.get, reverse=True)

    while edges:
        # Remove the longest edge from network (first item)
        migration.remove_edge(*edges[0])
        # Remove the longest edge from the edges list
        edges.pop(0)
        migration_list.append(migration.copy())

    return migration_list



def sort_edges(edges,nodes_order):
    edge_dict = {}
    result_edges = []
    for edge in edges:
        dict_key = nodes_order.index(edge[1])
        edge_dict.update({dict_key: edge})
    # Sort the dictionary by keys
    sorted_keys = sorted(edge_dict.keys())
    # Loop through the sorted keys and append the corresponding tuples to the list
    for key in sorted_keys:
        result_edges.append(edge_dict[key])
    return result_edges


def remove_edge_regressive(net: nx.Graph) -> list:
    """
    Remove edge from the side\edge of the migration network.
    Get the x coordinates of all nodes and sort them from 0 to 1.
    Remove all edges that are connected the first node (no specfic order).
    Parameters:NEED TO IMPROVE
    net (nx.Graph): The initial NetworkX network object.
    Returns:
    list: A list of networks after each edge removal - from initial network to isolated nodes.
    """

    migration = net.copy()
    pos = nx.get_node_attributes(migration, 'pos')
    nodes = list(pos.keys())
    migration_list = []

    # Get the x coordinates for each node
    distances = {}
    for node in nodes:
        dist = pos[node][0]
        distances.update({node: dist})
    # create a list of nodes ordered by distance from the edge
    nodes_order = sorted(distances,key=distances.get)
    nodes_list = sorted(distances,key=distances.get)
    while nx.number_of_edges(migration):
        #get the edges of the first node (the node at the edge)
        edges = list(migration.edges(nodes_list[0]))
        sorted_edges = sort_edges(edges=edges, nodes_order=nodes_order)
        # Remove the node from nodes list if it has no edges
        if not sorted_edges:
            nodes_list.pop(0)
            continue
        # Remove the edge from network
        migration.remove_edge(*sorted_edges[0])
        # Remove the edge from the edges list
        sorted_edges.pop(0)
        # add the resulting graph to the list
        migration_list.append(migration.copy())

    return migration_list


def find_intersection(x1, y1, x2, y2, x3, y3, x4, y4):
    """
    Find interesection between two line segments. if no interesection return None.
    used in remove edge divisive.
    """

    # Calculate the denominators
    denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if denominator == 0:
        return None  # The lines are parallel or coincident

    # Calculate the intersection point
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denominator
    u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denominator

    if 0 <= t <= 1 and 0 <= u <= 1:
        intersection_x = x1 + t * (x2 - x1)
        intersection_y = y1 + t * (y2 - y1)
        return (intersection_x, intersection_y)
    else:
        return None  # The intersection point is not within the bounds of the line segments


def make_edges_lines(net):
    """
    create edges segments from nodes positions.
    used in remove edge divisive.
    :param net:
    :return:
    """
    edges_pos = {}
    for edge in net.edges():
        p1 = tuple(net.nodes[edge[0]]['pos'])
        p2 = tuple(net.nodes[edge[1]]['pos'])
        edge_line = [p1, p2]
        edges_pos.update({edge: edge_line})

    return edges_pos


def generate_divisive_line():
    while True:
        # Line coordinates must be within 0 and 1 (the edges of the metric)
        coords = [0, 1, random.random(), random.random()]
        random.shuffle(coords)
        coord_x = tuple(coords[:2])
        coord_y = tuple(coords[2:])

        if coord_y == (1, 0) or coord_y == (0, 1) or coord_x == (1, 0) or coord_x == (0, 1):
            continue

        return [coord_x, coord_y]


def  remove_edge_divisive(net: nx.Graph) -> list:
    """
    Remove edges from a network based on a linear path crossing the network.
    genrate random line and make all edges line segments.
    find intersection between divisive line and edges.
    if there is interesection, remove edges based on their x coorrdinates.
    :return: list of networks with edges removed
    """

    migration = net.copy()
    migration_list = []

    edges_lines = make_edges_lines(migration)
    intersections = {}
    pos = nx.get_node_attributes(net, 'pos')

    while edges_lines:
        div_line = generate_divisive_line()
        for edge_nodes, edge in edges_lines.items():
            intersect = find_intersection(x1=div_line[0][0], y1=div_line[0][1],
                                          x2=div_line[1][0], y2=div_line[1][1],
                                          x3=edge[0][0], y3=edge[0][1],
                                          x4=edge[1][0], y4=edge[1][1])
            if intersect:
                intersections.update({edge_nodes:intersect})

        edges_to_remove = sorted(intersections,key=intersections.get)

        for edge in edges_to_remove:
            nx.draw_networkx(migration,pos=pos,with_labels=True)
            plt.show()
            # Remove the selected edge from the graph
            migration.remove_edge(*edge)
            # Append a copy of the current network to the list
            migration_list.append(migration.copy())
            # remove edge from edge dict
            edges_lines.pop(edge)

        intersections.clear()
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




# check for copy type if it overwrote the mutubale objects in fragmentatioon