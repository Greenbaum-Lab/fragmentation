
import statistics
import numpy as np
from matplotlib import pyplot as plt
from collections import OrderedDict
import networkx as nx
import random
from math import sqrt
def _remove_edge_pair(G: nx.Graph, u, v):
    """
    Remove the directed pair u→v and v→u *if* they exist.
    Works for Graph and DiGraph, raises no error if one
    direction is already absent.
    """
    if G.has_edge(u, v):
        G.remove_edge(u, v)
    if G.has_edge(v, u):
        G.remove_edge(v, u)


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
    pos = nx.get_node_attributes(migration, 'pos')
    # Continue removing edges until only two remain
    while edge_count:
        # nx.draw_networkx(net, pos)
        # plt.show()

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



import random
import networkx as nx
from typing import List, Set


def remove_edge_random(net: nx.DiGraph) -> List[nx.DiGraph]:
    """
    Iteratively removes a *pair* of opposite directed edges (u→v and v→u)
    chosen at random until the graph has no edges left.

    Parameters
    ----------
    net : nx.DiGraph
        The starting migration network.  Node and graph attributes are
        preserved in every snapshot.

    Returns
    -------
    List[nx.DiGraph]
        Snapshots of the network after each edge-pair removal, starting
        with the original graph and ending with isolated nodes.
    """
    migration = net.copy()
    all_migration = [migration.copy()]

    while migration.number_of_edges():
        # sample one existing *directed* edge
        u, v = random.choice(list(migration.edges()))
        # remove that direction
        migration.remove_edge(u, v)
        migration.remove_edge(v, u)

        all_migration.append(migration.copy())

    return all_migration


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
            _remove_edge_pair(migration, edge[0], edge[1])
            edges.remove(edge)
            migration_list.append(migration.copy())
            edge_count -= 1

    return migration_list


def get_connected_nodes(net: nx.Graph) -> Set:
    """
    Return every node that belongs to a size > 1 component.
    Works for:
      • nx.Graph / nx.MultiGraph  → uses nx.connected_components
      • nx.DiGraph / nx.MultiDiGraph → uses nx.weakly_connected_components
        (i.e. connectivity is evaluated while ignoring edge direction).
    """
    if net.is_directed():
        components = nx.weakly_connected_components(net)
    else:
        components = nx.connected_components(net)

    connected_nodes = {node
                       for comp in components if len(comp) > 1
                       for node in comp}

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
        # migration.remove_edge(*edge)
        _remove_edge_pair(migration, edge[0], edge[1])

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
        # migration.remove_edge(*edges[0])
        _remove_edge_pair(migration, edges[0][0], edges[0][1])

        # Remove the longest edge from the edges list
        edges.pop(0)
        migration_list.append(migration.copy())

    return migration_list


def sort_edges(edges,nodes_order):
    """Sort edges based on the order of nodes."""
    result_edges = []
    # Create a dictionary with the edges as values and the index of the second node as keys
    edge_dict = {nodes_order.index(edge[1]): edge for edge in edges}
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
    nodes_order = sorted(pos, key=lambda node: pos[node][0])
    migration_list = [migration.copy()]

    while migration.number_of_edges():
        for node in nodes_order:
            edges = list(migration.edges(node))
            if edges:
                sorted_edges = sort_edges(edges, nodes_order)
                # migration.remove_edge(*sorted_edges[0])
                _remove_edge_pair(migration, sorted_edges[0][0], sorted_edges[0][1])

                migration_list.append(migration.copy())
                break
        else:
            break

    return migration_list


def find_intersection(p1, p2, p3, p4):
    """
    Find intersection between two line segments. If no intersection, return None.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if denominator == 0:
        return None  # The lines are parallel or coincident

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

        if coord_y not in [(1, 0), (0, 1)] and coord_x not in [(1, 0), (0, 1)]:
            return [coord_x, coord_y]


def remove_edge_divisive(net: nx.Graph) -> list:
    """
    Remove edges from a network based on a linear path crossing the network.
    Generate a random line and make all edges line segments.
    Find intersection between divisive line and edges.
    If there is intersection, remove edges based on their x coordinates.
    :return: list of networks with edges removed
    """
    migration = net.copy()
    migration_list = []

    while migration.number_of_edges():
        div_line = generate_divisive_line()
        edges_lines = make_edges_lines(migration)
        intersections = {}

        for edge_nodes, edge in edges_lines.items():
            intersect = find_intersection(p1=edge[0], p2=edge[1],
                                          p3=div_line[0], p4=div_line[1])
            if intersect:
                intersections[edge_nodes] = intersect

        edges_to_remove = sorted(intersections, key=intersections.get)

        for edge in edges_to_remove:
            # migration.remove_edge(*edge)
            _remove_edge_pair(migration, edge[0], edge[1])

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

    while nx.number_of_edges(migration):  # Adjust the condition as needed
        # Compute edge betweenness centrality to determine the importance of edges
        centrality = nx.edge_betweenness_centrality(migration)

        # Find the edge with the lowest centrality (minimal impact on connectivity)
        edge_to_remove = min(centrality, key=centrality.get)

        # Remove the selected edge
        # migration.remove_edge(*edge_to_remove)
        _remove_edge_pair(edge_to_remove[0], edge_to_remove[1])
        # Add the resulting graph to the list
        migration_list.append(migration.copy())
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

    while nx.number_of_edges(migration):  # Adjust the condition as needed
        # Compute edge betweenness centrality to determine the importance of edges
        centrality = nx.edge_betweenness_centrality(migration)

        # Find the edge with the highest betweenness
        edge_to_remove = max(centrality, key=centrality.get)

        # Remove the selected edge
        # migration.remove_edge(*edge_to_remove)
        _remove_edge_pair(edge_to_remove[0], edge_to_remove[1])

        # Add the resulting graph to the list
        migration_list.append(migration.copy())
    return migration_list

