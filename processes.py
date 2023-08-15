import random
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from collections import OrderedDict


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



def remove_edge_intrusive(net: nx.Graph) -> list:
    """
    Remove edges from specific nodes until there are no more
    edges conneted to it. then choose anothe rrandom node ans islote him ...
    Keep track of all generated graphs in a list until only two nodes remain.

    :param net: initial networkx object
    :return: list of networkx objects
    """
    migration = net.copy()
    migration_list = []  # initialize list with the original network
    nodes = list(migration.nodes)

    while nx.number_of_edges(migration) > 2:  # stop when the network includes only two connected nodes

        # choose a random node
        node = random.choice(nodes)

        # get all the edges of the corresponding node
        edges = list(migration.edges(node))

        #update nodes list
        nodes.remove(node)

        for edge in range(len(edges)):

            edge = random.choice(edges)

            # update network and edges list
            migration.remove_edge(*edge)
            edges.remove(edge)

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




#
#
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


def intervals(lst):
    """
    take snapshots of the process
    :param lst:
    :return:
    """
    if len(lst) <= 50:
        return lst
    n = 19 # number of bins (-1)
    interval = max((len(lst) - 1) // n, 1)
    return lst[:n * interval:interval] + [lst[-1]]


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

# def create_networks(n, p_or_k, replicates=200):
#     er_edges = []
#     rgg_edges = []
#     ab_edges = []
#     sw_edges = []
#
#     er_density = []
#     rgg_density = []
#     ab_density = []
#     sw_density = []
#     connect = []
#     for _ in range(replicates):
#         # Erdos-Renyi (ER) Graph
#         er_graph = nx.erdos_renyi_graph(n, 0.2)
#         er_edges.append(er_graph.number_of_edges())
#         er_density.append(nx.density(er_graph))
#
#         # Random Geometric Graph (RGG)
#         rgg_graph = nx.random_geometric_graph(n, 0.3)
#
#         rgg_edges.append(rgg_graph.number_of_edges())
#         rgg_density.append(nx.density(rgg_graph))
#         connect.append(nx.is_connected(rgg_graph))
#         # Albert-Barabasi Graph
#         ab_graph = nx.barabasi_albert_graph(n, m=5)
#         ab_edges.append(ab_graph.number_of_edges())
#         ab_density.append(nx.density(ab_graph))
#
#         # Small-World Graph
#         sw_graph = nx.watts_strogatz_graph(n,k=9, p=0.1)
#         sw_edges.append(sw_graph.number_of_edges())
#         sw_density.append(nx.density(sw_graph))
#
#     print(sum(connect))
#     return er_edges, rgg_edges, ab_edges, sw_edges, er_density, rgg_density, ab_density, sw_density
# name="sdsdsds"
# def plot_distribution(edges, title):
#     plt.hist(edges, bins=20, alpha=0.5)
#     plt.title(title)
#     plt.xlabel('Number of Edges')
#     plt.ylabel('Frequency')
#     plt.legend(['ER;p=0.2', 'RGG;d=0.3', 'Albert-Barabasi;m=5', 'Small World;-k=9, p=0.1'])
#     plt.savefig(f'edges.png',format="png")
#     plt.show()

# def plot_density(density, title):
#     plt.hist(density, bins=20, alpha=0.5)
#     plt.title(title)
#     plt.xlabel('Density')
#     plt.ylabel('Frequency')
#     plt.legend(['ER', 'RGG', 'Albert-Barabasi', 'Small World'])
#     plt.show()

# # Parameters
# n = 50  # Number of nodes
# p_or_k = 0.1  # Probability for ER and RGG, and k (nearest neighbors) for Albert-Barabasi and Small World
#
# # Create networks
# er_edges, rgg_edges, ab_edges, sw_edges, er_density, rgg_density, ab_density, sw_density = create_networks(n, p_or_k)
#
# # Plot distribution of edges
# plot_distribution([er_edges, rgg_edges, ab_edges, sw_edges], 'Distribution of Edges')

# Plot density
# plot_density([er_density, rgg_density, ab_density, sw_density], 'Density Distribution')
#
# net = nx.newman_watts_strogatz_graph(50,10,0.1)
# nx.draw_networkx(net)
# plt.show()
# net = nx.random_geometric_graph(50,0.25)
# nx.draw_networkx(net)
# plt.show()
# net = nx.erdos_renyi_graph(50,0.15)
# nx.draw_networkx(net)
# plt.show()