import networkx as nx
import numpy as np

####
import networkx as nx

def build_graph(adj_matrix, debug: int = 0, **node_features):
    """
    Builds a NetworkX graph from an adjacency matrix and node features provided as kwargs.
    It might result in a different graph than using the Specie-class function

    Parameters
    ----------
    adj_matrix : (N, N) array
        Adjacency matrix
    node_features : keyword arguments
        Each value must be list-like of the same length as the adj_matrix.
        Example: labels=[...], mass=[...]

    Returns
    -------
    G : networkx.Graph
    """
    G = nx.Graph()
    N = len(adj_matrix)

    # Reads features and performs sanity checks
    for name, values in node_features.items():
        if len(values) != N: raise ValueError(f"Feature '{name}' has length {len(values)}, expected {N}")

    # Add nodes with attributes
    for i in range(N):
        attrs = {}
        for name, values in node_features.items():
            attrs[name] = values[i]
        G.add_node(i, **attrs)
        if debug > 0:
            print(f"BUILD_GRAPH: node {i} created with attrs {attrs}")

    # Add edges
    for i in range(N):
        for j in range(i + 1, N):
            if adj_matrix[i, j] > 0:
                G.add_edge(i, j)
                if debug > 0:
                    print(f"BUILD_GRAPH: edge {i}-{j} created from adjacency matrix")

    if debug > 0:
        print(f"BUILD_GRAPH: {G.number_of_nodes()} nodes created")
        print(f"BUILD_GRAPH: {G.number_of_edges()} edges created")
    return G

#def build_graph(adj_matrix, labels, debug: int=0):
#    # Builds a NetworkX graph from an adjacency matrix and node labels
#    # It might result in a different graph than using the Specie-class function
#    G = nx.Graph()
#    N = len(labels)
#    for i in range(N):
#        G.add_node(i, label=labels[i])#, feature=features[i] if features else None)
#        if debug > 0: print(f"BUILD_GRAPH: node created for {i}:{labels[i]}") 
#    for i in range(N):
#        for j in range(i+1, N):
#            if adj_matrix[i, j] > 0:
#                G.add_edge(i, j)
#                if debug > 0: print(f"BUILD_GRAPH: edge created between {i}:{labels[i]} and {j}:{labels[j]}")
#    if debug > 0:
#        print(f"BUILD_GRAPH: {G.number_of_nodes()} nodes created")
#        print(f"BUILD_GRAPH: {G.number_of_edges()} edges created")
#    return G

####
def get_permutation_from_isomorphism(G1, G2):
    from networkx.algorithms.isomorphism import GraphMatcher
    matcher = GraphMatcher(G1, G2, node_match=lambda n1, n2: (n1['label'] == n2['label']))
    if matcher.is_isomorphic():
        mapping = matcher.mapping
        return [mapping[i] for i in range(len(mapping))]
    else:
        return None

####
def resistance_matrix(G):
    """Computes the resistance-distance matrix for a graph."""
    L = nx.laplacian_matrix(G).astype(float).toarray()
    L_pinv = np.linalg.pinv(L)  # Moore–Penrose pseudoinverse
    n = L.shape[0]
    R = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            R[i, j] = L_pinv[i, i] + L_pinv[j, j] - 2 * L_pinv[i, j]
    return R

####
def compare_graphs(G1, G2, debug: int=0):

    # Compares whether two graphs are identical using resistances
    # 1) compares number of nodes and edges
    if G1.number_of_nodes() != G2.number_of_nodes(): 
        if debug > 0: print(f"COMPARE_GRAPHS: different number of nodes: {G1.number_of_nodes()} vs. {G2.number_of_nodes()}")
        return False
    if G1.number_of_edges() != G2.number_of_edges(): 
        if debug > 0: print(f"COMPARE_GRAPHS: different number of edges: {G1.number_of_edges()} vs. {G2.number_of_edges()}")
        return False

    # 2) compares topologies via resistance matrices:
    # 2.1) Compute resistance matrices
    R1 = resistance_matrix(G1)
    R2 = resistance_matrix(G2)

    # 2.2) Compare their eigenvalue spectra. If spectra differ, graphs are not isomorphic
    eig1 = np.sort(np.linalg.eigvalsh(R1))
    eig2 = np.sort(np.linalg.eigvalsh(R2))
    if not np.allclose(eig1, eig2, atol=1e-6):       
        if debug > 0: print(f"COMPARE_GRAPHS: different eigenvalue spectra for the resistance matrices")
        if debug > 0: print(f"\t {eig1}")
        if debug > 0: print(f"\t {eig2}")
        return False

    # 2.3) If spectra match, perform a graph isomorphism test
    GM = nx.isomorphism.GraphMatcher(G1, G2)
    if not GM.is_isomorphic(): 
        if debug > 0: print(f"COMPARE_GRAPHS: graphs are not isomorphic according to NX test")
        return False

    # 3) Compares Topologies via signatures
    _, n_neigh_x_layer1, sign1 = get_signatures(G1) 
    _, n_neigh_x_layer2, sign2 = get_signatures(G2) 
    # 3.1) Get neighbours per layer. If they differ, graphs are not isomorphic
    if n_neigh_x_layer1 != n_neigh_x_layer2: 
        if debug > 0: print(f"COMPARE_GRAPHS: different neighbours per layer")
        if debug > 0: print(f"\t {n_neigh_x_layer1}")
        if debug > 0: print(f"\t {n_neigh_x_layer2}")
        return False

    # 3.2) Compares signatures:
    if not compare_signatures(sign1, sign2): 
        if debug > 0: print(f"COMPARE_GRAPHS: different signatures")
        if debug > 0: print(f"\t {sign1}")
        if debug > 0: print(f"\t {sign2}")
        return False

    return True

####
def nth_neighbor_labels(G, source, n):
    # Returns a list including the label of the source, and the labels of all nodes at exactly distance n from source
    layers = {source: 0}
    queue = [source]

    while queue:
        curr = queue.pop(0)
        for neigh in G.neighbors(curr):
            if neigh not in layers:
                layers[neigh] = layers[curr] + 1
                queue.append(neigh)
    labels = []
    labels.append(G.nodes[source].get("label", None))
    for node, dist in layers.items():
        if dist == n:
            labels.append(G.nodes[node].get("label", None))
    return tuple(sorted(str(l) for l in labels))

####
def get_signatures(G, convergence_layers: int=2):
    n_neigh_x_layer = []
    signatures = {}

    n = 1
    while True:
        new_signatures = {}

        # Compute signatures at distance n for each node
        for node in G.nodes:
            new_signatures[node] = nth_neighbor_labels(G, node, n)

        # Count unique signatures
        Cn = len(set(new_signatures.values()))

        # Evaluates if any is empty
        has_empty = any(len(idict) == 0 for idict in new_signatures.values())

        # Checks Convergence
        if n_neigh_x_layer.count(Cn) >= convergence_layers or has_empty:
            break

        # Applies updates
        n_neigh_x_layer.append(Cn)
        signatures[n] = new_signatures
        n += 1

    # While is closed, get the results
    return n, n_neigh_x_layer, signatures

####
def compare_signatures(sign1: dict, sign2: dict):
    from collections import Counter
    """
    Checks if two specific types of dictionaries, the signatures obtained in scope.operations.graphs.get_signatures() are equivalent:
    """
    # 1) Compare that they have the same layers
    if not list(sign1.keys()) == list(sign2.keys()): return False
    # 2) Compares the counter
    for layer in list(sign1.keys()):
        if not Counter(sign1[layer].values()) == Counter(sign2[layer].values()): return False
    return True

####
def compute_topological_distances(G, ref_atom: int) -> dict:
    # Compute shortest path lengths from ref_atom
    dist_dict = nx.single_source_shortest_path_length(G, ref_atom)
    distances = np.full(G.number_of_nodes(), -1, dtype=int)
    for idx, d in dist_dict.items():
        distances[idx] = d
    return distances

####
def print_graph_info(G):
    print(f"Nodes: {G.number_of_nodes()}")
    print(f"Edges: {G.number_of_edges()}")
