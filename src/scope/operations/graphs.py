import networkx as nx
import numpy as np

####
def build_graph(adj_matrix, labels):
    G = nx.Graph()
    N = len(labels)
    for i in range(N):
        G.add_node(i, label=labels[i])#, feature=features[i] if features else None)
    for i in range(N):
        for j in range(i+1, N):
            if adj_matrix[i, j] > 0:
                G.add_edge(i, j)
    return G

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
def compare_graphs(G1, G2):

    # Compares whether two graphs are identical using resistances
    # 1) compares number of nodes and edges
    if G1.number_of_nodes() != G2.number_of_nodes(): return False
    if G1.number_of_edges() != G2.number_of_edges(): return False

    # 2) compares topologies via resistance matrices:
    # 2.1) Compute resistance matrices
    R1 = resistance_matrix(G1)
    R2 = resistance_matrix(G2)

    # 2.2) Compare their eigenvalue spectra. If spectra differ, graphs are not isomorphic
    eig1 = np.sort(np.linalg.eigvalsh(R1))
    eig2 = np.sort(np.linalg.eigvalsh(R2))
    if not np.allclose(eig1, eig2, atol=1e-6):       return False

    # 2.3) If spectra match, perform a graph isomorphism test
    GM = nx.isomorphism.GraphMatcher(G1, G2)
    if not GM.is_isomorphic(): return False

    # 3) Compares Topologies via signatures
    _, n_neigh_x_layer1, sign1 = get_signatures(G1) 
    _, n_neigh_x_layer2, sign2 = get_signatures(G2) 
    # 3.1) Get neighbours per layer. If they differ, graphs are not isomorphic
    if n_neigh_x_layer1 != n_neigh_x_layer2: return False

    # 3.2) Compares signatures:
    if not compare_signatures(sign1, sign2): return False

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