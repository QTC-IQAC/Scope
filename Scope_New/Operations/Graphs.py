####
def build_graph(adj_matrix, labels):
    import networkx as nx
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
    import networkx as nx
    from networkx.algorithms.isomorphism import GraphMatcher
    matcher = GraphMatcher(G1, G2, node_match=lambda n1, n2: (n1['label'] == n2['label']))
    if matcher.is_isomorphic():
        mapping = matcher.mapping
        return [mapping[i] for i in range(len(mapping))]
    else:
        return None
