from collections import Counter
import warnings
import numpy as np
from scope.operations.graphs import build_graph, get_signatures

########################
## Input Preparation  ##
########################
def _validate_matrix(matrix, natoms: int, debug: int=0):
    # Small helper to ensure symmetry and size.
    matrix = np.asarray(matrix)
    if debug > 1: print(f"VALIDATE_MATRIX: Inspecting matrix with shape={matrix.shape}; expected_shape=({natoms}, {natoms})")
    if matrix.shape != (natoms, natoms):
        raise ValueError(f"OVERLAP_MOLECULES: Inspected Matrix must have shape ({natoms}, {natoms})")
    if not np.allclose(matrix, matrix.T):
        raise ValueError(f"OVERLAP_MOLECULES: Inspected Matrix must be symmetric")
    if debug > 1: print("VALIDATE_MATRIX: Matrix shape and symmetry are valid")
    return matrix.copy()

def _get_center(labels, coords, center_method: str, debug: int=0):
    from scope.connectivity import compute_centroid
    from scope.other import get_metal_idxs

    if center_method == "centroid":
        center = np.asarray(compute_centroid(coords))
        if debug > 1: print(f"GET_CENTER: Using molecular centroid={center}")
        return center
    if center_method == "metal":
        metal_indices = get_metal_idxs(labels)
        if not metal_indices:
            raise ValueError("OVERLAP_MOLECULES: metal centering requires a metal atom")
        center = np.asarray(coords)[metal_indices[0]]
        if debug > 1: print(f"GET_CENTER: Using first metal atom at index={metal_indices[0]}, label={labels[metal_indices[0]]}, center={center}")
        return center
    raise ValueError(f"OVERLAP_MOLECULES: unknown center method: {center_method}")

def _get_connectivity(labels1, coords1, labels2, coords2, adjmat1, adjmat2, debug: int):
    from scope.connectivity import get_adjmatrix

    assert len(labels1) == len(labels2), "OVERLAP_MOLECULES: molecules must contain the same number of atoms"
    assert len(coords1) == len(labels1) and len(coords2) == len(labels2), "OVERLAP_MOLECULES: coordinates and labels must match their array size"
    natoms = len(labels1)
    if (adjmat1 is None) != (adjmat2 is None):
        raise ValueError("OVERLAP_MOLECULES: adjmat1 and adjmat2 must be supplied together")

    # Coordinate-only callers retain the previous connectivity fallback.
    if adjmat1 is None:
        if debug > 0: print("GET_CONNECTIVITY: Computing adjacency matrices from atomic labels and coordinates")
        isgood1, adjmat1, adjnum1 = get_adjmatrix(labels1, coords1, smart=True, debug=debug)
        isgood2, adjmat2, adjnum2 = get_adjmatrix(labels2, coords2, smart=True, debug=debug)
    else:
        if debug > 0: print("GET_CONNECTIVITY: Using the supplied adjacency matrices")
        adjmat1 = _validate_matrix(adjmat1, natoms, debug=debug)
        adjmat2 = _validate_matrix(adjmat2, natoms, debug=debug)
        adjnum1 = np.sum(adjmat1 > 0, axis=1)
        adjnum2 = np.sum(adjmat2 > 0, axis=1)
        isgood1 = True
        isgood2 = True

    if not isgood1 or np.any(adjnum1 == 0):
        raise ValueError("OVERLAP_MOLECULES: first molecule might be fragmented")
    if not isgood2 or np.any(adjnum2 == 0):
        raise ValueError("OVERLAP_MOLECULES: second molecule might be fragmented")
    if debug > 0: print(f"GET_CONNECTIVITY: Connectivity is valid for both molecules; edges1={int(np.sum(adjmat1 > 0)//2)}, edges2={int(np.sum(adjmat2 > 0)//2)}")
    if debug > 1: print(f"GET_CONNECTIVITY: {adjnum1=}; {adjnum2=}")
    return adjmat1, adjnum1, adjmat2, adjnum2

def _get_bond_orders(bond_orders1, bond_orders2, natoms: int, debug: int=0):
    if (bond_orders1 is None) != (bond_orders2 is None):
        raise ValueError("OVERLAP_MOLECULES: bond_orders1 and bond_orders2 must be supplied together")
    if bond_orders1 is None:
        if debug > 0: print("GET_BOND_ORDERS: Bond-order matrices were not supplied; graph matching will use connectivity only")
        return None, None
    bond_orders1 = _validate_matrix(bond_orders1, natoms, debug=debug)
    bond_orders2 = _validate_matrix(bond_orders2, natoms, debug=debug)
    if debug > 0: print("GET_BOND_ORDERS: Bond-order matrices are valid and will constrain atom mappings")
    if debug > 1: print(f"GET_BOND_ORDERS: orders1={sorted(set(bond_orders1.flatten()))}; orders2={sorted(set(bond_orders2.flatten()))}")
    return bond_orders1, bond_orders2

def get_extended_info(labels, adjmat, adjnum, debug: int=0):
    """Build atom labels decorated with their topological environments."""
    graph             = build_graph(adjmat, label=labels)
    _, _, signatures = get_signatures(graph)

    # get_signatures() is organized by layer and then by atom. Transpose it
    # into one stable signature string per atom for the Hungarian groups.
    layers          = sorted(signatures)
    atom_signatures = []
    for atom_index in range(len(labels)):
        parts = ["".join(signatures[layer][atom_index]) for layer in layers]
        atom_signatures.append("_".join(parts))

    information = np.asarray([f"{label}{neighbors}{signature}" for label, neighbors, signature in zip(labels, adjnum, atom_signatures)])
    if debug > 0: print(f"GET_EXTENDED_INFO: Built topological information for {len(information)} atoms across {len(layers)} signature layers")
    if debug > 1: print(f"GET_EXTENDED_INFO: {information=}")
    return information

def _get_atom_information(labels1, adjmat1, adjnum1, labels2, adjmat2, adjnum2, use_ext_info: bool, debug: int):
    if debug > 0: print(f"GET_ATOM_INFORMATION: Using {'topological signatures' if use_ext_info else 'element labels'} to classify atoms")
    if use_ext_info:
        data1 = get_extended_info(labels1, adjmat1, adjnum1, debug=debug)
        data2 = get_extended_info(labels2, adjmat2, adjnum2, debug=debug)
    else:
        data1 = np.asarray(labels1, dtype=str)
        data2 = np.asarray(labels2, dtype=str)

    # Counters preserve multiplicity; sets would miss repeated atom classes.
    if Counter(data1) != Counter(data2):
        if debug > 0: print(f"GET_ATOM_INFORMATION: Atom classes differ; only_in_molecule1={Counter(data1) - Counter(data2)}, only_in_molecule2={Counter(data2) - Counter(data1)}")
        raise ValueError("OVERLAP_MOLECULES: the molecules have different topological atom information")
    if debug > 0: print(f"GET_ATOM_INFORMATION: Atom information is compatible; unique_classes={len(Counter(data1))}")
    return data1, data2


############################
## Graph Mapping Search   ##
############################
def _get_terminal_atoms(graph, center_method: str, debug: int=0):
    from scope.other import get_metal_idxs

    metal_indices = set()
    if center_method == "metal":
        labels = [graph.nodes[index]["label"] for index in graph.nodes]
        metal_indices = set(get_metal_idxs(labels))

    terminal_atoms = set()
    for atom_index in graph.nodes:
        if atom_index in metal_indices:
            continue
        neighbors = list(graph.neighbors(atom_index))
        if len(neighbors) == 1 and graph.degree[neighbors[0]] > 1:
            terminal_atoms.add(atom_index)
    if debug > 1: print(f"GET_TERMINAL_ATOMS: Found {len(terminal_atoms)} terminal atoms; indices={sorted(terminal_atoms)}")
    return terminal_atoms


def _get_terminal_key(graph, parent_index: int, atom_index: int, use_bond_orders: bool):
    key = graph.nodes[atom_index]["label"]
    if use_bond_orders:
        key = (key, graph.edges[parent_index, atom_index]["order"])
    return key

def _get_core_graph(graph, terminal_atoms, use_bond_orders: bool, debug: int=0):
    core_indices = [index for index in graph.nodes if index not in terminal_atoms]
    core_graph   = graph.subgraph(core_indices).copy()

    # Store the removed terminal environments on their parent core atoms.
    for parent_index in core_graph.nodes:
        terminal_information = [_get_terminal_key(graph, parent_index, neighbor, use_bond_orders) for neighbor in graph.neighbors(parent_index) if neighbor in terminal_atoms]
        core_graph.nodes[parent_index]["terminal_atoms"] = tuple(sorted(terminal_information))
    if debug > 1: print(f"GET_CORE_GRAPH: Reduced graph from {graph.number_of_nodes()} to {core_graph.number_of_nodes()} nodes while retaining terminal-atom environments")
    return core_graph

def _group_terminal_atoms(graph, parent_index: int, terminal_atoms, use_bond_orders: bool):
    groups = {}
    for atom_index in terminal_atoms:
        key = _get_terminal_key(graph, parent_index, atom_index, use_bond_orders)
        groups.setdefault(key, []).append(atom_index)
    return groups

def _complete_core_mapping(core_mapping, graph1, graph2, terminal1, terminal2, use_bond_orders: bool, debug: int=0):
    indices         = np.full(graph1.number_of_nodes(), -1, dtype=int)
    terminal_groups = []

    for reference_index, mobile_index in core_mapping.items():
        indices[reference_index] = mobile_index
        reference_terminal = sorted(set(graph1.neighbors(reference_index)) & terminal1)
        mobile_terminal    = sorted(set(graph2.neighbors(mobile_index)) & terminal2)
        reference_groups = _group_terminal_atoms(graph1, reference_index, reference_terminal, use_bond_orders)
        mobile_groups    = _group_terminal_atoms(graph2, mobile_index, mobile_terminal, use_bond_orders)

        if set(reference_groups) != set(mobile_groups):
            raise ValueError("OVERLAP_MOLECULES: inconsistent terminal-atom groups")
        for key, reference_group in reference_groups.items():
            mobile_group = mobile_groups[key]
            if len(reference_group) != len(mobile_group):
                raise ValueError("OVERLAP_MOLECULES: inconsistent terminal-atom groups")
            indices[reference_group] = mobile_group
            terminal_groups.append((reference_group, mobile_group))

    if np.any(indices < 0):
        raise ValueError("OVERLAP_MOLECULES: graph mapping did not include every atom")
    if debug > 2: print(f"COMPLETE_CORE_MAPPING: Completed core mapping with {len(terminal_groups)} terminal groups; indices={indices.tolist()}")
    return indices, terminal_groups

def _assign_terminal_groups(indices, terminal_groups, rotation, translation, coords1, coords2, debug: int=0):
    from scope.reconstruct import hungarian
    for reference_group, mobile_group in terminal_groups:
        mobile_coords            = np.asarray(coords2)[mobile_group]
        aligned_mobile_coords    = (rotation @ mobile_coords.T).T + translation
        assignment               = hungarian(np.asarray(coords1)[reference_group], aligned_mobile_coords, debug=debug)
        indices[reference_group] = np.asarray(mobile_group)[assignment]
        if debug > 2: print(f"ASSIGN_TERMINAL_GROUPS: reference_group={reference_group}, mobile_group={mobile_group}, assignment={assignment.tolist()}")

def _optimize_terminal_groups(labels1, coords1, labels2, coords2, core_mapping, indices, terminal_groups, center_method: str, max_iter: int, debug: int=0):
    if terminal_groups:
        # Start from the core alignment so arbitrary terminal ordering has less
        # influence on the first complete-molecule rotation.
        reference_core                  = np.asarray(list(core_mapping.keys()))
        mobile_core                     = np.asarray([core_mapping[index] for index in reference_core])
        rotation, translation, _, _     = kabsch_align(np.asarray(labels2)[mobile_core], np.asarray(coords2)[mobile_core], np.asarray(labels1)[reference_core], np.asarray(coords1)[reference_core], center_method=center_method, debug=debug)
        _assign_terminal_groups(indices, terminal_groups, rotation, translation, coords1, coords2, debug=debug)
        if debug > 2: print(f"OPTIMIZE_TERMINAL_GROUPS: Initialized {len(terminal_groups)} terminal groups from the core alignment")

    # Alternate complete-molecule alignment and assignment inside each fixed
    # chemical terminal group until its permutation no longer changes.
    previous_indices = None
    iterations_run   = 0
    for iteration in range(max_iter):
        iterations_run               = iteration + 1
        rotation, translation, _, _  = kabsch_align(np.asarray(labels2)[indices], np.asarray(coords2)[indices], labels1, coords1, center_method=center_method, debug=debug)
        _assign_terminal_groups(indices, terminal_groups, rotation, translation, coords1, coords2, debug=debug)
        current_indices = tuple(indices)
        if current_indices == previous_indices:
            break
        previous_indices = current_indices

    _, _, aligned_coords, rmsd_value = kabsch_align(np.asarray(labels2)[indices], np.asarray(coords2)[indices], labels1, coords1, center_method=center_method, debug=debug)
    if debug > 2: print(f"OPTIMIZE_TERMINAL_GROUPS: Finished after {iterations_run} iterations with rmsd={rmsd_value:.6f}")
    return rmsd_value, indices, aligned_coords

def _search_graph_mappings(labels1, coords1, adjmat1, labels2, coords2, adjmat2, bond_orders1, bond_orders2, center_method: str, max_mappings: int, max_iter: int, debug: int=0):
    from networkx.algorithms.isomorphism import GraphMatcher, categorical_node_match, numerical_edge_match

    use_bond_orders = bond_orders1 is not None
    edge_features1  = None
    edge_features2  = None
    if use_bond_orders:
        edge_features1 = {"order": bond_orders1}
        edge_features2 = {"order": bond_orders2}
    if debug > 0: print(f"GRAPH_MAPPING_SEARCH: Starting graph search; max_mappings={max_mappings}, use_bond_orders={use_bond_orders}")
    graph1      = build_graph(adjmat1, edge_features=edge_features1, label=labels1)
    graph2      = build_graph(adjmat2, edge_features=edge_features2, label=labels2)
    terminal1   = _get_terminal_atoms(graph1, center_method, debug=debug)
    terminal2   = _get_terminal_atoms(graph2, center_method, debug=debug)
    core_graph1 = _get_core_graph(graph1, terminal1, use_bond_orders, debug=debug)
    core_graph2 = _get_core_graph(graph2, terminal2, use_bond_orders, debug=debug)

    node_match = categorical_node_match(["label", "terminal_atoms"], [None, ()])
    edge_match = numerical_edge_match("order", 1.0) if use_bond_orders else None

    matcher = GraphMatcher(core_graph1, core_graph2, node_match=node_match, edge_match=edge_match)
    if not matcher.is_isomorphic():
        if not use_bond_orders:
            if debug > 0: print("GRAPH_MAPPING_SEARCH: Core graphs are not isomorphic under the requested element and adjacency constraints")
            raise ValueError("OVERLAP_MOLECULES: molecular graphs are not isomorphic")

        if debug > 0: print("GRAPH_MAPPING_SEARCH: Strict bond-order mapping failed; retrying with element and adjacency constraints")
        # Resonance-equivalent structures can localize formal bond orders on
        # different atoms while retaining the same element-labelled adjacency.
        use_bond_orders = False
        core_graph1     = _get_core_graph(graph1, terminal1, use_bond_orders, debug=debug)
        core_graph2     = _get_core_graph(graph2, terminal2, use_bond_orders, debug=debug)
        edge_match      = None
        matcher         = GraphMatcher(core_graph1, core_graph2, node_match=node_match, edge_match=edge_match)
        if not matcher.is_isomorphic():
            if debug > 0: print("GRAPH_MAPPING_SEARCH: Core graphs are not isomorphic under the requested element and adjacency constraints")
            raise ValueError("OVERLAP_MOLECULES: molecular graphs are not isomorphic")
        warnings.warn("OVERLAP_MOLECULES: Molecular connectivity matches, but the stored bond orders cannot be mapped exactly. Continuing without bond-order constraints.", RuntimeWarning, stacklevel=2)
        if debug > 0: print("GRAPH_MAPPING_SEARCH: Relaxed element and adjacency mapping succeeded")

    best_rmsd        = np.inf
    best_indices     = None
    best_coords      = None
    mappings_checked = 0
    if debug > 0: print(f"GRAPH_MAPPING_SEARCH: Core graphs are isomorphic; core_nodes1={core_graph1.number_of_nodes()}, core_nodes2={core_graph2.number_of_nodes()}, terminal_atoms1={len(terminal1)}, terminal_atoms2={len(terminal2)}")

    # Evaluate chemically valid core mappings. Equivalent terminal atoms are
    # optimized separately, avoiding their factorial graph permutations.
    for core_mapping in matcher.isomorphisms_iter():
        indices, terminal_groups           = _complete_core_mapping(core_mapping, graph1, graph2, terminal1, terminal2, use_bond_orders, debug=debug)
        rmsd_value, indices, aligned_coords = _optimize_terminal_groups(labels1, coords1, labels2, coords2, core_mapping, indices, terminal_groups, center_method, max_iter, debug=debug)
        if rmsd_value < best_rmsd:
            best_rmsd    = rmsd_value
            best_indices = indices.copy()
            best_coords  = aligned_coords.copy()
            if debug > 1: print(f"GRAPH_MAPPING_SEARCH: New best mapping at mapping={mappings_checked + 1}; rmsd={best_rmsd:.6f}")

        mappings_checked += 1
        if mappings_checked >= max_mappings:
            if debug > 0:
                print(f"GRAPH_MAPPING_SEARCH: Reached max_mappings={max_mappings}")
            break
    if debug > 0: print(f"GRAPH_MAPPING_SEARCH: Finished after checking {mappings_checked} mappings; best_rmsd={best_rmsd:.6f}")
    return best_rmsd, best_indices, best_coords, mappings_checked


def _mapping_is_valid(adjmat1, adjmat2, bond_orders1, bond_orders2, debug: int=0):
    if not np.array_equal(adjmat1 > 0, adjmat2 > 0):
        if debug > 1: print("MAPPING_IS_VALID: Mapping rejected because adjacency matrices differ")
        return False
    if bond_orders1 is not None and not np.allclose(bond_orders1, bond_orders2):
        if debug > 1: print("MAPPING_IS_VALID: Mapping rejected because bond-order matrices differ")
        return False
    if debug > 1: print("MAPPING_IS_VALID: Mapping preserves connectivity and bond orders")
    return True


############################
## Fast Hungarian Search  ##
############################
def _run_hungarian_search(labels1, coords1, data1, adjmat1, bond_orders1, labels2, coords2, data2, adjmat2, bond_orders2, center_method: str, max_iter: int, rmsd_convergence: float, debug: int=0):
    from scope.reconstruct import reorder_hungarian

    natoms               = len(labels1)
    labels2_current      = labels2.copy()
    coords2_current      = coords2.copy()
    data2_current        = data2.copy()
    adjmat2_current      = adjmat2.copy()
    bond_orders2_current = None if bond_orders2 is None else bond_orders2.copy()
    indices2_current     = np.arange(natoms)

    best_rmsd         = np.inf
    best_indices      = None
    best_coords       = None
    valid_rmsd_values = []
    seen_mappings     = set()
    iterations_run    = 0
    if debug > 0: print(f"HUNGARIAN_SEARCH: Starting iterative alignment and assignment; max_iter={max_iter}, rmsd_convergence={rmsd_convergence}")

    # Fast path: Kabsch alignment with Hungarian atom assignment.
    for iteration in range(max_iter):
        iterations_run              = iteration + 1
        _, _, coords2_aligned, _    = kabsch_align(labels2_current, coords2_current, labels1, coords1, center_method=center_method, debug=debug)
        step_mapping                = reorder_hungarian(data1, data2_current, coords1, coords2_aligned, debug=debug)

        labels2_reordered  = labels2_current[step_mapping]
        coords2_reordered  = coords2_aligned[step_mapping]
        data2_reordered    = data2_current[step_mapping]
        adjmat2_reordered  = adjmat2_current[np.ix_(step_mapping, step_mapping)]
        indices2_reordered = indices2_current[step_mapping]
        if bond_orders2_current is None:
            bond_orders2_reordered = None
        else:
            bond_orders2_reordered = bond_orders2_current[np.ix_(step_mapping, step_mapping)]

        # Only retain mappings that preserve connectivity and formal bond order.
        candidate_rmsd = None
        graph_is_valid = _mapping_is_valid(adjmat1, adjmat2_reordered, bond_orders1, bond_orders2_reordered, debug=debug)
        if graph_is_valid:
            _, _, candidate_coords, candidate_rmsd = kabsch_align(labels2_reordered, coords2_reordered, labels1, coords1, center_method=center_method, debug=debug)
            valid_rmsd_values.append(candidate_rmsd)
            if candidate_rmsd < best_rmsd:
                best_rmsd    = candidate_rmsd
                best_indices = indices2_reordered.copy()
                best_coords  = candidate_coords.copy()
        candidate_text = f"{candidate_rmsd:.6f}" if candidate_rmsd is not None else "not available"
        if debug > 0: print(f"HUNGARIAN_SEARCH: iteration={iterations_run}, mapping_valid={graph_is_valid}, candidate_rmsd={candidate_text}")
        if debug > 1: print(f"HUNGARIAN_SEARCH: step_mapping={step_mapping.tolist()}, composed_mapping={indices2_reordered.tolist()}")

        # Stop on a stable local assignment, a cycle, or converged RMSD.
        composed_mapping = tuple(indices2_reordered)
        mapping_repeated = composed_mapping in seen_mappings
        no_reordering    = np.array_equal(step_mapping, np.arange(natoms))
        rmsd_converged   = len(valid_rmsd_values) >= 2 and abs(valid_rmsd_values[-2] - valid_rmsd_values[-1]) <= rmsd_convergence
        seen_mappings.add(composed_mapping)
        stop_reasons = []
        if no_reordering:    stop_reasons.append("assignment did not change")
        if mapping_repeated: stop_reasons.append("mapping cycle detected")
        if rmsd_converged:   stop_reasons.append("RMSD converged")
        if stop_reasons:
            if debug > 0: print(f"HUNGARIAN_SEARCH: Stopping because {', '.join(stop_reasons)}")
            break

        labels2_current      = labels2_reordered
        coords2_current      = coords2_reordered
        data2_current        = data2_reordered
        adjmat2_current      = adjmat2_reordered
        bond_orders2_current = bond_orders2_reordered
        indices2_current     = indices2_reordered

    if debug > 0 and best_indices is None: print(f"HUNGARIAN_SEARCH: Finished after {iterations_run} iterations without a chemically valid mapping")
    if debug > 0 and best_indices is not None: print(f"HUNGARIAN_SEARCH: Finished after {iterations_run} iterations; best_rmsd={best_rmsd:.6f}")
    return best_rmsd, best_indices, best_coords, iterations_run


######################
## Public Interface ##
######################
def overlap_molecules(labels1, coords1, labels2, coords2, center_method: str="centroid", use_ext_info: bool=True, translate_to_ref: bool=True, max_iter: int=5, adjmat1=None, adjmat2=None, bond_orders1=None, bond_orders2=None, max_graph_mappings: int=1000, rmsd_convergence: float=1e-8, debug: int=0):
    """Overlap equivalent molecules using topology-aware atom mappings.

    Parameters:
        adjmat1, adjmat2:               Optional stored adjacency matrices.
        bond_orders1, bond_orders2:     Optional stored bond-order matrices.
        max_graph_mappings:             Maximum graph-valid core mappings to test.
        debug:                          0 is silent; 1 reports decisions; 2 reports detailed mappings and coordinates; 3 reports low-level optimization details.

    Returns:
        tuple:                          Aligned labels, coordinates, and atom mapping.
    """
    from scope.read_write import print_xyz, write_xyz

    if len(labels1) != len(labels2):
        raise ValueError("OVERLAP_MOLECULES: molecules must contain the same number of atoms")
    if max_iter < 1:
        raise ValueError("OVERLAP_MOLECULES: max_iter must be at least 1")
    if max_graph_mappings < 0:
        raise ValueError("OVERLAP_MOLECULES: max_graph_mappings cannot be negative")

    labels1 = np.asarray(labels1)
    labels2 = np.asarray(labels2)
    coords1 = np.asarray(coords1)
    coords2 = np.asarray(coords2)
    natoms  = len(labels1)
    if debug > 0: print(f"OVERLAP_MOLECULES: Starting overlap; natoms={natoms}, center_method={center_method}, use_ext_info={use_ext_info}, max_iter={max_iter}, max_graph_mappings={max_graph_mappings}")

    # 1) Obtain and validate the chemical information used for matching.
    adjmat1, adjnum1, adjmat2, adjnum2 = _get_connectivity(labels1, coords1, labels2, coords2, adjmat1, adjmat2, debug=debug)
    bond_orders1, bond_orders2         = _get_bond_orders(bond_orders1, bond_orders2, natoms, debug=debug)
    data1, data2                       = _get_atom_information(labels1, adjmat1, adjnum1, labels2, adjmat2, adjnum2, use_ext_info, debug=debug)

    # 2) Center both geometries before assignment and alignment.
    center1         = _get_center(labels1, coords1, center_method, debug=debug)
    center2         = _get_center(labels2, coords2, center_method, debug=debug)
    coords1_centered = coords1 - center1
    coords2_centered = coords2 - center2

    # 3) Run the inexpensive iterative search and retain its best valid result.
    best_rmsd, best_indices, best_coords, iterations_run = _run_hungarian_search(labels1, coords1_centered, data1, adjmat1, bond_orders1, labels2, coords2_centered, data2, adjmat2, bond_orders2, center_method, max_iter, rmsd_convergence, debug=debug)
    best_method = "Hungarian"

    # 4) Search alternative graph-valid mappings. If the fast path found no
    # valid result, one graph mapping is still required as a repair step.
    mappings_to_check = max_graph_mappings
    if mappings_to_check == 0 and best_indices is None:
        mappings_to_check = 1
    if mappings_to_check > 0:
        graph_rmsd, graph_indices, graph_coords, mappings_checked = _search_graph_mappings(labels1, coords1_centered, adjmat1, labels2, coords2_centered, adjmat2, bond_orders1, bond_orders2, center_method, mappings_to_check, max_iter, debug=debug)
        if graph_rmsd < best_rmsd:
            best_rmsd    = graph_rmsd
            best_indices = graph_indices
            best_coords  = graph_coords
            best_method  = "Graph"
        if debug > 0: print(f"OVERLAP_MOLECULES: Graph search checked {mappings_checked} mappings; graph_rmsd={graph_rmsd:.6f}")

    if best_indices is None:
        raise ValueError("OVERLAP_MOLECULES: no chemically valid atom mapping was found")

    # 5) Restore the reference position and return the original mobile indices.
    labels2_reordered = labels2[best_indices]
    if translate_to_ref:
        coords1_final = coords1_centered + center1
        coords2_final = best_coords + center1
    else:
        coords1_final = coords1_centered
        coords2_final = best_coords

    if debug > 0: print(f"OVERLAP_MOLECULES: Finished successfully with method={best_method}, final_rmsd={best_rmsd:.6f}, iterative_passes={iterations_run}")
    if debug > 1:
        print("OVERLAP_MOLECULES: Final aligned coordinates")
        print_xyz(labels2_reordered, coords2_final)
    if debug > 1:
        output_folder = "./test_overlap/"
        import os
        os.makedirs(output_folder, exist_ok=True)
        write_xyz(output_folder + "final1.xyz", labels1, coords1_final)
        write_xyz(output_folder + "final2.xyz", labels2_reordered, coords2_final)

    return True, labels1, coords1_final, labels2_reordered, coords2_final, best_indices.tolist()


def rmsd(labels1, coords1, labels2, coords2, reorder: bool=False, center_method="centroid", atom_idxs: list=None, use_ext_info: bool=True, translate_to_ref: bool=True, max_iter: int=5, adjmat1=None, adjmat2=None, bond_orders1=None, bond_orders2=None, max_graph_mappings: int=1000, debug: int=0):
    """Compute RMSD, optionally reordering equivalent atoms before alignment."""
    if len(labels1) != len(labels2):
        raise ValueError("RMSD: molecules must contain the same number of atoms")
    coords1 = np.asarray(coords1)
    coords2 = np.asarray(coords2)
    if debug > 0: print(f"RMSD: Starting calculation; natoms={len(labels1)}, reorder={reorder}, center_method={center_method}")
    if atom_idxs is not None:
        atom_idxs = np.sort(np.asarray(atom_idxs))

    if reorder:
        isgood, labels1, coords1, labels2, coords2, _ = overlap_molecules(labels1, coords1, labels2, coords2, center_method=center_method, use_ext_info=use_ext_info, translate_to_ref=translate_to_ref, max_iter=max_iter,adjmat1=adjmat1, adjmat2=adjmat2, bond_orders1=bond_orders1, bond_orders2=bond_orders2, max_graph_mappings=max_graph_mappings, debug=debug)
        if not isgood:
            return None
        if not np.array_equal(labels1, labels2):
            raise ValueError("RMSD: reordered atom labels do not coincide")

    if atom_idxs is not None:
        coords1 = coords1[atom_idxs]
        coords2 = coords2[atom_idxs]
    value = np.sqrt(np.mean(np.sum((coords1 - coords2) ** 2, axis=1)))
    if debug > 0: print(f"RMSD: Finished calculation; compared_atoms={len(coords1)}, rmsd={value:.6f}")
    return np.round(value, 4)


######################
## Kabsch Alignment ##
######################
def kabsch_align(labels1, coords1, labels2, coords2, center_method: str="centroid", debug: int=0):
    """Align the first coordinate set onto the second with the Kabsch algorithm."""
    coords1          = np.asarray(coords1)
    coords2          = np.asarray(coords2)
    center1          = _get_center(labels1, coords1, center_method, debug=debug)
    center2          = _get_center(labels2, coords2, center_method, debug=debug)
    coords1_centered = coords1 - center1
    coords2_centered = coords2 - center2

    covariance                = coords1_centered.T @ coords2_centered
    left, _, right_transposed = np.linalg.svd(covariance)
    rotation                  = right_transposed.T @ left.T
    reflection_corrected      = np.linalg.det(rotation) < 0
    if reflection_corrected:
        right_transposed[2, :] *= -1
        rotation = right_transposed.T @ left.T

    translation   = center2 - rotation @ center1
    aligned_coords = (rotation @ coords1.T).T + translation
    rmsd_value     = np.sqrt(np.mean(np.sum((aligned_coords - coords2) ** 2, axis=1)))
    if debug > 1: print(f"KABSCH_ALIGN: Aligned {len(coords1)} atoms with center_method={center_method}; rmsd={rmsd_value:.6f}, reflection_corrected={reflection_corrected}")
    return rotation, translation, aligned_coords, rmsd_value


def kabsch_rotate(coords1, coords2, debug: int=0):
    """Return the Kabsch rotation that aligns coords2 onto coords1."""
    covariance       = np.dot(np.asarray(coords2).T, np.asarray(coords1))
    left, _, right   = np.linalg.svd(covariance)
    determinant_sign = np.sign(np.linalg.det(left) * np.linalg.det(right))
    correction       = np.diag([1, 1, determinant_sign])
    rotation         = np.dot(np.dot(left, correction), right)
    if debug > 0: print(f"KABSCH_ROTATE: Built rotation for {len(coords1)} atoms; determinant={np.linalg.det(rotation):.6f}")
    return rotation
