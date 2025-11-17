import warnings
import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

from typing import Tuple
from scope.operations.dicts_and_lists import extract_from_list  
from scope.elementdata import ElementData  
elemdatabase = ElementData()

######################################
#### Functions Related to Species ####
######################################
def labels2formula(labels: list):
    elems = elemdatabase.elementnr.keys()
    formula=[]
    for z in elems:
        nz = labels.count(z)
        if nz > 1:
            formula.append(f"{z}{nz}-")
        if nz == 1:
            formula.append(f"{z}-")
    formula = ''.join(formula)[:-1]
    return formula

######
def labels2ratio(labels):
    elems = elemdatabase.elementnr.keys()
    ratio=[]
    for z in elems:
        nz = labels.count(z)
        if nz > 0: ratio.append(nz)
    return ratio

######
def labels2electrons(labels):
    eleccount = 0
    for l in labels:
        eleccount += elemdatabase.elementnr[l]
    return eleccount

######
def get_element_count(labels: list, heavy_only: bool=False) -> np.ndarray:
    elems = list(elemdatabase.elementnr.keys())
    count = np.zeros((len(elems)),dtype=int)
    for l in labels:
        for jdx, elem in enumerate(elems):
            if l == elem:                             count[jdx] += 1
            if (l == 'H' or l == 'D') and heavy_only: count = 0
    return count

######
def get_adjacency_types(label: list, conmat: np.ndarray) -> np.ndarray:
    elems = elemdatabase.elementnr.keys()
    natoms = len(label)
    adjtypes = np.zeros((len(elems), len(elems)),dtype=int)
    found = np.zeros((natoms, natoms))

    for i in range(0, natoms):
        for j in range(i, natoms):
            if i != j:
                if (conmat[i, j] == 1) and (found[i, j] == 0):
                    for k, elem1 in enumerate(elems):
                        if label[i] == elem1:
                            for l, elem2 in enumerate(elems):
                                if label[j] == elem2:
                                    adjtypes[k, l] += 1
                                    if elem1 != elem2:
                                        adjtypes[l, k] += 1
                                    found[i, j] = 1
                                    found[j, i] = 1
                                    break
                            break
    return adjtypes

######
def get_radii(labels: list) -> np.ndarray:
    radii = []
    for l in labels:
        if l[-1].isdigit(): label = l[:-1]
        else: label = l
        radii.append(elemdatabase.CovalentRadius2[label])
    return np.array(radii)

######
def get_adjmatrix(labels: list, pos: list, cov_factor: float=1.3, metal_factor: float=1.0, adjust_factor: bool=False, radii="default", metal_only: bool=False, debug: int=0) -> Tuple[bool, np.array, np.array]:
    import numpy as np
    import warnings
    from scipy.spatial.distance import pdist, squareform
    """
    Compute the adjacency matrix of a molecule based on atomic labels and coordinates.

    Parameters
    ----------
    labels : list of str
        Atomic symbols (e.g., ['C', 'H', 'H', 'H', 'H']).
    pos : list of list of float or np.ndarray
        Atomic coordinates, shape (N, 3).
    cov_factor : float, optional
        Scaling factor applied to the sum of covalent radii when defining a bond (default 1.3).
    metal_factor : float, optional
        Multiplier for covalent radii of metal atoms (d- or f-block elements). Default is 1.0.
    radii : str, list, or np.ndarray, optional
        Covalent radii for atoms. If 'default', values are obtained from `get_radii(labels)`.
    metal_only : bool, optional
        If True, only metal–metal or metal–nonmetal bonds are included.
    adjust_factor : bool, optional
        If True, adaptively increases cov_factor until all atoms have at least one neighbor.
    debug : int, optional
        Print extra diagnostic information if > 0.
    max_factor : float, optional
        Upper bound for adaptive cov_factor search (default 2.5).
    factor_step : float, optional
        Step size for incrementing cov_factor (default 0.05).

    Returns
    -------
    isgood : bool
        False if atom–atom distance is below the clash threshold (molecule invalid).
    adjmat : np.ndarray of shape (N, N)
        Symmetric adjacency matrix (1 = bonded, 0 = not bonded).
    adjnum : np.ndarray of shape (N,)
        Number of bonded neighbors for each atom.
    cov_factor : float
        Final cov_factor used (may be higher than initial if adjust_factor=True).
    """
    isgood = True 
    clash_threshold = 0.3
    factor_step = 0.02
    max_factor = 1.5
    natoms = len(labels)
    adjmat = np.zeros((natoms, natoms))
    adjnum = np.zeros((natoms))

    # --- Load or validate covalent radii ---
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        if isinstance(radii, str) and radii == "default":
            radii = get_radii(labels)
        else:
            radii = np.array(radii, dtype=float)
    # --- Adjust radii for metal atoms ---
    if metal_factor != 1.0:
        for i, elem in enumerate(labels):
            block = elemdatabase.elementblock.get(elem, "")
            if block in ("d", "f"):  # transition or f-block metals
                radii[i] *= metal_factor
    # --- Precompute distance matrix ---
    distmat = squareform(pdist(pos))

    def compute_adjacency(cov_factor):
        """Helper to compute adjacency given a specific cov_factor."""
        # Compute the pairwise bonding distance threshold for each atom pair
        sum_radii = radii[:, None] + radii[None, :]    # all combinations of radii
        bond_threshold = cov_factor * sum_radii        # scale by cov_factor
        # Determine which atom pairs are bonded
        bonded = (distmat <= bond_threshold) & (distmat > clash_threshold)  ## Matrix of True/False
        adjmat = bonded.astype(int)                                         ## Converted into Integers
        # Handle metal-only case
        if metal_only:
            metal_mask = np.array([(elemdatabase.elementblock.get(lbl, "") in ("d", "f")) for lbl in labels],dtype=bool)
            # A bond is valid if either atom is a metal
            metal_bonds = np.logical_or.outer(metal_mask, metal_mask)
            adjmat *= metal_bonds
        # Zero diagonal
        np.fill_diagonal(adjmat, 0)
        # Check for clashes
        clash = np.any((distmat < clash_threshold) & (distmat > 0))
        return adjmat, clash

    # --- Adjust cov_factor if needed ---
    while True:
        adjmat, clash = compute_adjacency(cov_factor)
        adjnum = adjmat.sum(axis=1)
        if not clash:
            isgood = True
        else:
            isgood = False
            if debug > 0: print(f"GET_ADJMATRIX: Clash detected at {cov_factor=:.2f}")

        if adjust_factor and np.any(adjnum == 0) and cov_factor < max_factor:
            cov_factor += factor_step
            if debug > 0: print(f"GET_ADJMATRIX: Increasing cov_factor to {cov_factor:.2f} (some atoms have 0 neighbors)")
            continue
        break
    if debug and adjust_factor:
        if cov_factor >= max_factor: print(f"GET_ADJMATRIX: Reached max_factor={max_factor} but some atoms remain isolated.")
    return isgood, adjmat, adjnum.astype(int)

######
def inv(perm: list) -> list:
    inverse = [0] * len(perm)
    for i, p in enumerate(perm):
        inverse[p] = i
    return inverse

######
def get_blocks(matrix: np.ndarray) -> Tuple[list, list]:
    # retrieves the blocks from a diagonal block matrix
    startlist = []  # List including the starting atom for all blocks
    endlist = []    # List including the final atom for all blocks
    start = 1
    pos = start
    posold = 0
    blockcount = 0
    j = 1
    while j < len(matrix):
        if matrix[pos - 1, j] != 0.0: pos = j + 1
        if j == len(matrix) - 1:
            blockcount = blockcount + 1
            startlist.append(posold)
            endlist.append(pos - 1)
            posold = pos
            pos = pos + 1
            j = pos - 1
            continue
        j += 1
    if (blockcount == 0) and (len(matrix) == 1):  # if a 1x1 matrix is provided, it then finds 1 block
        startlist.append(0)
        endlist.append(0)
    return startlist, endlist

######
def compute_centroid(coord: list) -> list:
    natoms = len(coord)
    x = 0
    y = 0
    z = 0
    for idx in range(natoms):
        x += coord[idx][0]; y += coord[idx][1]; z += coord[idx][2]
    centroid = [float(x / natoms), float(y / natoms), float(z / natoms)]
    return centroid

######
def count_species(labels: list, pos: list, radii: list=None, indices: list=None, cov_factor: float=1.3, metal_factor: float=1.0, debug: int=0) -> Tuple[bool, list]:
    # Gets the covalent radii
    if radii is None:    radii = get_radii(labels)
    if indices is None:  indices = [*range(0,len(labels),1)]

    # Computes the adjacency matrix of what is received
    # isgood indicates whether the adjacency matrix could be built normally, or errors were detected. Typically, those errors are steric clashes
    isgood, adjmat, adjnum = get_adjmatrix(labels, pos, cov_factor, metal_factor, radii)
    if not isgood: return int(0)

    degree = np.diag(adjnum)  # creates a matrix with adjnum as diagonal values. Needed for the laplacian
    lap = adjmat - degree     # computes laplacian

    # creates block matrix
    graph = csr_matrix(lap)
    perm = reverse_cuthill_mckee(graph)
    gp1 = graph[perm, :]
    gp2 = gp1[:, perm]
    dense = gp2.toarray()

    # detects blocks in the block diagonal matrix called "dense"
    startlist, endlist = get_blocks(dense)

    nblocks = len(startlist)
    return nblocks

######
def split_species(labels: list, pos: list, radii: list=None, indices: list=None, cov_factor: float=1.3, metal_factor: float=1.0, debug: int=0) -> Tuple[bool, list]:
    ## Function that identifies connected groups of atoms from their atomic coordinates and labels.

    # Gets the covalent radii
    if radii is None:    radii = get_radii(labels)
    if indices is None:  indices = [*range(0,len(labels),1)]

    # Computes the adjacency matrix of what is received
    # isgood indicates whether the adjacency matrix could be built normally, or errors were detected. Typically, those errors are steric clashes
    isgood, adjmat, adjnum = get_adjmatrix(labels, pos, cov_factor, metal_factor, radii=radii)
    if not isgood: return None

    degree = np.diag(adjnum)  # creates a matrix with adjnum as diagonal values. Needed for the laplacian
    lap = adjmat - degree     # computes laplacian

    # creates block matrix
    graph = csr_matrix(lap)
    perm = reverse_cuthill_mckee(graph)
    gp1 = graph[perm, :]
    gp2 = gp1[:, perm]
    dense = gp2.toarray()

    # detects blocks in the block diagonal matrix called "dense"
    startlist, endlist = get_blocks(dense)

    nblocks = len(startlist)
    # keeps track of the atom movement within the matrix. Needed later
    atomlist = np.zeros((len(dense)))
    for b in range(0, nblocks):
        for i in range(0, len(dense)):
            if (i >= startlist[b]) and (i <= endlist[b]):
                atomlist[i] = b + 1
    invperm = inv(perm)
    atomlistperm = [int(atomlist[i]) for i in invperm]

    # assigns atoms to molecules
    blocklist = []
    for b in range(0, nblocks):
        atlist = []    # atom indices in the original ordering
        for i in range(0, len(atomlistperm)):
            if atomlistperm[i] == b + 1:
                atlist.append(indices[i])
        blocklist.append(atlist)
    return blocklist

#################################
def get_non_transition_metal_idxs(labels: list, debug: int=0):
    """ alkali metals (Group 1) and alkaline earth metals (Group 2)  """
    non_transition_metal_indices = []
    for idx, l in enumerate(labels):
        if elemdatabase.elementgroup[l]==1 and l != "H" and l != "D": # Alkali Metals
            non_transition_metal_indices.append(idx)
        elif elemdatabase.elementgroup[l]==2 :
            non_transition_metal_indices.append(idx) # Alkaline Earth Metals
        elif l in ["Al", "Ga", "Ge", "In", "Sn", "Tl", "Pb", "Bi", "Po", "At"]:
            non_transition_metal_indices.append(idx)
    return non_transition_metal_indices

#################################
def split_group(original_group, conn_idx, final_ligand_indices, debug: int=0):
    # Splits the "group" to obtain the groups connected to a specific metal
    from scope.classes_molecule import Group
    splitted_groups = []
    
    if debug > 1: print(f"GROUP.SPLIT_GROUP: {conn_idx=}")
    if debug > 1: print(f"GROUP.SPLIT_GROUP: {original_group.labels=}")
    if debug > 1: print(f"GROUP.SPLIT_GROUP: {original_group.coord=}")
    if debug > 1: print(f"GROUP.SPLIT_GROUP: {original_group.atoms=}")
    conn_labels       = extract_from_list(conn_idx, original_group.labels, dimension=1)
    conn_coord        = extract_from_list(conn_idx, original_group.coord, dimension=1)
    conn_frac_coord   = extract_from_list(conn_idx, original_group.frac_coord, dimension=1)
    conn_radii        = extract_from_list(conn_idx, original_group.radii, dimension=1)
    conn_atoms        = extract_from_list(conn_idx, original_group.atoms, dimension=1)
    if debug > 1: print(f"GROUP.SPLIT_GROUP: {conn_labels=}")

    cov_factor=original_group.get_parent("ligand").cov_factor
    blocklist = split_species(conn_labels, conn_coord, radii=conn_radii, cov_factor=cov_factor, debug=debug)      
    if debug > 0: print(f"GROUP.SPLIT_GROUP: {blocklist=}")

    ## Arranges Groups 
    for b in blocklist:
        if debug > 1: print(f"GROUP.SPLIT_GROUP: block={b}")
        gr_indices      = extract_from_list(b, conn_idx, dimension=1)
        ligand_idx      = extract_from_list(b, final_ligand_indices, dimension=1)
        gr_labels       = extract_from_list(b, conn_labels, dimension=1)
        gr_coord        = extract_from_list(b, conn_coord, dimension=1)
        gr_frac_coord   = extract_from_list(b, conn_frac_coord, dimension=1)
        gr_radii        = extract_from_list(b, conn_radii, dimension=1)
        gr_atoms        = extract_from_list(b, conn_atoms, dimension=1)
        # Create Group Object
        newgroup = Group(gr_labels, gr_coord, gr_frac_coord, radii=gr_radii)
        if debug > 1: print(f"GROUP.SPLIT_GROUP: {newgroup.labels=}")
        # For debugging
        newgroup.origin = "split_group"
        # Define the GROUP as parent of the group. Bottom-Up hierarchy
        newgroup.add_parent(original_group.get_parent("ligand"), indices=ligand_idx)
        # Pass the GROUP atoms to the groud
        newgroup.set_atoms(atomlist=gr_atoms)
        # Inherit the adjacencies from molecule
        newgroup.inherit_adjmatrix("ligand")
        # Associate the Groups with the Metals
        newgroup.get_connected_metals(debug=debug)
        newgroup.get_closest_metal(debug=debug)
        newgroup.get_hapticity(debug=debug)
        newgroup.get_denticity(debug=debug)
        # Top-down hierarchy
        splitted_groups.append(newgroup)
    return splitted_groups
