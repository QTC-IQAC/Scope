import warnings
import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
from typing import Tuple
from Scope.Elementdata import ElementData  
elemdatabase = ElementData()

#from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal


################################
def extract_from_list(entrylist: list, old_array: list, dimension: int=2, debug: int=0) -> list:
    if debug >= 1: print(f"EXTRACT_FROM_LIST. received: {entrylist=}")
    if debug >= 1: print(f"EXTRACT_FROM_LIST. received: {old_array=}")
    if debug >= 1: print(f"EXTRACT_FROM_LIST. maximum value received in entrylist: {np.max(entrylist)+1}")
    if debug >= 1: print(f"EXTRACT_FROM_LIST. length of old_array: {len(old_array)}")
    assert len(old_array) >= np.max(entrylist)+1
    length = len(entrylist)
    if dimension == 2:
        new_array = np.empty((length, length), dtype=object)
        for idx, row in enumerate(entrylist):
            for jdx, col in enumerate(entrylist):
                new_array[idx, jdx] = old_array[row][col]
    elif dimension == 1:
        new_array = np.empty((length), dtype=object)
        for idx, val in enumerate(entrylist):
            new_array[idx] = old_array[val]
    return list(new_array)

################################
def printxyz(labels, pos):
    print(len(labels))
    print("")
    for idx, l in enumerate(labels):
        print("%s  %.6f  %.6f  %.6f" % (l, pos[idx][0], pos[idx][1], pos[idx][2]))

################################
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

################################
def labels2ratio(labels):
    elems = elemdatabase.elementnr.keys()
    ratio=[]
    for z in elems:
        nz = labels.count(z)
        if nz > 0: ratio.append(nz)
    return ratio

################################
def labels2electrons(labels):
    eleccount = 0
    for l in labels:
        eleccount += elemdatabase.elementnr[l]
    return eleccount 

################################
def get_element_count(labels: list, heavy_only: bool=False) -> np.ndarray:
    elems = list(elemdatabase.elementnr.keys())
    count = np.zeros((len(elems)),dtype=int)
    for l in labels:
        for jdx, elem in enumerate(elems):
            if l == elem:                             count[jdx] += 1
            if (l == 'H' or l == 'D') and heavy_only: count = 0
    return count

################################
def get_adjacency_types(label: list, conmat: np.ndarray) -> np.ndarray:
    elems = elemdatabase.elementnr.keys()
    natoms = len(label)
    bondtypes = np.zeros((len(elems), len(elems)),dtype=int)
    found = np.zeros((natoms, natoms))

    for i in range(0, natoms):
        for j in range(i, natoms):
            if i != j:
                if (conmat[i, j] == 1) and (found[i, j] == 0):
                    for k, elem1 in enumerate(elems):
                        if label[i] == elem1:
                            for l, elem2 in enumerate(elems):
                                if label[j] == elem2:
                                    bondtypes[k, l] += 1
                                    if elem1 != elem2:
                                        bondtypes[l, k] += 1
                                    found[i, j] = 1
                                    found[j, i] = 1
                                    break
                            break
    return bondtypes

################################
def get_radii(labels: list) -> np.ndarray:
    radii = []
    for l in labels:
        if l[-1].isdigit(): label = l[:-1]
        else: label = l
        radii.append(elemdatabase.CovalentRadius2[label])
    return np.array(radii)

####################################
def get_adjmatrix(labels: list, pos: list, cov_factor: float=1.3, metal_factor: float=1.0, radii="default", metal_only: bool=False, debug: int=0) -> Tuple[int, list, list]:
    isgood = True 
    clash_threshold = 0.3
    natoms = len(labels)
    adjmat = np.zeros((natoms, natoms))
    adjnum = np.zeros((natoms))

    # Sometimes argument radii np.ndarry, or list
    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore", category=FutureWarning)
        if type(radii) == str:
            if radii == "default":
                radii = get_radii(labels)

    # Modifies the radii for metal atoms, if necessary
    if metal_factor != 1.0:
        for idx, r in enumerate(radii):
            if (elemdatabase.elementblock[labels[idx]] == "d" or elemdatabase.elementblock[labels[idx]] == "f"):
                radii[idx] = (r * metal_factor)  # the covalent radii of the metal is modified

    # Creates Adjacency Matrix
    for i in range(0, natoms - 1):
        for j in range(i, natoms):
            if i != j:
                a = np.array(pos[i])
                b = np.array(pos[j])
                dist = np.linalg.norm(a - b)
                thres = (radii[i] + radii[j]) * cov_factor
                if dist <= clash_threshold:
                    isgood = False # invalid molecule
                    if debug > 0: print("Adjacency Matrix: Distance", dist, "smaller than clash for atoms", i, j)
                elif dist <= thres:
                    if not metal_only:
                        adjmat[i, j] = 1
                        adjmat[j, i] = 1
                    if metal_only:
                        if (elemdatabase.elementblock[labels[i]] == "d"
                        or elemdatabase.elementblock[labels[i]] == "f"
                        or elemdatabase.elementblock[labels[j]] == "d"
                        or elemdatabase.elementblock[labels[j]] == "f"):
                            adjmat[i, j] = 1
                            adjmat[j, i] = 1

    # Sums the adjacencies of each atom to obtain "adjnum" 
    for i in range(0, natoms):
        adjnum[i] = np.sum(adjmat[i, :])

    adjmat = adjmat.astype(int)
    adjnum = adjnum.astype(int)
    return isgood, adjmat, adjnum

####################################
def inv(perm: list) -> list:
    inverse = [0] * len(perm)
    for i, p in enumerate(perm):
        inverse[p] = i
    return inverse

####################################
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

####################################
#def get_molecules(labels: list, pos: list, factor: float=1.3, debug: int=0) -> Tuple[bool, list]:
#    ## Function that identifies connected groups of atoms from their atomic coordinates and labels.
#
#    # Gets the covalent radii
#    radii = get_radii(labels)
#
#    # Computes the adjacency matrix of what is received
#    isgood, adjmat, adjnum = get_adjmatrix(labels, pos, factor, radii)
#
#    # isgood indicates whether the adjacency matrix could be built normally, or errors were detected. Typically, those errors are steric clashes
#    if isgood:
#        degree = np.diag(adjnum)  # creates a matrix with adjnum as diagonal values. Needed for the laplacian
#        lap = adjmat - degree     # computes laplacian
#
#        # creates block matrix
#        graph = csr_matrix(lap)
#        perm = reverse_cuthill_mckee(graph)
#        gp1 = graph[perm, :]
#        gp2 = gp1[:, perm]
#        dense = gp2.toarray()
#
#        # detects blocks in the block diagonal matrix called "dense"
#        startlist, endlist = get_blocks(dense)
#
#        nmolec = len(startlist)
#        # keeps track of the atom movement within the matrix. Needed later
#        atomlist = np.zeros((len(dense)))
#        for b in range(0, nmolec):
#            for i in range(0, len(dense)):
#                if (i >= startlist[b]) and (i <= endlist[b]):
#                    atomlist[i] = b + 1
#        invperm = inv(perm)
#        atomlistperm = [int(atomlist[i]) for i in invperm]
#
#        # assigns atoms to molecules
#        moleclist = []
#        for b in range(0, nmolec):
#            labelist = []
#            poslist = []
#            atlist = []    # atom indices in the original ordering
#            for i in range(0, len(atomlistperm)):
#                if atomlistperm[i] == b + 1:
#                    labelist.append(labels[i])
#                    poslist.append(pos[i])
#                    atlist.append(i)
#            radiilist = get_radii(labelist)
#
#            newmolec = molecule(str(b), atlist, labelist, poslist, radiilist)
#            newmolec.information(1.3,1.0)
#
#            moleclist.append(newmolec)
#            #moleclist.append(list([atlist, labelist, poslist, radiilist]))
#    return moleclist

################################
def compute_centroid(coord: list) -> list:
    natoms = len(coord)
    x = 0
    y = 0
    z = 0
    for idx in range(natoms):
        x += coord[idx][0]; y += coord[idx][1]; z += coord[idx][2]
    centroid = [float(x / natoms), float(y / natoms), float(z / natoms)]
    return np.array(centroid)
#########################

#########################
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
#########################

####################################
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

#####################
def merge_atoms(atoms):
    labels = [] 
    coords  = [] 
    for a in atoms:
        labels.append(a.label) 
        coords.append(a.coord) 
    return labels, coords

#################################
def compare_atoms(at1, at2, check_coordinates: bool=False, debug: int=0):
    if debug > 0: 
        print("Comparing Atoms")
        print(at1)
        print(at2)
    # Compares Species, Coordinates, Charge and Spin
    if (at1.label != at2.label): return False
    if check_coordinates:
        if (at1.coord[0] != at2.coord[0]): return False
        if (at1.coord[1] != at2.coord[1]): return False
        if (at1.coord[2] != at2.coord[2]): return False
    if hasattr(at1,"charge") and hasattr(at2,"charge"):
        if (at1.charge != at2.charge): return False
    if hasattr(at1,"spin") and hasattr(at2,"spin"):
        if (at1.spin != at2.spin): return False
    return True

#################################
def compare_metals (at1, at2, check_coordinates: bool=False, debug: int=0):
    if debug > 0: 
        print("COMPARE_METALS. Comparing:")
        print(at1.label)
        print(at2.label)

    if at1.subtype != "metal" or at2.subtype != "metal": 
        if debug > 0: print("COMPARE_METALS. Different subtype")
        if debug > 0: print(at1.subtype)
        if debug > 0: print(at1.subtype)
        return False

    if (at1.label != at2.label): 
        if debug > 0: print("COMPARE_METALS. Different label")
        return False

    if not hasattr(at1,"coord_sphere_formula"): at1.get_coord_sphere_formula()
    if not hasattr(at2,"coord_sphere_formula"): at2.get_coord_sphere_formula()
    if (at1.coord_sphere_formula != at2.coord_sphere_formula):
        if debug > 0: print("COMPARE_METALS. Different coordination sphere")
        if debug > 0: print(at1.coord_sphere_formula)
        if debug > 0: print(at2.coord_sphere_formula)
        return False
    
    if check_coordinates:
        if (at1.coord[0] != at2.coord[0]): return False
        if (at1.coord[1] != at2.coord[1]): return False
        if (at1.coord[2] != at2.coord[2]): return False
        
    return True

#################################
def compare_species(mol1, mol2, debug: int=0):
    if debug > 0: 
        print("Comparing Species")
        print(mol1)
        print(mol2)
    # a pair of species is compared on the basis of:
    # 1) the total number of atoms
    if (mol1.natoms != mol2.natoms): return False

    # 2) the total number of electrons (as sum of atomic number)
    if (mol1.eleccount != mol2.eleccount): return False

    # 3) the number of atoms of each type
    if not hasattr(mol1,"element_count"): mol1.set_element_count()
    if not hasattr(mol2,"element_count"): mol2.set_element_count()
    for kdx, elem in enumerate(mol1.element_count):
        if elem != mol2.element_count[kdx]: return False       

    # 4) the number of adjacencies between each pair of element types
    if not hasattr(mol1,"adj_types"):     mol1.set_adj_types()
    if not hasattr(mol2,"adj_types"):     mol2.set_adj_types()
    for kdx, elem in enumerate(mol1.adj_types):
        for ldx, elem2 in enumerate(elem):
            if elem2 != mol2.adj_types[kdx, ldx]: return False
    return True
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
    ## Yuri's function. I don't know why it is not in group class
    from Scope.Classes_Molecule import group
    # Split the "group" to obtain the groups connected to a specific metal
    splitted_groups = []
    
    if debug > 1: print(f"GROUP.SPLIT_GROUP: {conn_idx=}")
    if debug > 1: print(f"GROUP.SPLIT_GROUP: {original_group.labels=}")
    if debug > 1: print(f"GROUP.SPLIT_GROUP: {original_group.coord=}")
    if debug > 1: print(f"GROUP.SPLIT_GROUP: {original_group.atoms=}")
    conn_labels  = extract_from_list(conn_idx, original_group.labels, dimension=1)
    conn_coord   = extract_from_list(conn_idx, original_group.coord, dimension=1)
    conn_frac_coord   = extract_from_list(conn_idx, original_group.frac_coord, dimension=1)
    conn_radii   = extract_from_list(conn_idx, original_group.radii, dimension=1)
    conn_atoms   = extract_from_list(conn_idx, original_group.atoms, dimension=1)
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
        newgroup = group(gr_labels, gr_coord, gr_frac_coord, radii=gr_radii)
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
