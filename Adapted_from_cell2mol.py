import warnings
import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
from typing import Tuple
from Test_V3.Elementdata import ElementData  
elemdatabase = ElementData()

################################
def labels2formula(labels):
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
def get_elementcount(labels: list) -> np.ndarray:
    elems = elemdatabase.elementnr.keys()
    times = np.zeros((len(elems)),dtype=int)
    for l in labels:
        for jdx, elem in enumerate(elems):
            if l == elem:
                times[jdx] += 1
    return times

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
def get_adjmatrix(labels: list, pos: list, factor: float, radii="default") -> Tuple[int, list, list]:
    status = 1  # good molecule, no clashes yet
    clash_threshold = 0.3
    natoms = len(labels)
    adjmat = np.zeros((natoms, natoms))
    adjnum = np.zeros((natoms))

    # Sometimes argument radii np.ndarry, or list
    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore", category=FutureWarning)
        if radii == "default":
            radii = get_radii(labels)

    # Creates Adjacency Matrix
    for i in range(0, natoms - 1):
        for j in range(i, natoms):
            if i != j:
                a = np.array(pos[i])
                b = np.array(pos[j])
                dist = np.linalg.norm(a - b)
                thres = (radii[i] + radii[j]) * factor
                if dist <= clash_threshold:
                    status = 0  # invalid molecule
                    print("Adjacency Matrix: Distance", dist, "smaller than clash for atoms", i, j)
                elif dist <= thres:
                    adjmat[i, j] = 1
                    adjmat[j, i] = 1

    # Sums the adjacencies of each atom to obtain "adjnum" 
    for i in range(0, natoms):
        adjnum[i] = np.sum(adjmat[i, :])

    adjmat = adjmat.astype(int)
    adjnum = adjnum.astype(int)
    return status, adjmat, adjnum

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
def get_molecules(labels: list, pos: list, factor: float=1.3, debug: int=0) -> Tuple[bool, list]:
    ## Function that identifies connected groups of atoms from their atomic coordinates and labels.

    # Gets the covalent radii, and modifies that of the metal if necessary
    radii = get_radii(labels)

    # Computes the adjacency matrix of what is received
    status, adjmat, adjnum = get_adjmatrix(labels, pos, factor, radii)

    # status indicates whether the adjacency matrix could be built normally, or errors were detected. Typically, those errors are steric clashes
    if status == 1:
        Warning = False
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

        nmolec = len(startlist)
        # keeps track of the atom movement within the matrix. Needed later
        atomlist = np.zeros((len(dense)))
        for b in range(0, nmolec):
            for i in range(0, len(dense)):
                if (i >= startlist[b]) and (i <= endlist[b]):
                    atomlist[i] = b + 1
        invperm = inv(perm)
        atomlistperm = [int(atomlist[i]) for i in invperm]

        # assigns atoms to molecules
        moleclist = []
        for b in range(0, nmolec):
            labelist = []
            poslist = []
            atlist = []    # atom indices in the original ordering
            for i in range(0, len(atomlistperm)):
                if atomlistperm[i] == b + 1:
                    labelist.append(labels[i])
                    poslist.append(pos[i])
                    atlist.append(i)
            radiilist = get_radii(labelist)

            moleclist.append(list([atlist, labelist, poslist, radiilist]))
    return Warning, moleclist
