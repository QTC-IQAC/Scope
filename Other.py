import os, sys
import numpy as np

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

#############
def rmsd(labels1, coord1, labels2, coord2, reorder=True, debug: int=0):
    from Scope.Adapted_from_cell2mol import compute_centroid
    from Scope.Reconstruct import reorder_hungarian
    from Scope.Read_Write  import print_xyz

    """
    Compute the RMSD between two sets of coordinates
    Consider Expanding the data that is sent to hungarian. Currently only labels...
    ...but other data such as adjnum could be added 
    """
    ## Ensure both species have the same number of atoms
    if len(labels1) != len(labels2):
        raise ValueError("The number of atoms in the two lists must be the same.")

    ## Arranges data
    labels1 = np.array([f"{l}" for l in labels1])
    labels2 = np.array([f"{l}" for l in labels2])

    coords1 = np.array(coord1) - compute_centroid(coord1)
    coords2 = np.array(coord2) - compute_centroid(coord2)

    if reorder:
        idx = reorder_hungarian(labels1, labels2, coords1, coords2)
        coords2 = coords2[idx]
        labels2 = labels2[idx]

    U = kabsch_rotate(coords1, coords2)
    aligned_coords2 = np.dot(coords2, U.T)

    if debug > 0: print("Coordinates for self:")
    if debug > 0: print_xyz(labels1, coords1)
    if debug > 0: print("Aligned Coordinates for other:")
    if debug > 0: print_xyz(labels2, aligned_coords2)

    rmsd = np.sqrt(np.mean(np.sum((coords1 - aligned_coords2) ** 2, axis=1)))
    return np.round(rmsd, 4)

#############
def kabsch_rotate(coords1, coords2):
    """
    Compute the rotation matrix to align coords2 with coords1 using the Kabsch algorithm.
    """
    H = np.dot(coords2.T, coords1)
    V, S, W = np.linalg.svd(H)
    d = np.sign(np.linalg.det(V) * np.linalg.det(W))
    S = np.diag([1, 1, d])
    U = np.dot(np.dot(V, S), W)
    return U

###########
def correct_smiles_ligand(lig: object):
    ## Receives a ligand class object and constructs the smiles and the rdkit_mol object from scratch, using atoms and bond information

    Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    #### Creates an empty editable molecule
    rwlig = Chem.RWMol()

    # Adds atoms with their formal charge
    for jdx, a in enumerate(lig.atoms):
        rdkit_atom = Chem.Atom(lig.atnums[jdx])
        rdkit_atom.SetFormalCharge(int(a.charge))
        rdkit_atom.SetNoImplicit(True)
        rwlig.AddAtom(rdkit_atom)

    # Sets bond information and hybridization
    for jdx, a in enumerate(lig.atoms):
        nbonds = 0
        for bond in a.bond:
            nbonds += 1
            isaromatic = False
            if bond[2] == 1.0: btype = Chem.BondType.SINGLE
            elif bond[2] == 2.0: btype = Chem.BondType.DOUBLE
            elif bond[2] == 3.0: btype = Chem.BondType.TRIPLE
            elif bond[2] == 1.5:
                btype = Chem.BondType.AROMATIC
                rdkit_atom.SetIsAromatic(True)
            if bond[0] == jdx and bond[1] > jdx: rwlig.AddBond(bond[0], bond[1], btype)
    
        if nbonds == 1: hyb = Chem.HybridizationType.S
        elif nbonds == 2: hyb = Chem.HybridizationType.SP
        elif nbonds == 3: hyb = Chem.HybridizationType.SP2
        elif nbonds == 4: hyb = Chem.HybridizationType.SP3
        else: hyb = Chem.HybridizationType.UNSPECIFIED
        rdkit_atom.SetHybridization(hyb)

    # Creates Molecule
    obj = rwlig.GetMol()
    smiles = Chem.MolToSmiles(obj)
    
    Chem.SanitizeMol(obj)
    Chem.DetectBondStereochemistry(obj, -1)
    Chem.AssignStereochemistry(obj, flagPossibleStereoCenters=True, force=True)
    Chem.AssignAtomChiralTagsFromStructure(obj, -1)
    
    return smiles, obj

#######################################################
def get_dist (atom1_pos: list, atom2_pos: list) -> float :
    dist = np.linalg.norm(np.array(atom1_pos) - np.array(atom2_pos))
    return round(dist, 3)

def extract_from_list(entrylist: list, old_array: np.ndarray, dimension: int=2) -> np.ndarray:
    length = len(entrylist)
    if dimension == 2:
        new_array = np.empty((length, length))
        for idx, row in enumerate(entrylist):
            for jdx, col in enumerate(entrylist):
                new_array[idx, jdx] = old_array[row][col]
    elif dimension == 1:
        new_array = np.empty((length))
        for idx, val in enumerate(entrylist):
            new_array[idx] = old_array[val]
    return new_array

def where_in_array(array,condition) -> list:
    results = []
    for idx, a in enumerate(array):
        if a == condition: results.append(idx)
    return results

def mergelists(list1, list2, prop1, prop2):
    #print("Received", list1, list2)
    nitems=len(list1)+len(list2)
    mergedlist = []
    for idx in range(0,nitems):
        for jdx, at1 in enumerate(list1):
            if (idx == at1):
                mergedlist.append(prop1[jdx])
        for jdx, at2 in enumerate(list2):
            if (idx == at2):
                mergedlist.append(prop2[jdx])
    return mergedlist

def range2list(rang: range):
    lst = []
    for i in rang:
        lst.append(i)
    return lst

def get_metal_idxs(labels: list, debug: int=0):
    from Scope.Elementdata import ElementData
    elemdatabase = ElementData()
    metal_indices = []
    for idx, l in enumerate(labels):
        if (elemdatabase.elementblock[l] == 'd' or elemdatabase.elementblock[l] == 'f'): metal_indices.append(idx)
    return metal_indices

def get_metal_species(labels: list):
    from Scope.Elementdata import ElementData
    elemdatabase = ElementData()
    metal_species = []
    elems = list(set(labels))
    for idx, l in enumerate(elems):
        if l[-1].isdigit(): label = l[:-1]
        else: label = l
        if (elemdatabase.elementblock[label] == 'd' or elemdatabase.elementblock[label] == 'f') and l not in metal_species: metal_species.append(l)
    return metal_species

def pairwise(iterable):
    from itertools import tee
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def symmetrize_matrix(matrix, debug: int=0):
    assert np.shape(matrix)[0] == np.shape(matrix)[1]
    dim = np.shape(matrix)[0]
    sym_matrix = matrix.copy()
    for idx in range(dim):
        sym_matrix[idx,idx] = matrix[idx,idx]
        for jdx in range(idx+1,dim):
            if debug > 0: print(matrix[idx,jdx], matrix[jdx,idx])
            sym_matrix[idx,jdx] = (matrix[idx,jdx] + matrix[jdx,idx])/2.0
            sym_matrix[jdx,idx] = sym_matrix[idx,jdx]
    return sym_matrix 

def replace_zero(array): 
    
    for i in range(len(array)) :
        if array[i] == 0 : 
            array[i] = 1
    return array

def gram_schmidt(A,norm=True,row_vect=False):
    """Orthonormalizes vectors by gram-schmidt process
    
    Parameters
    -----------
    A : ndarray,
    Matrix having vectors in its columns
    
    norm : bool,
    Do you need Normalized vectors?
    
    row_vect: bool,
    Does Matrix A has vectors in its rows?
    
    Returns 
    -------
    G : ndarray,
    Matrix of orthogonal vectors 
    
    """
    if row_vect :
        # if true, transpose it to make column vector matrix
        A = A.T
    
    no_of_vectors = A.shape[1]
    G = A[:,0:1].copy() # copy the first vector in matrix
    # 0:1 is done to to be consistent with dimensions - [[1,2,3]]
    
    # iterate from 2nd vector to number of vectors
    for i in range(1,no_of_vectors):
        
        # calculates weights(coefficents) for every vector in G
        numerator = A[:,i].dot(G)
        denominator = np.diag(np.dot(G.T,G)) #to get elements in diagonal
        weights = np.squeeze(numerator/denominator)
        
        # projected vector onto subspace G 
        projected_vector = np.sum(weights * G,
                                  axis=1,
                                  keepdims=True)
        
        # orthogonal vector to subspace G
        orthogonalized_vector = A[:,i:i+1] - projected_vector
        
        # now add the orthogonal vector to our set 
        G = np.hstack((G,orthogonalized_vector))
        
    if norm :
        # to get orthoNORMAL vectors (unit orthogonal vectors)
        # replace zero to 1 to deal with division by 0 if matrix has 0 vector
        G = G/replace_zero(np.linalg.norm(G,axis=0))
    
    if row_vect:
        return G.T
    
    return G

