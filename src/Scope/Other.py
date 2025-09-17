import os, sys
import numpy as np

####
def check_convergence(values: list, current_step: int=None, thres: float=1e-5, debug: int=0):
    ## This function checks if a series of values has converged, based on the difference between the last two values
    if debug > 0: print(f"CHECK_CONVERGENCE: received {values=}")
    if current_step is None:
        ## None when you want the overall convergence, not that of a given step
        current_step = -1
        for v in values:
            if v != float(0.0): current_step += 1
        if debug > 0: print(f"CHECK_CONVERGENCE: {current_step=}")
        if debug > 0: print(f"CHECK_CONVERGENCE: difference: {values[current_step-1]-values[current_step]}")
    if np.abs(values[current_step-1]-values[current_step]) > thres: 
        if debug > 0: print(f"CHECK_CONVERGENCE: difference above threshold {thres}. Not converged")
        return False
    else: 
        if debug > 0: print(f"CHECK_CONVERGENCE: difference below threshold {thres}. Converged")
        return True

####
def overlap_molecules(labels1, coords1, labels2, coords2, center_method: str="centroid", use_ext_info: bool=True, translate_to_ref: bool=True, save_to_folder: str="/Users/sergivela/Documents/SCOPE/Program/tests/Overlap_Molecules/", debug: int=0):
    from .Connectivity import compute_centroid, get_adjmatrix
    from .Reconstruct import reorder_hungarian
    from .Read_Write  import print_xyz, write_xyz
    from collections import Counter

    ## Ensure both species have the same number of atoms
    if len(labels1) != len(labels2):
        raise ValueError("The number of atoms in the two lists must be the same.")

    ## Stores the original adjacency matrices
    isgood1, orig_adjmat1, orig_adjnum1 = get_adjmatrix(labels1, coords1)
    isgood2, orig_adjmat2, orig_adjnum2 = get_adjmatrix(labels2, coords2)

    if debug > 0: print(f"{isgood1}, {isgood2}")
    if debug > 0: print(f"original adjmat1={orig_adjmat1[0]}")
    if debug > 0: print(f"original adjmat2={orig_adjmat2[0]}")

    if center_method == "centroid":
        center1 = np.array(compute_centroid(coords1)) 
        center2 = np.array(compute_centroid(coords2))
    elif center_method == "metal":
        metal1 = get_metal_idxs(labels1)[0]
        metal2 = get_metal_idxs(labels2)[0]
        center1 = np.array(coords1)[metal1]
        center2 = np.array(coords2)[metal2]

    coords1c = np.array(coords1) - center1
    coords2c = np.array(coords2) - center2

    if debug > 0:
        print("---------------------")
        print("Centered Coords of 1:")
        print("---------------------")
        print_xyz(labels1, coords1c)
        write_xyz(save_to_folder+"centered1.xyz", labels1, coords1c)
        print("---------------------")
        print("Centered Coords of 2:")
        print("---------------------")
        print_xyz(labels2, coords2c)
        write_xyz(save_to_folder+"centered2.xyz", labels2, coords2c)

    ## We do a first alignment, to facilitate the hungarian
    _, _, coords2a, _ = kabsch_align(labels2, coords2c, labels1, coords1c, center_method=center_method, debug=debug)
    if debug > 0:
        print("--------------------")
        print("Aligned Coords of 2:")
        print("--------------------")
        print_xyz(labels2, coords2a)
        write_xyz(save_to_folder+"aligned2.xyz", labels2, coords2a)

    if use_ext_info:
        ## Decorate the labels to facilitate the reorder hungarian
        data1 = get_extended_info(labels1, coords1c, orig_adjmat1, orig_adjnum1, debug=debug)
        data2 = get_extended_info(labels2, coords2a, orig_adjmat2, orig_adjnum2, debug=debug)
        unique1 = [d for d in np.unique(data1) if Counter(data1)[d] == 1]    ## Data's which appear only once, and thus have no possibility of error
        unique2 = [d for d in np.unique(data2) if Counter(data2)[d] == 1]
        ## These are going to be all atoms whose data appears once in both structures
        safe = [index for index in range(len(data1)) if data1[index] in unique1 and data2[index] in unique2]
        if debug > 0: print(f"{safe=}")

        ## For each atom in SAFE, we add the topological distance of all atoms to their data
        for s in safe:
            if debug > 0: print(f"Taking atom {s=} with {data1[s]=} as reference when computing topological distances")
            dist1 = compute_topological_distances(orig_adjmat1, s)
            dist2 = compute_topological_distances(orig_adjmat2, s)
            data1 = np.array([str(f"{d1}{d2}") for d1, d2 in zip(data1, dist1)])
            data2 = np.array([str(f"{d1}{d2}") for d1, d2 in zip(data2, dist2)])
            labels1 = np.array([f"{l}" for l in labels1])
            labels2 = np.array([f"{l}" for l in labels2])
    else:
        data1 = np.array([f"{l}" for l in labels1])
        data2 = np.array([f"{l}" for l in labels2])
        labels1 = np.array([f"{l}" for l in labels1])
        labels2 = np.array([f"{l}" for l in labels2])

    if debug > 0: print(f"OVERLAP_MOLECULES: {data1=}")
    if debug > 0: print(f"OVERLAP_MOLECULES: {data2=}")
    # Verify that data1 and data2 contain the same items
    set1 = set(data1)
    set2 = set(data2)
    if set1 != set2 and debug > 0:
        print("OVERLAP_MOLECULES: WARNING: data1 and data2 do not contain the same items.")
        print("Items only in data1:", set1 - set2)
        print("Items only in data2:", set2 - set1)
    else:
        if debug > 0: print("data1 and data2 contain the same items.")

    ## Now we have the relevant data. We proceed to reorder
    if debug > 0: print("-------------------------")
    if debug > 0: print("## STARTING TO REORDER ##")
    if debug > 0: print("-------------------------")
    map12 = reorder_hungarian(data1, data2, coords1c, coords2a, debug=debug)
    if debug > 0: print(f"GET_RMSD: reorder map={map12}")
    labels2r = np.array(labels2)[map12]
    coords2r = coords2c[map12]
    data2    = data2[map12]
    if debug > 0:
        print("--------------------------")
        print("Coords of 1 after reorder:")
        print("--------------------------")
        print_xyz(labels1, coords1c)
        write_xyz(save_to_folder+"reorder1.xyz", labels1, coords1c)
        print("--------------------------")
        print("Coords of 2 after reorder:")
        print("--------------------------")
        print_xyz(labels2r, coords2r)
        write_xyz(save_to_folder+"reorder2.xyz", labels2r, coords2r)

    #### After this reorder, we check if any adjacency matrix entry is different
    adjmat2r = orig_adjmat2[np.ix_(map12, map12)]
    different = []
    isgood = True
    for idx in range(len(orig_adjmat1)):
        if any(orig_adjmat1[idx] != adjmat2r[idx]): different.append(idx)
    if len(different) > 0: 
        isgood = False
        if debug > 0: 
            print("OVERLAP_MOLECULES: WARNING! Different entries in the adjacency matrix:")
            for d in different:
                print(d, orig_adjmat1[d], adjmat2r[d])
    
    #### Then, we do the final alignment, which will be even better than the first
    _, _, coords2ra, _ = kabsch_align(labels2r, coords2r, labels1, coords1c, center_method=center_method, debug=debug)
    if debug > 0:
        print("----------------------------------")
        print("Coords of 2 after final alignment:")
        print("----------------------------------")
        print_xyz(labels2r, coords2ra)
        write_xyz(save_to_folder+"final2.xyz", labels2r, coords2ra)
    
    #### And finally, we translate to the original center of coordinates of the reference molecule 
    if translate_to_ref:
        coords1f = coords1c  + center1
        coords2f = coords2ra + center1
        write_xyz(save_to_folder+"final_translated1.xyz", labels1, coords1f)
        write_xyz(save_to_folder+"final_translated2.xyz", labels2r, coords2f)

    return isgood, labels1, coords1f, labels2r, coords2f, map12

####
def get_extended_info(labels, coords, adjmat, adjnum, debug: int=0):
    if debug > 0: print(f"GET_EXT_INFO, received {adjmat[0]} for first atom")
    ## Function to add decorators to labels
    bonded = []
    for jdx in range(len(labels)):
        labs = []
        for kdx in range(len(labels)):
           if jdx != kdx and adjmat[jdx, kdx] > 0: labs.append(labels[kdx])
        bonded.append(''.join(sorted(labs)))
    data = np.array([str(f"{l}{a}{b}") for l, a, b in zip(labels, adjnum, bonded)])
    if debug > 0: print(f"GET_EXT_INFO: {data=}")
    return data

#####
def compute_topological_distances(adj_matrix: np.ndarray, ref_atom: int) -> dict:
    import numpy as np
    import networkx as nx
    """
    Compute topological distances (graph distances) from a reference atom to all others.

    Parameters:
        adj_matrix: (N, N) binary numpy array
        ref_atom: index of the reference atom (0-based)

    Returns:
        distances: dict {atom_index: topological_distance}
    """
    G = nx.Graph()
    N = adj_matrix.shape[0]

    # Add edges based on adjacency matrix
    for i in range(N):
        for j in range(i+1, N):
            if adj_matrix[i, j] > 0:
                G.add_edge(i, j)
                #G.add_edge(j, i)

    # Compute shortest path lengths from ref_atom
    dist_dict = nx.single_source_shortest_path_length(G, ref_atom)

    distances = np.full(N, -1, dtype=int)
    for idx, d in dist_dict.items():
        distances[idx] = d

    return distances

#############
def rmsd(labels1, coords1, labels2, coords2, reorder: bool=False, center_method='centroid', atom_idxs: list=None, debug: int=0):
    assert len(labels1) == len(labels2)
    #coords1 = np.asarray(coords1)
    #coords2 = np.asarray(coords2)
    if atom_idxs is not None: atom_idxs = np.sort(np.asarray(atom_idxs))
    if reorder: 
        isgood, new_labels1, new_coords1, new_labels2, new_coords2, _ = overlap_molecules(labels1, coords1, labels2, coords2, center_method=center_method, debug=debug)
        if not isgood: return None
        assert all(new_labels1 == new_labels2)       ## In principle, one wants to compute the RMSD of atoms with the same ordering 
        if atom_idxs is None: rmsd = np.sqrt(np.mean(np.sum((new_coords1 - new_coords2)**2, axis=1)))
        else:                 rmsd = np.sqrt(np.mean(np.sum((new_coords1[atom_idxs] - new_coords2[atom_idxs])**2, axis=1))) 
    else:        
        if atom_idxs is None: rmsd = np.sqrt(np.mean(np.sum((coords1 - coords2)**2, axis=1)))
        else:                 rmsd = np.sqrt(np.mean(np.sum((coords1[atom_idxs] - coords2[atom_idxs])**2, axis=1))) 
    return np.round(rmsd, 4)

#############
def kabsch_align(LP, P, LQ, Q, center_method: str="centroid", debug: int=0):
    from .Connectivity import compute_centroid
    """
    Perform the Kabsch algorithm to find the optimal rotation and translation
    to align point set P (mobile) to Q (target/reference).

    Parameters:
        P: np.ndarray of shape (N, 3) - mobile coordinates
        Q: np.ndarray of shape (N, 3) - reference coordinates

    Returns:
        R: Rotation matrix (3x3)
        t: Translation vector (3,)
        P_aligned: Transformed coordinates of P (after alignment)
        rmsd: Root Mean Square Deviation after alignment
    """

    if center_method == "centroid":
        centroid_P = np.array(compute_centroid(P))
        centroid_Q = np.array(compute_centroid(Q))
        # Step 1: Center both point sets
        P_centered = P - centroid_P
        Q_centered = Q - centroid_Q
    elif center_method == "metal":
        metalP = get_metal_idxs(LP)[0]
        metalQ = get_metal_idxs(LQ)[0]
        centroid_P = P[metalP]
        centroid_Q = Q[metalQ]
        # Step 1: Center both point sets
        P_centered = P - centroid_P
        Q_centered = Q - centroid_Q

    # Step 2: Covariance matrix
    H = P_centered.T @ Q_centered

    # Step 3: SVD
    U, S, Vt = np.linalg.svd(H)

    # Step 4: Compute rotation matrix
    R = Vt.T @ U.T

    # Reflection correction (to ensure a right-handed coordinate system)
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = Vt.T @ U.T

    # Step 5: Compute translation
    t = centroid_Q - R @ centroid_P

    # Step 6: Apply transformation
    P_aligned = (R @ P.T).T + t

    # Step 7: Compute RMSD
    diff = P_aligned - Q
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    return R, t, P_aligned, rmsd

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

######
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

#def extract_from_list(entrylist: list, old_array: np.ndarray, dimension: int=2) -> np.ndarray:
#    length = len(entrylist)
#    if dimension == 2:
#        new_array = np.empty((length, length))
#        for idx, row in enumerate(entrylist):
#            for jdx, col in enumerate(entrylist):
#                new_array[idx, jdx] = old_array[row][col]
#    elif dimension == 1:
#        new_array = np.empty((length))
#        for idx, val in enumerate(entrylist):
#            new_array[idx] = old_array[val]
#    return new_array

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

####
def furthest_point_sampling(data, k, dist_func):
    from Scope.VNM_tools import pairwise_distance_matrix
    """
    Selects `k` furthest points from `data` using a custom distance function.

    Parameters:
    - data: np.ndarray of shape (n_samples, n_features)
    - k: number of samples to select
    - dist_func: function to compute distance between two samples

    Returns:
    - selected_indices: list of indices of selected points
    """
    data        = np.array(data)
    #freqs_cm1   = np.array(freqs_cm1)

    n_samples   = data.shape[0]
    #D = pairwise_distance_matrix(data, freqs_cm1, dist_func)
    D = pairwise_distance_matrix(data, dist_func)

    selected      = [0]  # Starts with 0
    min_distances = D[selected[0]].copy()         # Distance from selected[0] to all points

    for _ in range(1, k):
        # At each step, pick the point with the maximum minimum distance to the selected set
        for i in range(n_samples):
            min_distances[i] = min(min_distances[i], min(D[i, j] for j in selected))
        next_idx = np.argmax(min_distances)
        selected.append(next_idx)
    return selected
