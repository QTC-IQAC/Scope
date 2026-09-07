import numpy as np
from scope.overlap import get_extended_info, kabsch_align, kabsch_rotate, overlap_molecules, rmsd

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
def get_metal_idxs(labels: list, debug: int=0):
    from scope.elementdata import ElementData
    elemdatabase = ElementData()
    metal_indices = []
    for idx, l in enumerate(labels):
        if (elemdatabase.elementblock[l] == 'd' or elemdatabase.elementblock[l] == 'f'): metal_indices.append(idx)
    return metal_indices

####
def get_metal_species(labels: list):
    from scope.elementdata import ElementData
    elemdatabase = ElementData()
    metal_species = []
    elems = list(set(labels))
    for idx, l in enumerate(elems):
        if l[-1].isdigit(): label = l[:-1]
        else: label = l
        if (elemdatabase.elementblock[label] == 'd' or elemdatabase.elementblock[label] == 'f') and l not in metal_species: metal_species.append(l)
    return metal_species

####
def pairwise(iterable):
    from itertools import tee
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

####
def replace_zero(array): 
    for i in range(len(array)) :
        if array[i] == 0 : 
            array[i] = 1
    return array

####
def gram_schmidt(A,norm=True,row_vect=False):
    """
    Orthonormalize vectors with the Gram-Schmidt procedure.

    Parameters:
        A (ndarray):                   Input matrix of vectors.
        norm (bool):                   Whether to normalize the output vectors.
        row_vect (bool):               Whether vectors are stored by rows.

    Returns:
        ndarray: Orthogonalized vector matrix.
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
    from scope.vnm_tools import pairwise_distance_matrix
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
