import os, sys
import numpy as np

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

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

