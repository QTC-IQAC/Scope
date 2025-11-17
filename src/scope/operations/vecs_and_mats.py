import numpy as np

def add_item(item, vector):
    if item not in vector:
        vector.append(item)
    return vector

def normalize(arr):
    arr = np.asarray(arr)
    if arr.ndim == 1:
        norm = np.linalg.norm(arr)
        if norm == 0: norm = 1.0
        return arr / norm
    else:
        norms = np.linalg.norm(arr, axis=1)
        norms[norms == 0] = 1.0
        return arr / norms[:, None]

def determinant(matrix):
    return np.linalg.det(matrix)

def symmetrize(matrix, debug: int=0):
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
