import numpy as np

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
