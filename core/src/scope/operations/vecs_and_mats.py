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

def gcd(a, b):
    # gcd = Greatest common divisor
    while b != 0:
        a, b = b, a % b
    return a

def gcd_list(numbers):
    # returns the gcd for a list of numbers
    result = numbers[0]
    for n in numbers[1:]:
        result = gcd(result, n)
    return result

def gaussian(grid, center, sigma=0.2, normalize=True):
    """
    Parameters
    ----------
    grid : array_like           Grid/range where the Gaussian is evaluated (any units).
    center : float              Center of the Gaussian in the same units as E.
    sigma : float, optional     Standard deviation (width) in the same units as E. Default is 0.2.
    normalize : bool, optional  Whether to normalize the Gaussian to unit area. Default is True.

    Returns
    -------
    g: array_like               Gaussian profile evaluated on range.
    """
    if sigma <= 0: raise ValueError("sigma must be > 0")

    g = np.exp(-(grid - center) ** 2 / (2 * sigma ** 2))
    if normalize: g /= sigma * np.sqrt(2 * np.pi)
    return g