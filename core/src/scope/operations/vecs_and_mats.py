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

######
def gcd(a, b):
    # gcd = Greatest common divisor
    while b != 0:
        a, b = b, a % b
    return a

######
def gcd_list(numbers):
    # returns the gcd for a list of numbers
    result = numbers[0]
    for n in numbers[1:]:
        result = gcd(result, n)
    return result

###############
## Functions ##
###############
def gaussian(grid, center, sigma=0.2, normalize=False, debug: int=0):
    """
    Gaussian profile centered at `center`, and evaluated over 'grid'.
    If normalize, it normalizes the final profile to unit area. Default is False.
    """
    sigma = _resolve_sigma(sigma, grid)
    profile = np.exp(-(grid - center) ** 2 / (2 * sigma ** 2))
    if normalize: profile /= sigma * np.sqrt(2 * np.pi)
    return profile

######
def lorentzian(grid, center, sigma=0.2, normalize=False, debug: int=0):
    """
    Lorentzian profile centered at `center`, and evaluated over 'grid'.
    If normalize, it normalizes the final profile to unit area. Default is False.
    """
    sigma = _resolve_sigma(sigma, grid)
    profile = 1 / (1 + ((grid - center) / sigma) ** 2)
    if normalize: profile /= np.pi * sigma
    return profile

######
def laplacian(grid, center, sigma=0.2, normalize=False, debug: int=0):
    """
    Laplacian profile centered at `center`, and evaluated over 'grid'.
    If normalize, it normalizes the final profile to unit area. Default is False.
    """
    sigma = _resolve_sigma(sigma, grid)
    profile = np.exp(-np.abs(grid - center) / sigma)
    if normalize: profile /= 2 * sigma
    return profile

#################
## Convolution ##
#################
def _auto_sigma(values, default=0.2):
    """
    Estimate sigma as the median positive spacing of sorted input values.
    """
    vals = np.asarray(values, dtype=float).ravel()
    if vals.size < 2:
        return float(default)

    diffs = np.diff(np.sort(vals))
    diffs = diffs[diffs > 0]
    if diffs.size == 0:
        return float(default)
    return float(np.median(diffs))

######
def _resolve_sigma(sigma, values, default=0.2):
    """
    Resolve sigma from explicit value or auto mode.
    """
    if sigma is None or (isinstance(sigma, str) and sigma.lower() == "auto"):
        sigma = _auto_sigma(values, default=default)
    sigma = float(sigma)
    if sigma <= 0:
        raise ValueError("sigma must be > 0")
    return sigma

######
def build_spectrum(xrange, xvalues, yvalues, function: str='gaussian', sigma=0.2, normalize: bool=False, debug: int=0):
    """
    Build a broadened spectrum from discrete signal points.
    
    Parameters
    ----------
    xrange : array_like            Grid where the broadened signal is evaluated.
    xvalues : list of floats       Discrete x positions for each transition/peak.
    yvalues : list of floats       Discrete y intensities (weights) for each transition/peak.
    function : str, optional       Broadening function, either "gaussian", "lorentzian" or "laplacian". Default is "gaussian".
    sigma : float | str | None     Width parameter of the broadening function. Use "auto" (or None) to estimate sigma from median spacing in xvalues.
    normalize : bool, optional     Whether to normalize each kernel to unit area. Default is False.
    Returns
    -------
    tuple(array_like, array_like)  x grid and broadened spectrum.
    """
    assert len(xvalues) == len(yvalues), "BUILD_SPECTRUM: xvalues and yvalues must have the same length."
    
    sigma_resolved = _resolve_sigma(sigma, xvalues if len(xvalues) > 1 else xrange)
    if debug > 0 and (sigma is None or (isinstance(sigma, str) and sigma.lower() == "auto")):
        print(f'BUILD_SPECTRUM: auto sigma = {sigma_resolved}')

    if debug > 0: print(f"BUILD_SPECTRUM: Using {sigma_resolved=}")
    func = function.lower()
    spec = np.zeros_like(xrange, dtype=float)
    for x0, y0 in zip(xvalues, yvalues):
        if debug > 0: print(f'BUILD_SPECTRUM: x0: {x0} y0: {y0}')
        if   func.lower() == 'gaussian':     spec += y0 * gaussian(xrange, x0, sigma=sigma_resolved, normalize=normalize)
        elif func.lower() == 'laplacian':    spec += y0 * laplacian(xrange, x0, sigma=sigma_resolved, normalize=normalize)
        elif func.lower() == 'lorentzian':   spec += y0 * lorentzian(xrange, x0, sigma=sigma_resolved, normalize=normalize)
        else:                       raise ValueError("BUILD_SPECTRUM: the kernel function must be 'gaussian', 'lorentzian' or 'laplacian'")
    return xrange, spec
