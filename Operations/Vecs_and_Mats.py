def normalize(mat):
    norms = np.linalg.norm(mat, axis=1)
    norms[norms == 0] = 1.0
    return mat / norms[:, None]
