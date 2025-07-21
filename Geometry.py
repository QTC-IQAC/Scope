import numpy as np
import math

#########
def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def getangle(v1, v2):
    """
    Calculates the angle in radians between two vectors.

    Parameters:
        v1 (array-like): The first vector.
        v2 (array-like): The second vector.

    Returns:
        float: The angle in radians between vectors v1 and v2.

    Notes:
        - The function assumes that `unit_vector` and `np` (NumPy) are defined/imported in the scope.
        - The result is in the range [0, π].
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

#########
def getdihedral(V1, V2, V3, V4):
    """
    Calculate the dihedral angle (torsion angle) defined by four points in 3D space.

    Given four points V1, V2, V3, and V4, this function computes the signed dihedral angle
    between the planes formed by (V1, V2, V3) and (V2, V3, V4). The angle is returned in radians,
    ranging from -π to π.

    Parameters
    ----------
    V1, V2, V3, V4 : array-like
        The coordinates of the four points, each as a 1D array-like of length 3.

    Returns
    -------
    angle : float
        The signed dihedral angle in radians.

    Notes
    -----
    The sign of the angle follows the right-hand rule and is determined using the atan2 function.
    """
    V1, V2, V3, V4 = map(np.asarray, (V1, V2, V3, V4))
    # Bond vectors
    b1 = V2 - V1
    b2 = V3 - V2
    b3 = V4 - V3
    # Normalize b2 for projection
    b2_norm = b2 / np.linalg.norm(b2)
    # Normals to the planes
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    # Normalize normals
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    # Orthogonal vector to n1 in the plane of rotation
    m1 = np.cross(n1, b2_norm)
    # Compute angle using atan2 to get the sign
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    angle = np.arctan2(y, x)
    return angle

#########
# Rotation Matrices
#########
# 3D rotation matrix along x-axis
def rot_in_x(theta:float,input_geom:np.ndarray):
    c, s = math.cos(theta), math.sin(theta)
    rotx = np.array([[1.0, 0, 0], [0, c, -s], [0, s, c]])
    output_geom = np.round(np.dot(input_geom, rotx),8)
    return output_geom

# 3D rotation matrix along y-axis
def rot_in_y(theta:float,input_geom:np.ndarray):
    c, s = math.cos(theta), math.sin(theta)
    roty = np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])
    output_geom = np.round(np.dot(input_geom,roty),8)
    return output_geom

# 3D rotation matrix along z-axis
def rot_in_z(theta:float,input_geom:np.ndarray):
    c, s = math.cos(theta), math.sin(theta)
    rotz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
    output_geom = np.round(np.dot(input_geom, rotz),8)
    return output_geom

#########
def centercoords(coords, origin_atom):
    coords = np.array(coords)
    origin = coords[origin_atom]
    return coords - origin

#########
def displace_coords(coords, origin_atom, new_atom_coords):
    coords = np.array(coords)
    new_atom_coords = np.array(new_atom_coords)
    displacement_vector = new_atom_coords - coords[origin_atom]
    return coords + displacement_vector

#########
def set_dihedral(labels: list, coord: list, dih: float, atom1: int, atom2: int, atom3: int, atom4: int, debug: int=0):
    from copy import deepcopy
    from scipy import sparse
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee
    from typing import Tuple
    from Scope.Elementdata import ElementData
    from Scope.Adapted_from_cell2mol import get_adjmatrix, get_blocks, inv
    elemdatabase = ElementData()

    d0 = getdihedral(coord[atom1], coord[atom2], coord[atom3], coord[atom4])
    if debug > 0: print(f"initial dihedral {np.degrees(d0)=}")

    # centers coordinate to atom2
    c1 = centercoords(coord, atom2)
    if debug > 0: print(f"{c1=}")

    # we get rid of the X coord of atom3
    v1 = np.array(list([c1[atom3][0],c1[atom3][1],0])) 
    v2 = np.array([1,0,0])
    a1 = getangle(v1,v2) ## angle between xy projection of atom 3 and the x axis
    if debug > 0: print(f"{np.degrees(a1)=}")
    c2 = rot_in_z(a1, c1)
    if np.abs(c2[atom3][1]) > 0.0001: c2 = rot_in_z(-a1, c1)
    if debug > 0: print(f"{c2=}")
    # now we get rid of the Z coord of atom3
    v3 = np.array(list([c2[atom3][0],0,c2[atom3][2]]))
    v4 = np.array([1,0,0])
    a2 = getangle(v3,v4) ## angle between xz projection of atom 3 and the x axis
    if debug > 0: print(f"{a2=}")
    c3 = rot_in_y(a2, c2)
    if np.abs(c3[atom3][2]) > 0.0001: c3 = rot_in_y(-a2, c2)
    if debug > 0: print(f"{c3=}")
    
    ### Evaluates which atoms must be moved when applying dihedral. First, it removes connection between atoms2 and 3
    isgood, adjmat, adjnum = get_adjmatrix(labels, c3)
    adjnum[atom2] = adjnum[atom2]-1
    adjnum[atom3] = adjnum[atom3]-1
    adjmat[atom2][atom3] = 0
    adjmat[atom3][atom2] = 0

    ### This below is part of split_species function in Scope
    indices = [*range(0,len(labels),1)]
    degree = np.diag(adjnum)  # creates a matrix with adjnum as diagonal values. Needed for the laplacian
    lap = adjmat - degree

    # creates block matrix
    graph = csr_matrix(lap)
    perm = reverse_cuthill_mckee(graph)
    gp1 = graph[perm, :]
    gp2 = gp1[:, perm]
    dense = gp2.toarray()

    # detects blocks in the block diagonal matrix called "dense"
    startlist, endlist = get_blocks(dense)

    nblocks = len(startlist)
    if debug > 0: print("nblocks", nblocks)
    # keeps track of the atom movement within the matrix. Needed later
    atomlist = np.zeros((len(dense)))
    for b in range(0, nblocks):
        for i in range(0, len(dense)):
            if (i >= startlist[b]) and (i <= endlist[b]):
                atomlist[i] = b + 1
    invperm = inv(perm)
    atomlistperm = [int(atomlist[i]) for i in invperm]

    # assigns atoms to blocks
    blocklist = []
    for b in range(0, nblocks):
        atlist = []    # atom indices in the original ordering
        for i in range(0, len(atomlistperm)):
            if atomlistperm[i] == b + 1:
                atlist.append(indices[i])
        blocklist.append(atlist)
    for b in blocklist:
        if debug > 0: print("found block:", b, len(b))
        if atom3 in b and atom4 in b: atoms_to_move = deepcopy(b)
    
    if debug > 0: print("Atoms to move", atoms_to_move)
    
    d1 = getdihedral(c3[atom1], c3[atom2], c3[atom3], c3[atom4])
    if debug > 0: print(f"current dihedral {np.degrees(d1)=}")
    d2 = np.radians(dih) ## desired dihedral    
    if debug > 0: print(f"desired dihedral {np.degrees(d2)=}")

    d3 = d2 - d1 ## dihedral to be applied
    if debug > 0: print(f"rotation to apply {np.degrees(d3)=}")

    c4 = rot_in_x(d3, c3)  # dihedral is applied as a rotation along the x axis
    c5 = []
    for idx in range(len(labels)):
        if idx in atoms_to_move: c5.append(list(c4[idx]))          # takes rotated coordinates
        else:                    c5.append(list(c3[idx]))
    d4 = getdihedral(c5[atom1], c5[atom2], c5[atom3], c5[atom4])
    if debug > 0: print(f"current dihedral {np.degrees(d4)=}")

    if np.degrees(np.abs(d4 - d2)) > 5:
        if debug > 0: print(f"inverting rotation")
        c4 = rot_in_x(-d3, c3)  # invert rotation sign, as it means is counterclockwise

        c5 = []
        for idx in range(len(labels)):
            if idx in atoms_to_move: c5.append(list(c4[idx]))
            else:                    c5.append(list(c3[idx]))
        d4 = getdihedral(c5[atom1], c5[atom2], c5[atom3], c5[atom4])
    if debug > 0: print(f"final dihedral {np.degrees(d4)=}")

    c6 = displace_coords(c5, atom2, coord[atom2])
    return c6
