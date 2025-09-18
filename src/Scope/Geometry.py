import numpy as np
import math
from Scope.Operations.Vecs_and_Mats import normalize, determinant

######
def get_dist(coord1: list, coord2: list) -> float:
    dist = np.linalg.norm(np.array(coord1) - np.array(coord2))
    return dist

#########
def get_angle(v1, v2) -> float:
    """
    Calculates the angle in radians between two vectors.

    Parameters:
        v1 (array-like): The first vector.
        v2 (array-like): The second vector.

    Returns:
        float: The angle in radians between vectors v1 and v2.

    Notes:
        - The result is in the range [0, π].
    """
    v1_u = normalize(v1)
    v2_u = normalize(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

#########
def get_dihedral(P1, P2, P3, P4) -> float:
    """
    Calculate the dihedral angle (torsion angle) defined by four points in 3D space.

    Given four points P1, P2, P3, and P4, this function computes the signed dihedral angle
    between the planes formed by (P1, P2, P3) and (P2, P3, P4). The angle is returned in radians,
    ranging from -π to π.

    Parameters
    ----------
    P1, P2, P3, P4 : array-like
        The coordinates of the four points, each as a 1D array-like of length 3.

    Returns
    -------
    angle : float
        The signed dihedral angle in radians.

    Notes
    -----
    The sign of the angle follows the right-hand rule and is determined using the atan2 function.
    """
    P1, P2, P3, P4 = map(np.asarray, (P1, P2, P3, P4))
    # Bond vectors
    b1 = P2 - P1
    b2 = P3 - P2
    b3 = P4 - P3
    # Normalize b2 for projection
    b2_norm = normalize(b2) 
    # Normals to the planes
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    # Normalize normals
    n1 = normalize(n1)
    n2 = normalize(n1)
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
    from .Elementdata import ElementData
    from .Connectivity import get_adjmatrix, get_blocks, inv
    elemdatabase = ElementData()

    d0 = get_dihedral(coord[atom1], coord[atom2], coord[atom3], coord[atom4])
    if debug > 0: print(f"initial dihedral {np.degrees(d0)=}")

    # centers coordinate to atom2
    c1 = centercoords(coord, atom2)
    if debug > 0: print(f"{c1=}")

    # we get rid of the X coord of atom3
    v1 = np.array(list([c1[atom3][0],c1[atom3][1],0])) 
    v2 = np.array([1,0,0])
    a1 = get_angle(v1,v2) ## angle between xy projection of atom 3 and the x axis
    if debug > 0: print(f"{np.degrees(a1)=}")
    c2 = rot_in_z(a1, c1)
    if np.abs(c2[atom3][1]) > 0.0001: c2 = rot_in_z(-a1, c1)
    if debug > 0: print(f"{c2=}")
    # now we get rid of the Z coord of atom3
    v3 = np.array(list([c2[atom3][0],0,c2[atom3][2]]))
    v4 = np.array([1,0,0])
    a2 = get_angle(v3,v4) ## angle between xz projection of atom 3 and the x axis
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
    
    d1 = get_dihedral(c3[atom1], c3[atom2], c3[atom3], c3[atom4])
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
    d4 = get_dihedral(c5[atom1], c5[atom2], c5[atom3], c5[atom4])
    if debug > 0: print(f"current dihedral {np.degrees(d4)=}")

    if np.degrees(np.abs(d4 - d2)) > 5:
        if debug > 0: print(f"inverting rotation")
        c4 = rot_in_x(-d3, c3)  # invert rotation sign, as it means is counterclockwise

        c5 = []
        for idx in range(len(labels)):
            if idx in atoms_to_move: c5.append(list(c4[idx]))
            else:                    c5.append(list(c3[idx]))
        d4 = get_dihedral(c5[atom1], c5[atom2], c5[atom3], c5[atom4])
    if debug > 0: print(f"final dihedral {np.degrees(d4)=}")

    c6 = displace_coords(c5, atom2, coord[atom2])
    return c6

##############################
### Former Unit Cell Tools ###
##############################
def get_unit_cell_volume(a, b, c, alpha, beta, gamma):
    return float(a*b*c*np.sin(np.deg2rad(beta)))

######
def cellparam_2_cellvec(*args) -> list:
    """
    Convert unit cell parameters into three Cartesian cell vectors.

    Accepts either:
      - six separate arguments: (a, b, c, alpha, beta, gamma)
      - a single list/tuple with six elements.

    Returns
    -------
    vectors : list of np.ndarray
        The three cell vectors [v1, v2, v3].
    """
    # Handle list or tuple input
    if len(args) == 1 and isinstance(args[0], (list, tuple, np.ndarray)):
        a, b, c, alpha, beta, gamma = args[0]
    else:
        a, b, c, alpha, beta, gamma = args

    # Convert angles to radians
    alpha = np.radians(alpha)
    beta  = np.radians(beta)
    gamma = np.radians(gamma)

    # Vectors
    v1 = np.array([a, 0.0, 0.0])
    v2 = np.array([b * np.cos(gamma), b * np.sin(gamma), 0.0])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = c * np.sqrt(1 - np.cos(beta)**2 - ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma))**2)
    v3 = np.array([cx, cy, cz])
    return list([v1, v2, v3])

def cellvec_2_cellparam(cellvec):
    a = np.linalg.norm(cellvec[0])
    b = np.linalg.norm(cellvec[1])
    c = np.linalg.norm(cellvec[2])
    
    tmp1 = cellvec[0][0]*cellvec[1][0]+cellvec[0][1]*cellvec[1][1]+cellvec[0][2]*cellvec[1][2]
    tmp2 = cellvec[0][0]*cellvec[2][0]+cellvec[0][1]*cellvec[2][1]+cellvec[0][2]*cellvec[2][2]
    tmp3 = cellvec[1][0]*cellvec[2][0]+cellvec[1][1]*cellvec[2][1]+cellvec[1][2]*cellvec[2][2]    
    tmp4 = a * b
    tmp5 = a * c
    tmp6 = b * c
    tmp7 = tmp1/tmp4
    tmp8 = tmp2/tmp5
    tmp9 = tmp3/tmp6
    
    alpha = np.arccos(tmp7)/(2*np.pi/360)
    beta = np.arccos(tmp8)/(2*np.pi/360)
    gamma = np.arccos(tmp9)/(2*np.pi/360)
    
    return list([a, b, c, alpha, beta, gamma])

######
def frac2cart_fromparam(frac_coord, cellparam):

    a = cellparam[0]
    b = cellparam[1]
    c = cellparam[2]
    alpha = np.radians(cellparam[3])
    beta = np.radians(cellparam[4])
    gamma = np.radians(cellparam[5])

    volume = get_unit_cell_volume(a, b, c, alpha, beta, gamma)

    m = np.zeros((3, 3))
    m[0][0] = a
    m[0][1] = b * np.cos(gamma)
    m[0][2] = c * np.cos(beta)
    m[1][0] = 0
    m[1][1] = b * np.sin(gamma)
    m[1][2] = c * ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma))
    m[2][0] = 0
    m[2][1] = 0
    m[2][2] = volume / (a * b * np.sin(gamma))

    cartesian = []
    for idx, frac in enumerate(frac_coord):
        xcar = frac[0] * m[0][0] + frac[1] * m[0][1] + frac[2] * m[0][2]
        ycar = frac[0] * m[1][0] + frac[1] * m[1][1] + frac[2] * m[1][2]
        zcar = frac[0] * m[2][0] + frac[1] * m[2][1] + frac[2] * m[2][2]
        cartesian.append([float(xcar), float(ycar), float(zcar)])
    return cartesian

######
def frac2cart_fromcellvec(frac_coord, cellvec):
    cartesian = []
    for idx, frac in enumerate(frac_coord):
        xcar = (frac[0] * cellvec[0][0] + frac[1] * cellvec[1][0] + frac[2] * cellvec[2][0])
        ycar = (frac[0] * cellvec[0][1] + frac[1] * cellvec[1][1] + frac[2] * cellvec[2][1])
        zcar = (frac[0] * cellvec[0][2] + frac[1] * cellvec[1][2] + frac[2] * cellvec[2][2])
        cartesian.append([float(xcar), float(ycar), float(zcar)])
    return cartesian

######
def cart2frac(cartCoords, cellvec):
    latCnt = [x[:] for x in [[None] * 3] * 3]
    for a in range(3):
        for b in range(3):
            latCnt[a][b] = cellvec[b][a]
    fracCoords = []
    detLatCnt = determinant(latCnt)
    for i in cartCoords:
        aPos = (determinant([[i[0], latCnt[0][1], latCnt[0][2]],[i[1], latCnt[1][1], latCnt[1][2]],[i[2], latCnt[2][1], latCnt[2][2]]])) / detLatCnt
        bPos = (determinant([[latCnt[0][0], i[0], latCnt[0][2]],[latCnt[1][0], i[1], latCnt[1][2]],[latCnt[2][0], i[2], latCnt[2][2]]])) / detLatCnt
        cPos = (determinant([[latCnt[0][0], latCnt[0][1], i[0]],[latCnt[1][0], latCnt[1][1], i[1]],[latCnt[2][0], latCnt[2][1], i[2]]])) / detLatCnt
        fracCoords.append([aPos, bPos, cPos])
    return fracCoords

######
def translate(vector, coords, cellvec):
    newcoord = []
    for idx, coord in enumerate(coords):
        newx = (
            coord[0]
            + vector[0] * cellvec[0][0]
            + vector[1] * cellvec[1][0]
            + vector[2] * cellvec[2][0]
        )
        newy = (
            coord[1]
            + vector[0] * cellvec[0][1]
            + vector[1] * cellvec[1][1]
            + vector[2] * cellvec[2][1]
        )
        newz = (
            coord[2]
            + vector[0] * cellvec[0][2]
            + vector[1] * cellvec[1][2]
            + vector[2] * cellvec[2][2]
        )
        newcoord.append([float(newx), float(newy), float(newz)])
    return newcoord
