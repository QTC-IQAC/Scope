import numpy as np
import math
from scope.operations.vecs_and_mats import normalize, determinant

#################
## Translation ##
#################
def centercoords(coords: list, atom_idx: int):
    """
    Translate coordinates so a selected atom lies at the origin.

    Parameters:
        coords (list):                  Cartesian coordinates.
        atom_idx (int):                 Index of the reference atom.

    Returns:
        ndarray: Shifted coordinates.
    """
    coords = np.array(coords)
    origin = coords[atom_idx]
    return coords - origin

def displace_coords(coords: list, atom_idx: int, point: list=[0,0,0]):
    """
    Translate coordinates so a selected atom reaches a target point.

    Parameters:
        coords (list):                  Cartesian coordinates.
        atom_idx (int):                 Index of the atom to move.
        point (list):                   Target point.

    Returns:
        ndarray: Shifted coordinates.
    """
    coords = np.array(coords)
    point  = np.array(point)
    disp   = point - coords[atom_idx] # Displacement vector
    return coords + disp

##############
## Distance ##
##############
def get_dist(coord1: list, coord2: list) -> float:
    dist = np.linalg.norm(np.array(coord1) - np.array(coord2))
    return float(dist)

###########
## Angle ##
###########
def get_angle_vectors(v1, v2, eps: float=1e-8) -> float:
    """
    Calculate the angle in radians between two vectors.

    Parameters:
        v1 (array-like):                First vector.
        v2 (array-like):                Second vector.
        eps (float):                    Small cutoff for numerical stability.

    Returns:
        float: Angle in radians.
    """
    v1 = normalize(v1)
    v2 = normalize(v2)
    ## To avoid numerical instabilities after normalization
    if np.linalg.norm(v1) < eps or np.linalg.norm(v2) < eps: return 0.0 
    return float(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))

def get_angle_points(P1, P2, P3, eps: float=1e-8) -> float:
    """
    Calculate the angle in radians defined by three points.

    Parameters:
        P1 (array-like):                First point.
        P2 (array-like):                Vertex point.
        P3 (array-like):                Third point.
        eps (float):                    Small cutoff for numerical stability.

    Returns:
        float: Angle in radians.
    """
    P1, P2, P3 = map(np.asarray, (P1, P2, P3))
    v1 = P1 - P2
    v2 = P3 - P2
    return get_angle_vectors(v1, v2, eps=eps)

def get_angle(a1, a2, a3=None, eps: float=1e-8) -> float:
    """
    Compute an angle from either two vectors or three points.

    Parameters:
        a1:                            First vector or point.
        a2:                            Second vector or vertex point.
        a3:                            Optional third point.
        eps (float):                   Small cutoff for numerical stability.

    Returns:
        float: Angle in radians.
    """
    if a3 is None:
        return get_angle_vectors(a1, a2, eps=eps)
    return get_angle_points(a1, a2, a3, eps=eps)

#########
def set_angle(labels, coord, target_angle: float, atom1: int, atom2: int, atom3: int, debug=0):
    from copy import deepcopy
    from scipy.sparse         import csr_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee
    from scope.connectivity   import get_adjmatrix, get_blocks, inv

    """
    Adjust the angle formed by three atoms to a target value.

    Parameters:
        labels:                        Atomic symbols.
        coord:                         Cartesian coordinates.
        target_angle (float):          Target angle in degrees.
        atom1 (int):                   First atom index.
        atom2 (int):                   Vertex atom index.
        atom3 (int):                   Third atom index.
        debug (int):                   Verbosity level.

    Returns:
        list: Updated coordinates.
    """
    xy_plane = put_atoms_on_xy(coord, atom1, atom2, atom3, debug)
    isgood, adjmat, adjnum = get_adjmatrix(labels, xy_plane)
    
    if not isgood:
        print("SET_ANGLE: Adjacency matrix is not good, returning original coordinates")
        return coord

    # Remove connection between atoms 2 and 3
    adjnum[atom2] = adjnum[atom2]-1
    adjnum[atom3] = adjnum[atom3]-1
    adjmat[atom2][atom3] = 0
    adjmat[atom3][atom2] = 0

    ### Split species 
    indices = [*range(0, len(labels), 1)]
    degree = np.diag(adjnum)
    lap = adjmat - degree
    graph = csr_matrix(lap)
    perm = reverse_cuthill_mckee(graph)
    gp1 = graph[perm, :]
    gp2 = gp1[:, perm]
    dense = gp2.toarray()
    
    startlist, endlist = get_blocks(dense)
    nblocks = len(startlist)

    if debug > 0: print("SET_ANGLE: nblocks: ", nblocks)
    
    atomlist = np.zeros((len(dense)))
    for b in range(0, nblocks):
        for i in range(0, len(dense)):
            if (i >= startlist[b]) and (i <= endlist[b]):
                atomlist[i] = b + 1
                
    invperm = inv(perm)
    atomlistperm = [int(atomlist[i]) for i in invperm]

    blocklist = []
    for b in range(0, nblocks):
        atlist = []
        for i in range(0, len(atomlistperm)):
            if atomlistperm[i] == b + 1:
                atlist.append(indices[i])
        blocklist.append(atlist)
        
    for b in blocklist:
        if debug > 0: print("SET_ANGLE: found block: ", b, len(b))
        if atom3 in b: atoms_to_move = deepcopy(b)

    if debug > 0: print("SET_ANGLE: atoms_to_move: ", atoms_to_move)

    v4 = xy_plane[atom1] - xy_plane[atom2]
    v5 = xy_plane[atom3] - xy_plane[atom2]
    target_angle_rad = target_angle * (np.pi / 180)
    a5 = get_angle(v4, v5)
    if debug > 0: print("SET_ANGLE: target_angle: ", target_angle)

    theta = a5 - target_angle_rad
    c4 = []
    for idx in range(len(xy_plane)):
        if idx in atoms_to_move:
            c4.append(xy_plane[idx])

    c5 = rot_in_z(theta, c4)
    coord_final = deepcopy(xy_plane)
    for i, idx in enumerate(atoms_to_move):
        if idx in atoms_to_move:
            coord_final[idx] = c5[i]
    return coord_final

####################
## Dihedral Angle ##
####################
def get_dihedral(P1, P2, P3, P4, eps: float=1e-8) -> float:
    """
    Calculate the signed dihedral angle defined by four points.

    Parameters:
        P1:                            First point.
        P2:                            Second point.
        P3:                            Third point.
        P4:                            Fourth point.
        eps (float):                   Small cutoff for numerical stability.

    Returns:
        float: Signed dihedral angle in radians.
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
    n2 = normalize(n2)
    ## To avoid numerical instabilities after normalization, or cases of collinearity:
    if np.linalg.norm(n1) < eps or np.linalg.norm(n2) < eps or np.linalg.norm(b2_norm) < eps: return 0.0 
    # Orthogonal vector to n1 in the plane of rotation
    m1 = np.cross(n1, b2_norm)
    # Compute angle using atan2 to get the sign
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    angle = float(np.arctan2(y, x))
    return angle       # In radians

#########
def set_dihedral(labels: list, coord: list, dih: float, atom1: int, atom2: int, atom3: int, atom4: int, adjmat=None, adjnum=None,  debug: int=0):
    from copy import deepcopy
    from scipy.sparse         import csr_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee
    from scope.connectivity   import get_adjmatrix, get_blocks, inv

    d0 = get_dihedral(coord[atom1], coord[atom2], coord[atom3], coord[atom4])
    if debug > 0: print(f"SET_DIHEDRAL: initial dihedral {np.degrees(d0)=}")

    ### Evaluates which atoms must be moved later, when applying dihedral. 
    # First, it commputes the adjacency matrix
    if adjmat is None or adjnum is None:
        isgood, adjmat, adjnum = get_adjmatrix(labels, coord)
    else:
        adjmat = adjmat.copy()
        adjnum = adjnum.copy()
    # Second, it removes connection between atoms2 and 3
    adjnum[atom2] = adjnum[atom2]-1
    adjnum[atom3] = adjnum[atom3]-1
    adjmat[atom2][atom3] = 0
    adjmat[atom3][atom2] = 0

    # centers coordinate to atom2
    c1 = centercoords(coord, atom2)
    if debug > 0: print(f"{c1=}")

    # we get rid of the X coord of atom3
    v1 = np.array(list([c1[atom3][0],c1[atom3][1],0])) 
    v2 = np.array([1,0,0])
    a1 = get_angle(v1,v2) ## angle between xy projection of atom 3 and the x axis
    if debug > 0: print(f"SET_DIHEDRAL: {np.degrees(a1)=}")
    c2 = rot_in_z(a1, c1)
    if np.abs(c2[atom3][1]) > 0.0001: c2 = rot_in_z(-a1, c1)
    if debug > 0: print(f"SET_DIHEDRAL: {c2=}")
    # now we get rid of the Z coord of atom3
    v3 = np.array(list([c2[atom3][0],0,c2[atom3][2]]))
    v4 = np.array([1,0,0])
    a2 = get_angle(v3,v4) ## angle between xz projection of atom 3 and the x axis
    if debug > 0: print(f"SET_DIHEDRAL: {a2=}")
    c3 = rot_in_y(a2, c2)
    if np.abs(c3[atom3][2]) > 0.0001: c3 = rot_in_y(-a2, c2)
    if debug > 0: print(f"SET_DIHEDRAL: {c3=}")

    ### This below is part of split_species function in scope
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
        if debug > 0: print("SET_DIHEDRAL: found block:", b, len(b))
        if atom3 in b and atom4 in b: atoms_to_move = deepcopy(b)
    
    if debug > 0: print("SET_DIHEDRAL: Atoms to move", atoms_to_move)
    
    d1 = get_dihedral(c3[atom1], c3[atom2], c3[atom3], c3[atom4])
    if debug > 0: print(f"SET_DIHEDRAL: current dihedral {np.degrees(d1)=}")
    d2 = np.radians(dih) ## desired dihedral    
    if debug > 0: print(f"SET_DIHEDRAL: desired dihedral {np.degrees(d2)=}")

    d3 = d2 - d1 ## dihedral to be applied
    if debug > 0: print(f"SET_DIHEDRAL: rotation to apply {np.degrees(d3)=}")

    c4 = rot_in_x(d3, c3)  # dihedral is applied as a rotation along the x axis
    c5 = []
    for idx in range(len(labels)):
        if idx in atoms_to_move: c5.append(list(c4[idx]))          # takes rotated coordinates
        else:                    c5.append(list(c3[idx]))
    d4 = get_dihedral(c5[atom1], c5[atom2], c5[atom3], c5[atom4])
    if debug > 0: print(f"SET_DIHEDRAL: current dihedral {np.degrees(d4)=}")

    if np.degrees(np.abs(d4 - d2)) > 5:
        if debug > 0: print(f"SET_DIHEDRAL: inverting rotation")
        c4 = rot_in_x(-d3, c3)  # invert rotation sign, as it means is counterclockwise

        c5 = []
        for idx in range(len(labels)):
            if idx in atoms_to_move: c5.append(list(c4[idx]))
            else:                    c5.append(list(c3[idx]))
        d4 = get_dihedral(c5[atom1], c5[atom2], c5[atom3], c5[atom4])
    if debug > 0: print(f"SET_DIHEDRAL: final dihedral {np.degrees(d4)=}")

    c6 = displace_coords(c5, atom2, coord[atom2])
    return c6

######
def solve_dihedral(labels: list, coord: list, at0: int, at1: int, at2: int, at3: int, at4: int, at5: int, adjmat_ref, adjnum_ref, debug: int=0):
    """
    Adjust adjacent dihedrals to resolve steric clashes after a torsion change.

    Parameters:
        labels (list):                  Atomic symbols.
        coord (list):                   Cartesian coordinates.
        at0, at1, at2, at3, at4, at5:  Atom indices defining the main and adjacent dihedrals.
        adjmat_ref:                     Reference adjacency matrix.
        adjnum_ref:                     Reference coordination numbers.
        debug (int):                    Verbosity level.

    Returns:
        tuple: `(worked, coordinates)` with the best geometry found.
    """
    from copy import deepcopy
    from itertools import chain
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee
    from scope.connectivity import get_adjmatrix, get_blocks

    def get_nblocks_for_adjmat(adjmat, adjnum):
        degree  = np.diag(adjnum)
        lap     = adjmat - degree
        graph   = csr_matrix(lap)
        perm    = reverse_cuthill_mckee(graph)
        dense   = graph[perm, :][:, perm].toarray()
        startlist, endlist = get_blocks(dense)
        return len(startlist)

    def can_rotate_dihedral(adjmat, adjnum, bond_a, bond_b):
        adjmat_try = deepcopy(adjmat)
        adjnum_try = deepcopy(adjnum)

        if adjmat_try[bond_a][bond_b] == 0 or adjmat_try[bond_b][bond_a] == 0:
            return False
        adjmat_try[bond_a][bond_b] = 0
        adjmat_try[bond_b][bond_a] = 0
        adjnum_try[bond_a] -= 1
        adjnum_try[bond_b] -= 1
        return get_nblocks_for_adjmat(adjmat_try, adjnum_try) == 2

    def can_rotate_adjacent_dihedrals(adjmat, adjnum):
        adjmat_try = deepcopy(adjmat)
        adjnum_try = deepcopy(adjnum)
        for bond_a, bond_b in ((at1, at2), (at3, at4)):
            if adjmat_try[bond_a][bond_b] == 0 or adjmat_try[bond_b][bond_a] == 0:
                return False
            adjmat_try[bond_a][bond_b] = 0
            adjmat_try[bond_b][bond_a] = 0
            adjnum_try[bond_a] -= 1
            adjnum_try[bond_b] -= 1
        return get_nblocks_for_adjmat(adjmat_try, adjnum_try) == 3

    def get_rot_combinations(rot_steps_left, rot_steps_right):
        nleft   = len(rot_steps_left)
        nright  = len(rot_steps_right)
        nlayers = max(nleft, nright)

        for layer in range(nlayers):
            if layer < nright:
                yield rot_steps_left[0], rot_steps_right[layer]
            if layer > 0 and layer < nleft:
                yield rot_steps_left[layer], rot_steps_right[0]
            for idx in range(1, layer):
                if idx < nleft and layer < nright:
                    yield rot_steps_left[idx], rot_steps_right[layer]
                if layer < nleft and idx < nright:
                    yield rot_steps_left[layer], rot_steps_right[idx]
            if layer > 0 and layer < nleft and layer < nright:
                yield rot_steps_left[layer], rot_steps_right[layer]

    ## Computes left and right dihedral angles
    rot_steps = np.linspace(-180,180,64).astype(int)
    left_dih  = np.degrees(get_dihedral(coord[at0], coord[at1], coord[at2], coord[at3]))
    right_dih = np.degrees(get_dihedral(coord[at2], coord[at3], coord[at4], coord[at5]))
    if debug > 0: print(f'SOLVE_DIHEDRAL: Original adjacent dihedral angles {left_dih=}, {right_dih=}')

    rot_steps_ordered_left  = rot_steps[np.argsort(np.abs(((rot_steps - left_dih + 180) % 360) - 180))]
    rot_steps_ordered_right = rot_steps[np.argsort(np.abs(((rot_steps - right_dih + 180) % 360) - 180))]

    can_rotate_left  = can_rotate_dihedral(adjmat_ref, adjnum_ref, at1, at2)
    can_rotate_right = can_rotate_dihedral(adjmat_ref, adjnum_ref, at3, at4)
    can_rotate_both  = can_rotate_adjacent_dihedrals(adjmat_ref, adjnum_ref)
    if debug > 0: print(f'SOLVE_DIHEDRAL: {can_rotate_left=}, {can_rotate_right=}, {can_rotate_both=}')

    if can_rotate_left and can_rotate_right and can_rotate_both:
        rot_combinations = get_rot_combinations(rot_steps_ordered_left, rot_steps_ordered_right)
    elif can_rotate_left and can_rotate_right:
        rot_combinations = chain(
            get_rot_combinations(rot_steps_ordered_left, np.array([rot_steps_ordered_right[0]])),
            get_rot_combinations(np.array([rot_steps_ordered_left[0]]), rot_steps_ordered_right))
    elif can_rotate_left:
        rot_combinations = get_rot_combinations(rot_steps_ordered_left, np.array([rot_steps_ordered_right[0]]))
    elif can_rotate_right:
        rot_combinations = get_rot_combinations(np.array([rot_steps_ordered_left[0]]), rot_steps_ordered_right)
    else:
        print("SOLVE_DIHEDRAL: Adjacent dihedrals are not rotatable. Returning original coordinates.")
        return False, coord

    worked = False
    for angle1, angle2 in rot_combinations:
        if debug>0: print(f'SOLVE_DIHEDRAL: Trying {angle1=} {angle2=}')
        new_coord = set_dihedral(labels, coord, angle1, at0,at1,at2,at3, adjmat=adjmat_ref, adjnum=adjnum_ref)      # Sending original coords
        new_coord = set_dihedral(labels, new_coord, angle2, at2,at3,at4,at5, adjmat=adjmat_ref, adjnum=adjnum_ref)  # Sending coords after first modification above

        # Gets the adjacency matrix of the new coordinates, to check wether the original connectivity is preserved (i.e., no steric clashes)
        _, adjmat_try, adjnum_try = get_adjmatrix(labels,new_coord)
        is_same = np.array_equal(adjmat_try, adjmat_ref) and np.array_equal(adjnum_try, adjnum_ref)
        if is_same:
            if debug != 0: print(f'SOLVE_DIHEDRAL: Found good geometry by rotating adjacent dihedrals to {angle1=} {angle2=}')
            worked = True
            return worked, new_coord
    print(f"SOLVE_DIHEDRAL: Couldn't find good geometry. Returning original coordinates.")
    return worked, coord

#########
def get_planar_distortion(theta: float) -> float:
    """
    Compute the minimum angular distance to a planar orientation.

    Parameters:
        theta (float | array-like):     Input angle in radians.

    Returns:
        float: Minimum angular distance in radians.
    """
    mod_theta = np.mod(theta, np.pi)
    return float(np.minimum(mod_theta, np.pi - mod_theta))

#######################
## Rotation Matrices ##
#######################
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
def put_atoms_on_xy(coord: list, at1: int, at2: int, at3: int, debug=0):
    """
    Rotate coordinates so one bond lies on the x axis and a third atom lies on the xy plane.

    Parameters:
        coord (list):                   Cartesian coordinates.
        at1 (int):                      First atom used to define the bond.
        at2 (int):                      Central atom moved to the origin.
        at3 (int):                      Third atom used to define the plane.
        debug (int):                    Verbosity level.

    Returns:
        ndarray: Rotated coordinates.
    """
    # 1. Center geometry at at2 at (0,0,0)
    c1 = centercoords(coord, at2) 
    
    # 2. Rotate around Z-axis to eliminate the Y component of at1
    v1 = c1[at1]
    theta_z = math.atan2(v1[1], v1[0]) 
    c2 = rot_in_z(theta_z, c1)
    
    # 3. Rotate around Y-axis to eliminate the Z component of at1
    v1_z = c2[at1]
    theta_y = math.atan2(-v1_z[2], v1_z[0])
    c3 = rot_in_y(theta_y, c2)
    
    # 4. Rotate around X-axis to eliminate the Z component of at3
    v3_y = c3[at3]
    theta_x = math.atan2(v3_y[2], v3_y[1])
    c4 = rot_in_x(theta_x, c3)

    if debug > 0: 
        print(f"PUT_ATOMS_ON_XY: Applied rotations (degrees) -> Z: {np.degrees(theta_z):.2f}, Y: {np.degrees(theta_y):.2f}, X: {np.degrees(theta_x):.2f}")
    return c4

###############
## Unit Cell ##
###############
def get_unit_cell_volume(a: float, b: float, c: float, alpha: float, beta: float, gamma: float):
    # I know alpha and gamma are not used, but that way it is just simpler to execute as:
    # get_unit_cell_volume(*cell_parameters)
    return float(a*b*c*np.sin(np.deg2rad(beta)))

######
def cellparam_2_cellvec(*args) -> list:
    """
    Convert unit-cell parameters into Cartesian cell vectors.

    Parameters:
        *args:                          Either six scalar values or one 6-item sequence.

    Returns:
        list: Three Cartesian cell vectors.
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

######
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
