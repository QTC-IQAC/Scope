import numpy as np
from scope import constants

####
def displace_neg_freqs(ini_coord, VNMs: object, debug: int=0) -> list:
    neg_VNMs = list([vnm for vnm in VNMs if vnm.freq_cm < 0.0])
    disp_coord = displace_coords_with_vnm(neg_VNMs, ini_coord)
    return disp_coord

####
def map_vnms(vnmsA,vnmsB,labelsA,coordsA,labelsB,coordsB,f_weight: float=1.00,f_scale: float=100.0,center_method: str="centroid",use_ext_info: bool=True,adjmatA=None,adjmatB=None,bond_ordersA=None,bond_ordersB=None,max_graph_mappings: int=1000,debug: int=0):
    """
    Align two molecular geometries and map their vibrational normal modes.

    Geometry A defines the reference atom ordering and Cartesian frame.
    Geometry B and its modes are reordered and rotated into that frame before
    calculating mass-weighted eigenvector overlaps.

    The assignment cost combines eigenvector disagreement and frequency
    separation:

        cost = (1 - overlap) + f_weight * (delta_freq/f_scale)^2

    Parameters:
        vnmsA, vnmsB (list):            VNM collections to map.
        labelsA, labelsB (list):        Atomic labels in geometry order.
        coordsA, coordsB (array):       Cartesian coordinates with shape `(natoms, 3)`.
        f_weight (float):               Strength of the frequency penalty.
        f_scale (float):                Characteristic frequency difference in cm-1.
        center_method (str):            Centering method used for molecular alignment.
        use_ext_info (bool):            Whether extended chemical data are used for atom mapping.
        adjmatA, adjmatB:               Optional adjacency matrices.
        bond_ordersA, bond_ordersB:     Optional bond-order matrices.
        max_graph_mappings (int):       Maximum graph mappings checked during alignment.
        debug (int):                    Verbosity level.

    Returns:
        tuple:
            vnmsA_ordered (list):       Modes from A ordered by their assignment rows.
            vnmsB_ordered (list):       Aligned copies of the matched modes from B.
            results (dict):             Mapping, alignment, and assignment matrices.
    """
    from copy import deepcopy
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    from scope.overlap  import overlap_molecules

    # 0) Checks and arranges arrays
    if len(vnmsA) == 0 or len(vnmsB) == 0: raise ValueError("MAP_VNMS: VNM collections cannot be empty")
    if len(vnmsA) != len(vnmsB):           raise ValueError("MAP_VNMS: VNM collections must contain the same number of modes")
    if f_weight < 0.0:                     raise ValueError("MAP_VNMS: f_weight cannot be negative")
    if f_scale <= 0.0:                     raise ValueError("MAP_VNMS: f_scale must be positive")
    if any(not vnm.has_mode for vnm in [*vnmsA, *vnmsB]):  raise ValueError("MAP_VNMS: All VNMs must contain eigenvectors")

    labelsA = np.asarray(labelsA)
    labelsB = np.asarray(labelsB)
    coordsA = np.asarray(coordsA, dtype=float)
    coordsB = np.asarray(coordsB, dtype=float)
    natoms = len(labelsA)

    if len(labelsB) != natoms:                                         raise ValueError("MAP_VNMS: Geometries must contain the same number of atoms")
    if coordsA.shape != (natoms, 3) or coordsB.shape != (natoms, 3):   raise ValueError("MAP_VNMS: Coordinates must have shape (natoms, 3)")
    if any(vnm.mode.shape != (natoms, 3) for vnm in [*vnmsA, *vnmsB]): raise ValueError("MAP_VNMS: VNM dimensions do not match the geometries")
    if any(list(vnm.labels) != list(labelsA) for vnm in vnmsA):        raise ValueError("MAP_VNMS: VNM set A does not follow geometry A's atom order")
    if any(list(vnm.labels) != list(labelsB) for vnm in vnmsB):        raise ValueError("MAP_VNMS: VNM set B does not follow geometry B's atom order")

    # 1) Overlaps molecules. Needed to get the rotation
    overlap_debug = max(debug - 1, 0)
    isgood, labelsA_aligned, coordsA_aligned, labelsB_aligned, coordsB_aligned, atom_mapping, rotation = overlap_molecules(labelsA,coordsA,labelsB,coordsB,center_method=center_method,use_ext_info=use_ext_info,translate_to_ref=True,adjmat1=adjmatA,adjmat2=adjmatB,bond_orders1=bond_ordersA,bond_orders2=bond_ordersB,max_graph_mappings=max_graph_mappings,return_rotation=True,debug=overlap_debug)
    if not isgood: raise ValueError("MAP_VNMS: Molecular alignment failed")

    atom_mapping = np.asarray(atom_mapping, dtype=int)
    rotation     = np.asarray(rotation, dtype=float)

    if atom_mapping.shape != (natoms,):                        raise ValueError("MAP_VNMS: Invalid atom mapping")
    if rotation.shape != (3, 3):                               raise ValueError("MAP_VNMS: Invalid rotation matrix")
    if not np.array_equal(labelsA_aligned, labelsB_aligned):   raise ValueError("MAP_VNMS: Reordered atom labels do not match")

    # 2) Mass-weights modes and expresses set B in A's atom order and frame
    # A already defines the reference atom ordering and Cartesian frame.
    modesA = np.asarray([vnm.mass_weight_mode(permanent=False).reshape(-1) for vnm in vnmsA])
    # Apply the molecular atom mapping and rotation to every mode from B.
    modesB = []
    for vnm in vnmsB:
        mode = vnm.mass_weight_mode(permanent=False)
        mode = mode[atom_mapping]
        mode = (rotation @ mode.T).T
        modesB.append(mode.reshape(-1))
    modesB = np.asarray(modesB)
    if not np.all(np.isfinite(modesA)) or not np.all(np.isfinite(modesB)):  raise ValueError("MAP_VNMS: Mode vectors contain non-finite values")

    normsA = np.linalg.norm(modesA, axis=1)
    normsB = np.linalg.norm(modesB, axis=1)
    if np.any(normsA <= 0.0) or np.any(normsB <= 0.0):   raise ValueError("MAP_VNMS: Mode vectors must have nonzero norms")

    modesA_normalized = modesA / normsA[:, np.newaxis]
    modesB_normalized = modesB / normsB[:, np.newaxis]

    # 3) Calculates mode overlaps and frequency penalties
    # Absolute value removes the arbitrary sign of each eigenvector.
    overlap_matrix = np.abs(modesA_normalized @ modesB_normalized.T)

    freqsA = np.asarray([vnm.freq_cm for vnm in vnmsA], dtype=float)
    freqsB = np.asarray([vnm.freq_cm for vnm in vnmsB], dtype=float)

    frequency_difference = np.abs(freqsA[:, np.newaxis] - freqsB[np.newaxis, :])
    frequency_penalty    = f_weight * (frequency_difference / f_scale)**2

    # 4) Solves the one-to-one mode assignment
    cost_matrix      = (1.0 - overlap_matrix) + frequency_penalty
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    assignment = sorted(zip(row_ind, col_ind), key=lambda pair: pair[0])

    # 5) Stores mapping and molecular-alignment information
    mappings = []
    for i, j in assignment:
        mappings.append({
            # Positions in the supplied lists
            "modeA": int(i),
            "modeB": int(j),

            # Stored VNM indices
            "indexA": int(vnmsA[i].index),
            "indexB": int(vnmsB[j].index),

            # Assignment information
            "overlap": float(overlap_matrix[i, j]),
            "freqA": float(freqsA[i]),
            "freqB": float(freqsB[j]),
            "frequency_difference": float(frequency_difference[i, j]),
            "frequency_penalty": float(frequency_penalty[i, j]),
            "cost": float(cost_matrix[i, j]),
        })

    alignment_rmsd = np.sqrt(np.mean(np.sum((coordsA_aligned - coordsB_aligned)**2, axis=1)))
    alignment = {
        # Each reference position in A maps to this original position in B.
        "atom_mapping": atom_mapping.tolist(),
        "rotation": rotation,
        "labelsA": labelsA_aligned,
        "labelsB": labelsB_aligned,
        "coordsA": coordsA_aligned,
        "coordsB": coordsB_aligned,
        "rmsd": float(alignment_rmsd),
    }

    # 6) Builds ordered VNMs without modifying the input objects
    # B modes are independent copies expressed in A's atom ordering and frame.
    vnmsA_ordered = []
    vnmsB_ordered = []
    for i, j in assignment:
        vnmB_aligned  = deepcopy(vnmsB[j])
        modeB_aligned = modesB[j].reshape(natoms, 3)
        vnmB_aligned.set_mode(vnmsA[i].atomidxs, vnmsA[i].atnums, modeB_aligned[:, 0], modeB_aligned[:, 1], modeB_aligned[:, 2], is_mass_weighted=True)
        vnmsA_ordered.append(vnmsA[i])
        vnmsB_ordered.append(vnmB_aligned)

    # 7) Groups supplementary mapping and alignment results
    results = {
        "mappings": mappings,
        "alignment": alignment,
        "overlap_matrix": overlap_matrix,
        "frequency_difference_matrix": frequency_difference,
        "frequency_penalty_matrix": frequency_penalty,
        "cost_matrix": cost_matrix,
    }

    # 8) Prints diagnostics
    if debug > 0:
        assigned_overlaps         = np.asarray([mapping["overlap"] for mapping in mappings])
        assigned_freq_differences = np.asarray([mapping["frequency_difference"] for mapping in mappings])
        assigned_costs            = np.asarray([mapping["cost"] for mapping in mappings])

        print("MAP_VNMS: Mapping completed")
        print(f"  atoms                 = {natoms}")
        print(f"  modes                 = {len(mappings)}")
        print(f"  alignment RMSD        = {alignment_rmsd:.6f}")
        print(f"  overlap mean/median   = {np.mean(assigned_overlaps):.4f} / {np.median(assigned_overlaps):.4f}")
        print(f"  overlap min/max       = {np.min(assigned_overlaps):.4f} / {np.max(assigned_overlaps):.4f}")
        print(f"  overlap >= 0.8        = {np.mean(assigned_overlaps >= 0.8):.1%}")
        print(f"  frequency MAE/RMSE    = {np.mean(assigned_freq_differences):.2f} / {np.sqrt(np.mean(assigned_freq_differences**2)):.2f} cm-1")
        print(f"  frequency max error   = {np.max(assigned_freq_differences):.2f} cm-1")
        print(f"  assignment mean cost  = {np.mean(assigned_costs):.4f}")

    if debug > 1:
        weakest_mappings = sorted(mappings, key=lambda mapping: mapping["cost"], reverse=True)[:10]
        print("MAP_VNMS: Weakest assignments")
        print("  modeA  modeB    overlap    freqA    freqB    delta     cost")
        for mapping in weakest_mappings:
            print(f"  {mapping['indexA']:5d}  {mapping['indexB']:5d}    {mapping['overlap']:7.4f}  {mapping['freqA']:7.2f}  {mapping['freqB']:7.2f}  {mapping['frequency_difference']:7.2f}  {mapping['cost']:7.4f}")

    if debug > 2:
        print("MAP_VNMS: Atom mapping")
        print(atom_mapping.tolist())
        print("MAP_VNMS: Rotation matrix")
        print(np.round(rotation, 6))
        print("MAP_VNMS: Overlap matrix")
        print(np.round(overlap_matrix, 4))
        print("MAP_VNMS: Frequency-difference matrix (cm-1)")
        print(np.round(frequency_difference, 2))
        print("MAP_VNMS: Cost matrix")
        print(np.round(cost_matrix, 4))

    return vnmsA_ordered, vnmsB_ordered, results

######
def displace_coords_with_vnm(VNMs: list, initial_coord: list, which: list=[], which_side: str='positive', amplitude: int=6, debug: int=0):
    ### This function applies a displacement from the initial geometry 
    ### Using either 'all' VNMs provided, or only those whose index is in 'which'

    ## Frequencies must have eigenvectors (normal modes) stored
    if any(not vnm.has_mode for vnm in VNMs): 
        if debug > 0: print("One or more VNMs do not have mode. Stopping")
        return None

    ## Store Initial Coord
    new_coord = initial_coord.copy()    
    new_coord = np.array(new_coord)
    natoms    = len(initial_coord)
    if debug >= 1: print("initial coord:", new_coord[0])

    ## If empty which, then takes all
    if len(which) == 0: which = [vnm.index for vnm in VNMs]

    ## Applies Displacement
    for vnm in VNMs:
        assert vnm.has_mode and len(vnm.labels) == natoms
        mode = vnm.mass_weight_mode(permanent=False)
        ## The actual amount of displacement depends on the amplitude defined by the user
        ## ... and the freq_factor, which depends on the frequency
        if   vnm.freq < 20:   freq_factor = 0.1
        else:                 freq_factor = 0.2
        if vnm.index in which:
            if debug >= 1: print("displacing VNM with frequency:", vnm.freq_cm)
            if debug >= 1: print(f"using factor={amplitude*freq_factor}")
            
            for idx in range(natoms):
                ## Evaluate a temporary mass-weighted vector without modifying the VNM
                vector = mode[idx]

                ## Apply displacement to coordinates
                if   which_side.lower() == 'positive': displacement = vector*amplitude*freq_factor
                elif which_side.lower() == 'negative': displacement = -vector*amplitude*freq_factor
                new_coord[idx] = new_coord[idx]+displacement
                if idx == 0 and debug >= 1: print("displaced coord:", new_coord[idx])
        new_coord.reshape(natoms,3)
    return new_coord

####
def geom_sampling_from_vnm(labels, coord, freqs, qini: list=None, T: float=0.0, n_samples: int=10, sigma_damp_factor: float=1, freq_bottom_limit: float=80, check_adjacencies: bool=True, debug: int=0):
    from scope.connectivity             import get_adjmatrix
    from scope.operations.vecs_and_mats import normalize
    """
    Generates a set of geometries by sampling along vibrational normal modes (VNM) of a molecule.
    This function perturbs the input geometry along its vibrational normal modes, producing a set of 
    geometries that represent possible thermal fluctuations. Optionally, it checks that 
    the molecular connectivity (adjacency matrix) is preserved in the sampled geometries.
    Parameters
    ----------
    labels : list of str
        Atomic labels (element symbols) for each atom in the molecule.
    coord : array-like, shape (N_atoms, 3)
        Cartesian coordinates of the atoms (in Bohr).
    freqs : list
        List of frequency objects, each with attributes:
            - freq_cm: frequency in cm^-1
            - freq: frequency in atomic units
            - mode: eigenvector with shape `(N_atoms, 3)`
            - has_mode: boolean indicating presence of the eigenvector
            - is_mass_weighted: boolean describing the stored eigenvector
    qini  : list
        initial Q coordinates associated with the cartesian coordinates provided as coord
    T : float, optional
        Temperature in Kelvin for thermal sampling. If 0, only zero-point motion is considered.
    n_samples : int, optional
        Number of geometries to generate.
    sigma_damp_factor: float, optional
        Reduces large unphysical displacements. Specially important for low-frequency modes
    freq_bottom_limit : float, optional
        Minimum frequency (in cm^-1) to include in sampling. Modes below this are ignored.
    check_adjacencies : bool, optional
        If True, only accept geometries that preserve the original adjacency matrix.
    debug : int, optional
        Debug verbosity level. Higher values print more information.
    Returns
    -------
    geometries : np.ndarray, shape (n_samples, N_atoms, 3)
        Array of sampled geometries (Cartesian coordinates in Bohr).
    q_coords : np.ndarray, shape (n_samples, N_modes)
        Array of sampled normal mode displacements for each geometry.
    energies : np.ndarray, shape (n_samples,)
        Array of harmonic energies (in atomic units) for each sampled geometry.
    Notes
    -----
    - Modes are converted to mass-weighted coordinates locally without modifying the VNM objects.
    - The sampling is performed using a normal distribution for each mode, with width determined by thermal fluctuations.
    - If `check_adjacencies` is True, geometries that change the molecular connectivity are discarded.
    - Large unphysical displacements are discarded, which breaks the expected energy distribution of the resulting geometries. (np.mean(energies) should approach ZPE/2) 
    - The function may return fewer than `n_samples` geometries if the adjacency check fails frequently.
    """
    def taper_weight(freq_cm, f0=50.0, width=10.0):
        x = (freq_cm - f0) / width
        return 0.5 * (1 + np.tanh(x))

    ## Frequencies must have eigenvectors stored
    if any(not freq.has_mode for freq in freqs): 
        raise ValueError("One or more VNMs do not have eigenvectors. Stopping")

    # Stores the original adjacency matrix
    if check_adjacencies: 
        isgood, original_adjmat, original_adjnum = get_adjmatrix(labels, coord, smart=True, debug=debug)
        if not isgood: print("Warning: Initial adjacency matrix might have an issue")

    # Extract and manage data from input
    coord   = np.array(coord) 
    modes   = np.asarray([freq.mass_weight_mode(permanent=False).reshape(-1) for freq in freqs])
    N_atoms = len(coord)    
    N_modes = len(freqs)
    if qini is None: qini = np.zeros((N_modes))
    else:            qini = np.asarray(qini)
    if debug > 1: print("Number of modes:", len(modes))
    if debug > 1: print("First mode:", modes[0])
    if debug > 1: print("First mode norm:", normalize(modes[0]))

    # Initializes and Runs Main While loop
    geometries = []
    q_coords   = []
    energies   = []
    maxcount   = n_samples * 100
    count      = 0
    while len(geometries) < n_samples or count >= maxcount:
        displacement = np.zeros(3 * N_atoms)

        q_tot = qini.copy()
        q_vec = []
        e_harm = 0.0
        for i in range(N_modes):
            if freqs[i].freq_cm <= 0: continue  # skip imaginary modes

            omega = freqs[i].freq                               
            sigma_q = np.sqrt(constants.hbar / (2 * omega))     

            if T > 0:
                coth = 1 / np.tanh(constants.hbar * omega / (2 * constants.boltz_au * T))
                sigma_q *= np.sqrt(coth)
            else: coth = 1.0

            # Apply smooth tapering based on frequency in cm^-1
            if freq_bottom_limit > 0: weight = taper_weight(freqs[i].freq_cm, f0=freq_bottom_limit)
            else:                     weight = 1.0

            # Applies weight and limitations
            sigma_q *= np.sqrt(weight)
            sigma_q  = sigma_q/(freqs[i].freq_cm**sigma_damp_factor) ## a factor tunes down the displacement in the ith-VNM depending on its frequency
            sigma_q  = np.min((np.abs(sigma_q), 10))                 ## displacement is also limited to Q=10 

            # Sample a displacement
            q_i = np.random.normal(loc=0.0, scale=sigma_q)          ## Coordinates of Q for this VNM
            q_vec.append(q_i)                                       ## Collection of Q for all VNM to be applied to this sample
            q_tot[i] += q_i                                         ## Collection of Accumulated Q (including the initial ones) 
            e_harm += 0.5 * q_tot[i]**2 * omega**2
            displacement += q_i * modes[i]                 

            if debug > 0 and i < 10: print(f"  Mode {i}: freq_cm = {freqs[i].freq_cm}, weight = {weight:.3f}, q_i= {q_i:.3f}, sigma_q = {sigma_q:.3f}")
            if debug > 1: print(f"  max eigenvector component (mass-weighted)    = {np.max(np.abs(modes[i])):.3e}")
            if debug > 1: print(f"  max displacement contribution from this mode = {np.max(np.abs(q_i * modes[i])):.3e}")
            if debug > 1: print(f"  {i=} displacement: {np.round(displacement[i],3)}")

        if debug > 1: print(f"Coord:      {np.round(coord[0],3)}")
        if debug > 1: print(f"modes[i]    {np.round(modes[0][0:3],5)}")

        # Convert to Cartesian displacement
        if debug > 1: print(f"Cart_Disp:  {np.round(displacement[0:3],3)}")
        displaced_coord   = coord.flatten() + displacement
        if debug > 1: print(f"Disp_Coord1 [bohr]: {np.round(displaced_coord[0:3],3)}")
        displaced_coord   = displaced_coord.reshape((N_atoms, 3))

        if check_adjacencies: 
            # Apply the same smart construction because its internal covalent factor is intentionally not exposed
            isgood, new_adjmat, new_adjnum = get_adjmatrix(labels, displaced_coord, smart=True)
            if isgood and np.array_equal(original_adjmat, new_adjmat):
                if debug > 0: print(f"Geometry {count} ACCEPTED")
                geometries.append(displaced_coord)
                q_coords.append(q_tot)
                energies.append(e_harm)
            else:
                if debug > 0: print(f"Geometry {count} DISCARDED due to adjacency mismatch")
                if debug > 1: print(f"{new_adjnum=}")
                if debug > 1: print(f"{original_adjnum=}")
                if debug > 1: print(f"{new_adjnum-original_adjnum}")
        else:
            geometries.append(displaced_coord)
            q_coords.append(q_tot)
            energies.append(e_harm)
        count += 1
        if count >= maxcount: 
            print(f"Warning: Reached maximum count of {maxcount} with only {len(geometries)}/{n_samples} samples created.")
            break
    if debug > 0: print(f"------------------------------------------------------------------")
    if debug > 0: print(f"Sampling Produced {len(geometries)} structures in {count} attempts")
    if debug > 0: print(f"------------------------------------------------------------------")
    return np.array(geometries), np.array(q_coords), np.array(energies)

####
def apply_q_displacement(x_ref, q_coords, freqs):
    """
    Applies normal mode displacements to a reference geometry using given mode amplitudes.

    Parameters
    ----------
    x_ref : array-like, shape (natoms, 3)
        Reference geometry coordinates of the molecule (Cartesian coordinates).
    q_coords : array-like, shape (nmodes,)
        Amplitudes for each normal mode (dimensionless normal coordinates).
    freqs : list
        List of objects representing vibrational modes. Each object must have an attribute
        `eigenvec_format2` containing the mode eigenvector (flattened, length 3*natoms).

    Returns
    -------
    x_new : ndarray, shape (natoms, 3)
        New geometry coordinates after applying the normal mode displacements.

    Raises
    ------
    AssertionError
        If the length of `q_coords` does not match the number of modes in `freqs`.
    """

    ## Prepares initial data and reviews length
    natoms  = len(x_ref)                # Number of atoms
    x_ref   = np.array(x_ref)           # Reference geometry
    modes   = np.array([freq.eigenvec_format2 for freq in freqs])  
    assert(len(q_coords) == len(modes)), "Length of q_coords must match number of modes"

    ## Main loop, accumulates displacement for each mode
    displacement = np.zeros(3 * natoms)
    for i in range(len(modes)):
        displacement += q_coords[i] * modes[i]        
    x_new = x_ref.flatten() + displacement

    return x_new.reshape((natoms, 3))

####
def project_to_normal_modes(l1, x1, l2, x2, freqs, debug: int=0):
    from scope.overlap import rmsd
    """
    Aligns x_geom to x_ref and projects the displacement onto normal modes.

    Parameters:
    - l1, x1: (N_atoms, 3) labels, and reference geometry (usually minimum energy structure)
    - l2, x2: (N_atoms, 3) labels, and target geometry
    - freqs: (N_modes, 3*N_atoms) normal mode eigenvectors of the REFERENCE geometry in mass-weighted coordinates (Gaussian16-like)

    Returns:
    - q_coords: (N_modes,) generalized coordinates q_i for each normal mode
    """

    x1      = np.array(x1)    # Reference geometry
    l1      = np.array(l1)    # Reference geometry
    x2      = np.array(x2)    # Target geometry
    l2      = np.array(l2)    # Target geometry
    modes   = np.array([freq.eigenvec_format2 for freq in freqs])  # Extract Eigenvectors from frequencies. [units? Assuming Bohr/sqrt(amu)]

    # Step 2: Displacement 
    dx_cart = x2.flatten() - x1.flatten()

    # Step 3: Projection 
    q_coords = [i for i in np.dot(modes, dx_cart)]

    # Step 4: checks error:
    x3 = apply_q_displacement(x1, q_coords, freqs)
    for i in range(len(x1)):
        print(i, l1[i], rmsd(l1, x1, l1, x3, atom_idxs=[i], reorder=False))
    print("Error of the projection (RMSD):", rmsd(l1, x1, l1, x3, reorder=False))

    return q_coords

#############
## For FPS ##
#############
def beta_distance(Q1, Q2, freq_cm, T: float=300, debug: int=0):
    ### This is a similarity function to compare two sets of Q values. Meant to be used with the furthest point sampling function in Other
    omegas = np.array(freq_cm)*constants.cm2har
    dist = 0
    beta   = 1/(constants.boltz_au*T)
    ps     = np.exp(-beta*omegas)
    ps     = ps/np.sum(ps)
    for a, b, f, o in zip(Q1, Q2, ps, omegas):
        dist += np.abs(a-b) * f
    return dist

def custom_q_distance(Q1, Q2, freq_cm):
    ### This is a similarity function to compare two sets of Q values. Meant to be used with the furthest point sampling function in Other
    dist = 0
    for a, b, f in zip(Q1, Q2, freq_cm):
        dist += np.abs(a-b) * (1/f)
    return dist

def euclidean_q_distance(Q1, Q2):
    ### This is a similarity function to compare two sets of Q values. Meant to be used with the furthest point sampling function in Other
    dist = np.linalg.norm(Q1-Q2)
    return dist

def pairwise_distance_matrix(data, dist_func):
    ### This function computes the similarity matrix using one of the metrics above. Meant to be used with the furthest point sampling function in Other
    """Computes the pairwise distance matrix using a custom distance function."""
    n = data.shape[0]
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dist = dist_func(data[i], data[j]) #, freqs_cm)
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    return dist_matrix
