import numpy as np
from Scope import Constants

####
def map_vnms(vnmsA, vnmsB, debug: int=0):
    """
    Map vibrational normal modes (VNMs) between two sets, without
    dividing into degenerate blocks.

    Parameters
    ----------
    vnmsA, vnmsB : arrays (n_modes)
        Lists of VNM-class objects

    Returns
    -------
    mappings : list of dict
        Each dict has:
        - modeA: index in A
        - modeB: matched index in B
        - overlap: absolute dot product after weighting/normalization
        - freqA, freqB: the corresponding frequencies
    """

    # Extracts data from VNMs
    modesA = [v.mode_format2 for v in vnmsA]
    modesB = [v.mode_format2 for v in vnmsB] 
    freqsA = [v.freq_cm for v in vnmsA]
    freqsB = [v.freq_cm for v in vnmsB] 

    # Normalize (mass-weight should not be needed)
    from Scope.Operations.Vecs_and_Mats import normalize
    modesA_proc = normalize(modesA)
    modesB_proc = normalize(modesB)

    nA, nB = modesA_proc.shape[0], modesB_proc.shape[0]
    if nA != nB: raise ValueError("Number of modes in A and B must be equal.")

    # Overlap matrix (absolute values)
    S = np.abs(modesA_proc @ modesB_proc.T)

    # Hungarian algorithm (maximize total overlap → minimize -S)
    from scipy.optimize import linear_sum_assignment
    row_ind, col_ind = linear_sum_assignment(-S)

    # Build mapping list
    mappings = []
    for i, j in zip(row_ind, col_ind):
        mappings.append({"modeA": int(i),"modeB": int(j),"overlap": float(S[i, j]),"freqA": float(freqsA[i]),"freqB": float(freqsB[j])})
    return mappings

######
def displace_coords_with_vnm(VNMs: list, initial_coord: list, which: list=[], which_side: str='positive', amplitude: int=6, debug: int=0):
    ### This function applies a displacement from the initial geometry 
    ### Using either 'all' VNMs provided, or only those whose index is in 'which'

    ## Frequencies must have eigenvectors stored
    if any(not vnm.has_mode for vnm in VNMs): 
        if debug > 0: print("One or more VNMs do not have eigenvectors. Stopping")
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
        assert vnm.haseigenvec and len(vnm.xs) == natoms
        ## The actual amount of displacement depends on the amplitude defined by the user
        ## ... and the freq_factor, which depends on the frequency
        if   vnm.freq < 20:   freq_factor = 0.1
        else:                 freq_factor = 0.2
        if vnm.index in which:
            if debug >= 1: print("displacing VNM with frequency:", vnm.freq_cm)
            if debug >= 1: print(f"using factor={amplitude*freq_factor}")
            
            for idx in range(natoms):
                ## Evaluate the vector, here it expects to receive a mass-weighted eigenvector, as in Gaussian
                vector = vnm.eigenvec_format1[idx]

                ## Apply displacement to coordinates
                if   which_side.lower() == 'positive': displacement = vector*amplitude*freq_factor
                elif which_side.lower() == 'negative': displacement = -vector*amplitude*freq_factor
                new_coord[idx] = new_coord[idx]+displacement
                if idx == 0 and debug >= 1: print("displaced coord:", new_coord[idx])
        new_coord.reshape(natoms,3)
    return new_coord

####
def geom_sampling_from_vnm(labels, coord, freqs, qini: list=None, T: float=0.0, n_samples: int=10, freq_bottom_limit: float=0, check_adjacencies: bool=True, debug: int=0):
    from Scope.Adapted_from_cell2mol    import get_adjmatrix
    from Scope.Operations.Vecs_and_Mats import normalize
    """
    Generate a set of geometries by sampling along vibrational normal modes (VNM) of a molecule.
    This function perturbs the input geometry along its vibrational normal modes, producing a set of 
    geometries that represent possible thermal or quantum fluctuations. Optionally, it checks that 
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
            - eigenvec_format2: mass-weighted eigenvector (array)
            - haseigenvec: boolean indicating presence of eigenvector
    qini  : list
        initial Q coordinates associated with the cartesian coordinates provided as coord
    T : float, optional
        Temperature in Kelvin for thermal sampling. If 0, only zero-point motion is considered.
    n_samples : int, optional
        Number of geometries to generate.
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
    - The function assumes that the input frequencies and eigenvectors are mass-weighted and in atomic units.
    - Sampling is performed using a normal distribution for each mode, with width determined by quantum and/or thermal fluctuations.
    - If `check_adjacencies` is True, geometries that change the molecular connectivity are discarded.
    - The function may return fewer than `n_samples` geometries if the adjacency check fails frequently.
    """
    def taper_weight(freq_cm, f0=50.0, width=10.0):
        x = (freq_cm - f0) / width
        return 0.5 * (1 + np.tanh(x))

    ## Frequencies must have eigenvectors stored
    if any(not freq.has_mode for freq in freqs): 
        if debug > 0: print("One or more VNMs do not have eigenvectors. Stopping")
        return None

    # Stores the original adjacency matrix
    if check_adjacencies: 
        isgood, original_adjmat, original_adjnum = get_adjmatrix(labels, coord, cov_factor=1.3, metal_factor=1.0)
        if not isgood: print("Warning: Initial adjacency matrix might have an issue")

    # Extract and manage data from input
    coord   = np.array(coord) 
    modes   = np.array([freq.mode_format2 for freq in freqs])   # Extract Eigenvectors from frequencies. [units? Assuming Bohr/sqrt(amu)]
    N_atoms = len(coord)    
    N_modes = len(freqs)
    if qini is None: qini = np.zeros((N_modes))
    else:            qini = np.asarray(qini)
    if debug > 0: print("First mode norm:", normalize(modes[0].mode_format2))

    # Initializes and Runs Main While loop
    geometries = []
    q_coords   = []
    energies   = []
    maxcount   = n_samples * 100
    count      = 0
    while len(geometries) < n_samples or count >= maxcount: # for _ in range(n_samples):
        displacement = np.zeros(3 * N_atoms)

        q_tot = qini.copy()
        q_vec = []
        e_harm = 0.0
        for i in range(N_modes):
            if freqs[i].freq_cm <= freq_bottom_limit:
                continue  # skip high frequency or imaginary modes

            omega = freqs[i].freq                               
            sigma_q = np.sqrt(Constants.hbar / (2 * omega))     

            if T > 0:
                coth = 1 / np.tanh(Constants.hbar * omega / (2 * Constants.boltz_au * T))
                sigma_q *= np.sqrt(coth)
            else: coth = 1.0

            # Apply smooth tapering based on frequency in cm^-1
            weight = taper_weight(freqs[i].freq_cm, f0=freq_bottom_limit)
            sigma_q *= np.sqrt(weight)

            # Sample a displacement
            sigma_q = sigma_q/(freqs[i].freq_cm**1.8)
            sigma_q = np.min((np.abs(sigma_q), 10))
            q_i = np.random.normal(loc=0.0, scale=sigma_q)      ## Coordinates of Q for this VNM
            q_vec.append(q_i)                                   ## Collection of Q for all VNM to be applied to this sample
            q_tot[i] += q_i                                     ## Collection of Accumulated Q (including the initial ones) 
            e_harm += 0.5 * q_tot[i]**2 * omega**2
            displacement += q_i * modes[i]                 

            if debug > 0 and i < 10: print(f"  Mode {i}: freq_cm = {freqs[i].freq_cm}, q_ini[i]= {q_ini[i]:.3f}, q_vec[i]= {q_vec[i]:.3f}, sigma_q = {sigma_q:.3f}")
            if debug > 1: print(f"  max eigenvector component (mass-weighted)    = {np.max(np.abs(modes[i])):.3e}")
            if debug > 1: print(f"  max displacement contribution from this mode = {np.max(np.abs(q_i * modes[i])):.3e}")
            if debug > 1: print(f"  {i=} displacement[0]: {np.round(displacement[0],3)}")

        if debug > 1: print(f"Coord:      {np.round(coord[0],3)}")
        if debug > 1: print(f"modes[i]    {np.round(modes[0][0:3],5)}")

        # Convert to Cartesian displacement
        if debug > 0: print(f"Cart_Disp:  {np.round(displacement[0:3],3)}")
        displaced_coord   = coord.flatten() + displacement
        if debug > 0: print(f"Disp_Coord1 [bohr]: {np.round(displaced_coord[0:3],3)}")
        displaced_coord   = displaced_coord.reshape((N_atoms, 3))

        if check_adjacencies: 
            isgood, new_adjmat, new_adjnum = get_adjmatrix(labels, displaced_coord, cov_factor=1.3, metal_factor=1.0)
            if isgood and np.array_equal(original_adjmat, new_adjmat):
                geometries.append(displaced_coord)
                q_coords.append(q_tot)
                energies.append(e_harm)
            else:
                if debug > 0: print("Discarded geometry due to adjacency mismatch")
                if debug > 0: print(f"{new_adjnum=}")
                if debug > 0: print(f"{original_adjnum=}")
                if debug > 0: print(f"{new_adjnum-original_adjnum}")
        else:
            geometries.append(displaced_coord)
            q_coords.append(q_tot)
            energies.append(e_harm)
        count += 1
        if count >= maxcount: 
            print(f"Warning: Reached maximum count of {maxcount} with only {len(geometries)}/{n_samples} samples created.")
            break
    if debug > 0: print(f"Produced {len(geometries)} structures in {count} attempts")
    return np.array(geometries), np.array(q_coords), np.array(energies)

####
def expected_vibrational_energy(omega_au, T):
    ### Simple way to compute the vibrational energy of a given VNM at temperature T. 
    omega_au = np.array(omega_au)
    if T == 0:
        # Pure zero-point energy
        return 0.5 * Constants.hbar * np.sum(omega_au)
    else:
        # Finite temperature using coth expression
        x = (Constants.hbar * omega_au) / (2 * Constants.boltz_au * T)
        coth_x = (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))  # coth(x) = (e^x + e^(-x)) / (e^x - e^(-x))
        E_vib = 0.5 * Constants.hbar * np.sum(omega_au * coth_x)
        return E_vib
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
    from Scope.Other import rmsd
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
## For FPS ##
#############
def beta_distance(Q1, Q2, freq_cm, T: float=300, debug: int=0):
    ### This is a similarity function to compare two sets of Q values. Meant to be used with the furthest point sampling function in Other
    omegas = np.array(freq_cm)*Constants.cm2har
    dist = 0
    beta   = 1/(Constants.boltz_au*T)
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