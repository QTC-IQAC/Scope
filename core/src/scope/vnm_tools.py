import numpy as np
from scope import constants

####
def displace_neg_freqs(ini_coord, VNMs: object, amplitude: float=10, debug: int=0) -> np.ndarray:
    """Displace each imaginary mode by amplitude in sqrt(electron mass)*Bohr (default 10); return Angstrom coordinates."""
    neg_VNMs = [vnm for vnm in VNMs if vnm.freq_cm < 0.0]
    if debug > 0: print(f"DISPLACE_NEG_FREQS: Displacing {len(neg_VNMs)} imaginary modes with Q={amplitude:.5g}")
    return displace_coords_with_vnm(neg_VNMs, ini_coord, amplitude=amplitude, debug=max(debug - 1, 0))

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
def displace_coords_with_vnm(VNMs: list, initial_coord: list, which: list=None, which_side: str='positive', amplitude: float=10, debug: int=0):
    """Displace a geometry by a fixed atomic-unit Q along each selected mode.

    Parameters:
        VNMs (list):                   Modes in the geometry's atom order and frame, including imaginary modes.
        initial_coord (array):         Cartesian coordinates in Angstrom with shape (natoms, 3).
        which (list | None):           VNM indices to displace; None or an empty list selects all.
        which_side (str):              'positive' or 'negative', relative to the stored eigenvector phase.
        amplitude (float):             Q per selected mode in sqrt(electron mass)*Bohr, without frequency scaling.
        debug (int):                   Verbosity level.

    Returns:
        np.ndarray: New floating-point Cartesian coordinates in Angstrom.

    Notes:
        Both stored mode representations are normalized locally in mass-weighted coordinates,
        divided by sqrt(mass in electron masses), and converted from Bohr to Angstrom.
        Input coordinates and VNMs are not modified. Modes must share atom order and masses.
        No thermal distribution or harmonic energy is assigned to this fixed displacement.
        The former amplitude*frequency_factor*sqrt(mass)*mode rule had no consistent Q convention.
    """
    # 0) Checks the geometry and mode selection
    coord = np.array(initial_coord, dtype=float, copy=True)
    VNMs  = list(VNMs)
    if coord.ndim != 2 or coord.shape[1] != 3 or len(coord) == 0 or not np.all(np.isfinite(coord)): raise ValueError("DISPLACE_COORDS_WITH_VNM: Expected finite coordinates with shape (natoms, 3)")
    if not isinstance(which_side, str) or which_side.lower() not in ("positive", "negative"): raise ValueError("DISPLACE_COORDS_WITH_VNM: which_side must be 'positive' or 'negative'")
    if not np.isfinite(amplitude): raise ValueError("DISPLACE_COORDS_WITH_VNM: amplitude must be finite")
    available_indices = [vnm.index for vnm in VNMs]
    if which is not None and any(index not in available_indices for index in which): raise ValueError("DISPLACE_COORDS_WITH_VNM: A requested VNM index was not found")
    selected = VNMs if which is None or len(which) == 0 else [vnm for vnm in VNMs if vnm.index in which]
    sign     = 1.0 if which_side.lower() == "positive" else -1.0

    # 1) Converts each selected direction to Angstrom per atomic-unit Q
    displacement      = np.zeros_like(coord)
    ref_labels  = None
    ref_indices = None
    ref_masses  = None
    for vnm in selected:
        if not getattr(vnm, "has_mode", False): 
            raise ValueError("DISPLACE_COORDS_WITH_VNM: A selected VNM has no eigenvector")
        if len(vnm.labels) != len(coord):       
            raise ValueError("DISPLACE_COORDS_WITH_VNM: VNM atom count does not match the geometry")
        if ref_labels is not None and (list(vnm.labels) != ref_labels or list(vnm.atomidxs) != ref_indices): 
            raise ValueError("DISPLACE_COORDS_WITH_VNM: VNMs must share the same atom order")
        masses = np.asarray(vnm.masses, dtype=float)
        if masses.shape != (len(coord),) or not np.all(np.isfinite(masses)) or np.any(masses <= 0): 
            raise ValueError("DISPLACE_COORDS_WITH_VNM: Atomic masses must be finite and positive")
        if ref_masses is not None and not np.allclose(masses, ref_masses, rtol=1e-12, atol=0): 
            raise ValueError("DISPLACE_COORDS_WITH_VNM: VNMs must use the same atomic masses")

        ref_labels  = list(vnm.labels)
        ref_indices = list(vnm.atomidxs)
        ref_masses  = masses
        mode_mw = vnm.mass_weight_mode(permanent=False)
        if mode_mw.shape != coord.shape or not np.all(np.isfinite(mode_mw)): 
            raise ValueError("DISPLACE_COORDS_WITH_VNM: Invalid mode shape or components")
        mode_norm = np.linalg.norm(mode_mw)
        if not np.isfinite(mode_norm) or mode_norm <= 0: 
            raise ValueError("DISPLACE_COORDS_WITH_VNM: Mode norm must be finite and positive")
        cart_mode = constants.bohr2angs * (mode_mw / mode_norm) / np.sqrt(masses[:, None] * constants.amu2au)
        displacement += sign * amplitude * cart_mode
        if debug > 1: 
            print(f"DISPLACE_COORDS_WITH_VNM: Mode {vnm.index}: Q={sign * amplitude:.5g}")
    if debug > 0: 
        print(f"DISPLACE_COORDS_WITH_VNM: Displaced {len(selected)} modes with Q={sign * amplitude:.5g} each")
    
    return coord + displacement

####
def geom_sampling_from_vnm(labels, coord, freqs, qini: list=None, T: float=0.0, n_samples: int=10, sigma_damp_factor: float=0, freq_bottom_limit: float=0, freq_top_limit: float=None, check_adjacencies: bool=True, debug: int=0):
    """Sample Cartesian geometries from harmonic normal coordinates.

    Parameters:
        labels, coord:                 Atom symbols and coordinates in Angstrom, in mode order and frame.
        freqs (list):                  VNMs; only positive frequencies are sampled.
        qini (list):                   Parent normal coordinates in sqrt(electron mass)*Bohr, one per VNM.
                                       coord must already include these displacements from the reference geometry.
        T (float):                     Temperature in K; zero includes zero-point position fluctuations.
        n_samples (int):               Number of accepted geometries requested.
        sigma_damp_factor (float):      Optional exponent dividing sigma by freq_cm**exponent; zero disables it.
        freq_bottom_limit (float):      Optional smooth taper midpoint in cm^-1; nonpositive disables it.
        check_adjacencies (bool):       Reject changed connectivity; at most 100*n_samples attempts.
        freq_top_limit (float | None):  Smooth upper taper midpoint in cm^-1; None disables it.
        debug (int):                   Verbosity; child functions receive debug - 1.

    Returns:
        tuple: Geometries (n, natoms, 3) in Angstrom, cumulative Q (n, nmodes), and
               harmonic potential increments (n,) in Hartree, including the parent Q.

    Notes:
        Mass-weighted mode directions are normalized locally, divided by sqrt(mass in
        electron masses), and converted from Bohr to Angstrom. VNMs are not modified.
        Exact quantitative conversion requires the masses used by the frequency calculation.
        No Q or sigma cap is applied. The quantum variance is coth(omega/(2*kB*T))/(2*omega).
        Without attenuation or rejection, the mean potential at T=0 is half the zero-point energy.
        The taper multiplies variance by 0.5*(1+tanh((freq_cm-midpoint)/10)); it is not a hard cutoff.
        Frequency damping, tapering, connectivity rejection and subsequent FPS change the distribution.
        Skipped nonpositive modes must have zero qini. Seeding numpy.random reproduces draws.
        The upper taper multiplies variance by 0.5*(1-tanh((freq_cm-midpoint)/10)).
        Both taper weights multiply; parent Q and its cumulative energy remain included.
        Labels detect element-order mismatches; the caller must preserve repeated-atom order and frame.
    """
    from scope.connectivity import get_adjmatrix

    # 0) Checks and arranges arrays
    labels  = list(labels)
    coord   = np.asarray(coord, dtype=float)
    freqs   = list(freqs)
    N_atoms = len(labels)
    N_modes = len(freqs)
    if coord.shape != (N_atoms, 3) or N_atoms == 0 or not np.all(np.isfinite(coord)): raise ValueError("GEOM_SAMPLING_FROM_VNM: Expected finite coordinates with shape (natoms, 3)")
    if not isinstance(n_samples, (int, np.integer)) or isinstance(n_samples, (bool, np.bool_)) or n_samples < 0: raise ValueError("GEOM_SAMPLING_FROM_VNM: n_samples must be a nonnegative integer")
    if not np.isfinite(T) or T < 0: raise ValueError("GEOM_SAMPLING_FROM_VNM: T must be finite and nonnegative")
    if not np.isfinite(sigma_damp_factor) or not np.isfinite(freq_bottom_limit): raise ValueError("GEOM_SAMPLING_FROM_VNM: Attenuation settings must be finite")
    if freq_top_limit is not None and (not np.isfinite(freq_top_limit) or freq_top_limit < 0): raise ValueError("GEOM_SAMPLING_FROM_VNM: freq_top_limit must be None or finite and nonnegative")
    if N_modes == 0: raise ValueError("GEOM_SAMPLING_FROM_VNM: No VNMs were provided")

    qini    = np.zeros(N_modes) if qini is None else np.array(qini, dtype=float, copy=True)
    freq_cm = np.asarray([freq.freq_cm for freq in freqs], dtype=float)
    omegas  = np.asarray([freq.freq for freq in freqs], dtype=float)
    if qini.shape != (N_modes,) or not np.all(np.isfinite(qini)): raise ValueError("GEOM_SAMPLING_FROM_VNM: qini must contain one finite amplitude per VNM")
    if not np.all(np.isfinite(freq_cm)) or not np.all(np.isfinite(omegas)): raise ValueError("GEOM_SAMPLING_FROM_VNM: Frequencies must be finite")
    positive = np.flatnonzero(freq_cm > 0)
    if np.any(omegas[positive] <= 0): raise ValueError("GEOM_SAMPLING_FROM_VNM: Positive modes require positive atomic-unit frequencies")
    if np.any(qini[freq_cm <= 0] != 0): raise ValueError("GEOM_SAMPLING_FROM_VNM: Skipped nonpositive modes require zero qini")
    if n_samples == 0: return np.empty((0, N_atoms, 3)), np.empty((0, N_modes)), np.empty(0)

    # 1) Converts normalized mass-weighted eigenvectors to Angstrom per atomic-unit Q
    modes             = np.zeros((N_modes, 3 * N_atoms))
    sigmas            = np.zeros(N_modes)
    reference_masses  = None
    reference_indices = None
    for i in positive:
        freq = freqs[i]
        if not getattr(freq, "has_mode", False): raise ValueError("GEOM_SAMPLING_FROM_VNM: A sampled VNM has no eigenvector")
        if list(freq.labels) != labels: raise ValueError("GEOM_SAMPLING_FROM_VNM: VNM atom labels must match the geometry order")
        masses = np.asarray(freq.masses, dtype=float)
        if masses.shape != (N_atoms,) or not np.all(np.isfinite(masses)) or np.any(masses <= 0): raise ValueError("GEOM_SAMPLING_FROM_VNM: Atomic masses must be finite and positive")
        if reference_masses is not None and not np.allclose(masses, reference_masses, rtol=1e-12, atol=0): raise ValueError("GEOM_SAMPLING_FROM_VNM: VNMs must use the same atomic masses")
        if reference_indices is not None and list(freq.atomidxs) != reference_indices: raise ValueError("GEOM_SAMPLING_FROM_VNM: VNMs must use the same atom indices")
        reference_masses  = masses
        reference_indices = list(freq.atomidxs)
        mode_mw = freq.mass_weight_mode(permanent=False)
        if mode_mw.shape != (N_atoms, 3) or not np.all(np.isfinite(mode_mw)): raise ValueError("GEOM_SAMPLING_FROM_VNM: Invalid mode shape or components")
        mode_norm = np.linalg.norm(mode_mw)
        if not np.isfinite(mode_norm) or mode_norm <= 0: raise ValueError("GEOM_SAMPLING_FROM_VNM: Mode norm must be finite and positive")
        cart_mode = constants.bohr2angs * (mode_mw / mode_norm) / np.sqrt(masses[:, None] * constants.amu2au)
        modes[i]  = cart_mode.reshape(-1)

        # 2) Sets the thermal width and any explicitly requested empirical attenuation
        variance = 1.0 / (2.0 * omegas[i])
        if T > 0: variance /= np.tanh(omegas[i] / (2.0 * constants.boltz_au * T))
        weight = 0.5 * (1 + np.tanh((freq_cm[i] - freq_bottom_limit) / 10.0)) if freq_bottom_limit > 0 else 1.0
        if freq_top_limit is not None: weight *= 0.5 * (1 - np.tanh((freq_cm[i] - freq_top_limit) / 10.0))
        sigmas[i] = np.sqrt(variance * weight) / freq_cm[i]**sigma_damp_factor
        if not np.isfinite(sigmas[i]): raise ValueError("GEOM_SAMPLING_FROM_VNM: Sampling width is not finite")
        if debug > 1: print(f"GEOM_SAMPLING_FROM_VNM: Mode {i + 1}: frequency={freq_cm[i]:.3f} cm^-1, sigma={sigmas[i]:.5g}, taper={weight:.3f}")

    # 3) Samples increments around the supplied parent geometry and cumulative Q
    if check_adjacencies:
        isgood, original_adjmat, _ = get_adjmatrix(labels, coord, smart=True, debug=max(debug - 1, 0))
        if not isgood: print("GEOM_SAMPLING_FROM_VNM: Warning: Initial adjacency matrix might have an issue")
    geometries, q_coords, energies = [], [], []
    maxcount = n_samples * 100
    count    = 0
    while len(geometries) < n_samples and count < maxcount:
        q_tot = qini.copy()
        for i in positive: q_tot[i] += np.random.normal(loc=0.0, scale=sigmas[i])
        displacement    = (q_tot - qini) @ modes
        displaced_coord = coord + displacement.reshape((N_atoms, 3))
        e_harm          = 0.5 * np.sum((omegas[positive] * q_tot[positive])**2)
        count += 1

        if check_adjacencies:
            isgood, new_adjmat, _ = get_adjmatrix(labels, displaced_coord, smart=True, debug=max(debug - 1, 0))
            if not isgood or not np.array_equal(original_adjmat, new_adjmat): continue
        geometries.append(displaced_coord)
        q_coords.append(q_tot)
        energies.append(e_harm)

    if len(geometries) < n_samples: print(f"GEOM_SAMPLING_FROM_VNM: Reached {maxcount} attempts; accepted {len(geometries)}/{n_samples} geometries")
    if debug > 0: print(f"GEOM_SAMPLING_FROM_VNM: Accepted {len(geometries)}/{n_samples} geometries in {count} attempts")
    return np.asarray(geometries, dtype=float).reshape((-1, N_atoms, 3)), np.asarray(q_coords, dtype=float).reshape((-1, N_modes)), np.asarray(energies, dtype=float)

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
