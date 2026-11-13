"""Finite differences vibrational analysis.

Coordinates are in Angstrom; Cartesian Hessians are in Hartree/Bohr^2 and
atomic masses in amu.
"""
from scope import constants
import numpy as np


def apply_coord_displacement(coord, atom: int, axis: int, displacement: float=0.01, units: str='angstrom'):
    """Returns Angstrom coordinates with one zero-based atom/axis displaced.

    "units" describes the displacement units, not those of the "coord".
    The input is unchanged; the returned array has shape (natoms, 3).
    """
    # Check inputs
    if np.iscomplexobj(coord):
        raise ValueError("APPLY_COORD_DISPLACEMENT: Coordinates must be real")
    mod_coord = np.array(coord, dtype=float, copy=True)
    if mod_coord.ndim != 2 or mod_coord.shape[1] != 3 or len(mod_coord) == 0 or not np.all(np.isfinite(mod_coord)):
        raise ValueError("APPLY_COORD_DISPLACEMENT: Expected finite (natoms, 3) coordinates")
    if not isinstance(atom, (int, np.integer)) or isinstance(atom, (bool, np.bool_)) or not 0 <= atom < len(mod_coord):
        raise ValueError("APPLY_COORD_DISPLACEMENT: Atom index is out of range")
    if not isinstance(axis, (int, np.integer)) or isinstance(axis, (bool, np.bool_)) or axis not in [0,1,2]:
        raise ValueError("APPLY_COORD_DISPLACEMENT: Axis must be 0, 1 or 2")
    if not np.isfinite(displacement):
        raise ValueError("APPLY_COORD_DISPLACEMENT: Displacement must be finite")

    # Apply the displacement in Angstrom, converting from Bohr if requested
    if units == 'bohr':
        displacement *= constants.bohr2angs
    elif units != 'angstrom':
        raise ValueError("APPLY_COORD_DISPLACEMENT: units must be 'bohr' or 'angstrom'")
    mod_coord[atom,axis] += displacement
    return mod_coord


def get_central_difference(f1, f2, displacement: float=0.01, units: str='bohr', force_units: str='hartree/bohr', debug: int=0):
    """Returns -(F_plus - F_minus)/(2*displacement) in Hartree/Bohr^2.

    f1 and f2 are forces at positive and negative displacements, respectively, with matching shapes;
    displacement is the positive one-sided step;
    units is 'bohr' (default) or 'angstrom';
    force_units is 'hartree/bohr' or 'ry/bohr'; select 'ry/bohr' for forces returned by the QE parser.
    The returned array has the force-array shape; inputs are not modified.
    """
    # Check inputs
    if np.iscomplexobj(f1) or np.iscomplexobj(f2):
        raise ValueError("GET_CENTRAL_DIFFERENCE: Forces must be real")
    f1 = np.asarray(f1, dtype=float)
    f2 = np.asarray(f2, dtype=float)
    if f1.shape != f2.shape or f1.size == 0 or not np.all(np.isfinite(f1)) or not np.all(np.isfinite(f2)): raise ValueError("GET_CENTRAL_DIFFERENCE: Forces must have matching finite nonempty shapes")
    if not np.isfinite(displacement) or displacement <= 0:
        raise ValueError("GET_CENTRAL_DIFFERENCE: Displacement must be finite and positive")
    if units == 'angstrom': displacement *= constants.angs2bohr
    elif units != 'bohr':
        raise ValueError("GET_CENTRAL_DIFFERENCE: units must be 'bohr' or 'angstrom'")

    # Convert the force difference to Hartree/Bohr^2, if necessary
    if force_units == 'ry/bohr':
        factor = constants.ry2har
    elif force_units == 'hartree/bohr':
        factor = 1.0
    else:
        raise ValueError("GET_CENTRAL_DIFFERENCE: force_units must be 'hartree/bohr' or 'ry/bohr'")
    if debug > 0:
        print(f"GET_CENTRAL_DIFFERENCE: step={displacement:g} Bohr, forces={force_units}")

    # Returns the central-difference Hessian entries in Hartree/Bohr^2
    return -(f1-f2) * factor / (2*displacement)


def findiff_displacements(coord, displacement: float=0.01, units: str='bohr'):
    """Returns the 6N displaced Angstrom geometries and their coordinate names.

    coord is a finite (natoms, 3) Angstrom array. displacement is a positive
    one-sided step in units ('bohr' or 'angstrom'). Each zero-based atom/axis
    produces a positive then a negative geometry, named _atom_axis_pos/neg.
    """

    # Check inputs
    if np.iscomplexobj(coord):
        raise ValueError("FINDIFF_DISPLACEMENTS: Coordinates must be real")
    coord = np.asarray(coord, dtype=float)
    if coord.ndim != 2 or coord.shape[1] != 3 or len(coord) == 0 or not np.all(np.isfinite(coord)):
        raise ValueError("FINDIFF_DISPLACEMENTS: Expected finite (natoms, 3) coordinates")
    if not np.isfinite(displacement) or displacement <= 0:
        raise ValueError("FINDIFF_DISPLACEMENTS: Displacement must be finite and positive")
    if units not in ['bohr','angstrom']:
        raise ValueError("FINDIFF_DISPLACEMENTS: units must be 'bohr' or 'angstrom'")
    geoms, names = [], []
    for atom in range(len(coord)):
        for axis in range(3):
            # Positive Displacement
            geoms.append(apply_coord_displacement(coord, atom, axis, displacement, units))
            names.append(f'_{atom}_{axis}_pos')
            # Negative Displacement
            geoms.append(apply_coord_displacement(coord, atom, axis, -displacement, units))
            names.append(f'_{atom}_{axis}_neg')

    # Returns the 6N displaced geometries and their names
    return geoms, names


def mass_weight_hessian(hessian, masses, debug: int=0):
    """Divides Cartesian Hessian entries by sqrt(m_i*m_j), without modifying inputs.

    hessian is a finite real (3N, 3N) array in Hartree/Bohr^2, ordered by atom
    then x/y/z. masses contains N finite positive masses in amu. The result is
    in Hartree/(Bohr^2*amu); it is not an electron-mass atomic-unit matrix.
    """
    # Check inputs
    if np.iscomplexobj(hessian) or np.iscomplexobj(masses):
        raise ValueError("MASS_WEIGHT_HESSIAN: Hessian and masses must be real")
    hessian = np.asarray(hessian, dtype=float)
    masses  = np.asarray(masses, dtype=float)
    if masses.ndim != 1 or len(masses) == 0 or not np.all(np.isfinite(masses)) or np.any(masses <= 0):
        raise ValueError("MASS_WEIGHT_HESSIAN: Masses must be a finite positive vector")
    if hessian.shape != (3*len(masses),3*len(masses)) or not np.all(np.isfinite(hessian)):
        raise ValueError("MASS_WEIGHT_HESSIAN: Expected a finite (3N, 3N) Hessian")

    # Apply the mass weighting by dividing each Hessian entry by sqrt(m_i*m_j)
    mass_components = np.repeat(masses, 3)
    if debug > 0:
        print(f"MASS_WEIGHT_HESSIAN: {len(masses)} atoms, masses in amu")

    # Returns the mass-weighted Hessian in Hartree/(Bohr^2*amu)
    return hessian / np.sqrt(mass_components[:,None]*mass_components[None,:])


def project_out(hessian, proj_rot: bool=False, proj_tra: bool=True, debug: int=0):
    """Projects uniform translations from a Cartesian Hessian before mass weighting.

    Returns P.T @ hessian @ P, where P removes the three rigid translations.
    With proj_tra=False, returns an independent copy. Rotational projection needs
    geometry information and is not implemented; proj_rot=True raises an error.
    hessian must be finite, real, square and ordered as atom-wise x/y/z.
    """
    # Check inputs
    if proj_rot:
        raise NotImplementedError("PROJECT_OUT: Rotational projection is not implemented")
    if np.iscomplexobj(hessian):
        raise ValueError("PROJECT_OUT: Hessian must be real")
    hessian = np.asarray(hessian, dtype=float)
    if hessian.ndim != 2 or hessian.shape[0] == 0 or hessian.shape[0] != hessian.shape[1] or hessian.shape[0] % 3 or not np.all(np.isfinite(hessian)):
        raise ValueError("PROJECT_OUT: Expected a finite (3N, 3N) Hessian")
    if not proj_tra: return hessian.copy()
    # Project out the three Cartesian translations with a simple projector
    natoms       = len(hessian)//3
    translations = np.tile(np.eye(3), (natoms,1))/np.sqrt(natoms)
    projector    = np.eye(3*natoms) - translations@translations.T
    if debug > 0:
        print("PROJECT_OUT: Removing three Cartesian translations")
    # Return the projected Hessian in Hartree/Bohr^2, without modifying the input
    return projector.T @ hessian @ projector


def get_VNM_from_hessian(hessian, labels: list, masses: list, proj_rot: bool=False, proj_tra: bool=True, discard_translations: bool=False, debug: int=0):
    """Returns mass-weighted VNMs from a Cartesian Hessian in Hartree/Bohr^2.

    labels contains elemental symbols in Hessian atom order; masses supplies the
    calculation's masses in amu. Symmetrizes, optionally projects translations in
    Cartesian space, mass-weights and diagonalizes with a symmetric eigensolver.
    Rotational projection is not implemented. Inputs are not modified.
    With discard_translations=True, requires proj_tra=True and diagonalizes only
    the translation-free subspace, returning 3N-3 modes. Otherwise returns 3N.
    Raw relative Hessian antisymmetry above 1% prints a warning, independently
    of debug. This warning does not prevent diagonalization.

    Frequencies are ascending in cm^-1, with negative values for imaginary modes
    and any near-zero roundoff retained. Modes have unit mass-weighted norm and
    is_mass_weighted=True. Reduced masses are in amu and each VNM retains the
    supplied atomic masses. Intensities and spectroscopic force constants are
    not calculated. This function does not collect or register Job results.
    """
    from scope.classes_qc import VNM
    from scope.connectivity import labels2atnums
    from scope.operations.vecs_and_mats import symmetrize

    # 0) Check atom order and the Cartesian Hessian dimensions
    if np.iscomplexobj(hessian):
        raise ValueError("GET_VNM_FROM_HESSIAN: Hessian must be real")
    hessian = np.asarray(hessian, dtype=float)
    masses  = np.asarray(masses)
    natoms  = len(labels)
    if natoms == 0 or hessian.shape != (3*natoms,3*natoms):
        raise ValueError("GET_VNM_FROM_HESSIAN: Labels and Hessian dimensions must agree")
    atnums = labels2atnums(labels)

    if discard_translations and not proj_tra: raise ValueError("GET_VNM_FROM_HESSIAN: Discarding translations requires proj_tra=True")
    norm      = np.linalg.norm(hessian)
    asymmetry = np.linalg.norm(hessian-hessian.T)/norm if norm else 0.0
    if asymmetry > 0.01: print(f"GET_VNM_FROM_HESSIAN: WARNING: Hessian asymmetry {asymmetry:.2%} exceeds 1%")

    # 1) Symmetrize and project before applying the calculation's atomic masses
    sym_hessian  = symmetrize(hessian)
    proj_hessian = project_out(sym_hessian, proj_rot=proj_rot, proj_tra=proj_tra, debug=max(debug-1,0))
    mw_hessian   = mass_weight_hessian(proj_hessian, masses, debug=max(debug-1,0))
    masses       = np.asarray(masses, dtype=float)

    # 2) Diagonalize and convert eigenvalues to signed wavenumbers
    if discard_translations:
        # Mass-weighted translations span the null space of the projected matrix.
        translations = np.tile(np.eye(3),(natoms,1))*np.sqrt(np.repeat(masses,3))[:,None]
        basis        = np.linalg.qr(translations,mode='complete')[0][:,3:]
        eigenval, eigenvec = np.linalg.eigh(basis.T@mw_hessian@basis)
        eigenvec = basis@eigenvec
    else:
        eigenval, eigenvec = np.linalg.eigh(mw_hessian)
    factor_cminv      = 1/(np.sqrt(constants.amu2au)*constants.cm2har)
    frequencies       = np.sign(eigenval)*np.sqrt(np.abs(eigenval))*factor_cminv
    cartesian_modes   = eigenvec/np.sqrt(np.repeat(masses,3))[:,None]
    reduced_masses    = 1/np.sum(cartesian_modes**2, axis=0)

    # 3) Store the mass-weighted eigenvectors without substituting element-table masses
    atomidxs = list(range(1,natoms+1))
    VNMs = []
    for idx, freq in enumerate(frequencies):
        mode = eigenvec[:,idx].reshape(natoms,3)
        vnm = VNM(idx+1, float(freq), red_mass=float(reduced_masses[idx]))
        vnm.set_mode(atomidxs, atnums, mode[:,0], mode[:,1], mode[:,2], is_mass_weighted=True)
        vnm.masses = masses.tolist()
        VNMs.append(vnm)
    if debug > 0:
        print(f"GET_VNM_FROM_HESSIAN: {len(VNMs)} modes; relative Hessian antisymmetry={asymmetry:.3e}")
    return VNMs


def get_VNM_from_findiff(job: object, proj_rot: bool=False, proj_tra: bool=True, debug: int=0):
    """Collects periodic QE forces and returns 3N-3 VNMs when projecting translations.

    Selects the latest (step, run_number) for every _atom_axis_pos/neg keyword,
    independent of list order. Returns None while any point is missing, pending,
    unsuccessful or lacks forces. Invalid stencils raise ValueError. Inputs must
    match the reference State's atom order, cell and symmetric displacements.
    Uses submitted atomic masses and actual rounded steps, converting QE forces
    from Ry/Bohr. A 'ref' computation, when present, must supply converged energy
    and forces at the undisplaced geometry; Job registration requires it. The
    returned VNMs exclude translations when proj_tra=True. The Job and its States
    are not modified; rotations are unsupported.
    """
    from pathlib import Path
    from scope.read_write import read_lines_file
    from scope.software.quantum_espresso.qe_parse import parse_findiff_input

    # 0) Identify the reference geometry and the latest attempt at each point
    if proj_rot: raise NotImplementedError("GET_VNM_FROM_FINDIFF: Rotational projection is not implemented")
    found, reference = job.source.find_state(job.istate)
    if not found or not hasattr(reference, 'coord'): raise ValueError("GET_VNM_FROM_FINDIFF: Reference State has no geometry")
    if job.source.object_type != 'cell': raise ValueError("GET_VNM_FROM_FINDIFF: Only fixed-cell periodic QE calculations are supported")
    step  = getattr(job.job_data, 'findiff_displacement', 0.01)
    units = getattr(job.job_data, 'findiff_units', 'bohr')
    geoms, names = findiff_displacements(reference.coord, step, units)
    if any(comp.keyword == 'ref' for comp in job.computations):
        geoms.append(np.asarray(reference.coord,dtype=float))
        names.append('ref')
    latest = {}
    for comp in job.computations:
        if comp.keyword not in names: raise ValueError(f"GET_VNM_FROM_FINDIFF: Unexpected displacement keyword {comp.keyword!r}")
        previous = latest.get(comp.keyword)
        if previous is not None and (comp.step,comp.run_number) == (previous.step,previous.run_number): raise ValueError(f"GET_VNM_FROM_FINDIFF: Duplicate attempt for {comp.keyword}")
        if previous is None or (comp.step,comp.run_number) > (previous.step,previous.run_number): latest[comp.keyword] = comp
    if len(latest) != len(names):
        if debug > 0: print(f"GET_VNM_FROM_FINDIFF: {len(latest)}/{len(names)} displacement points available")
        return None

    # 1) Check submitted geometries and retain forces from the current attempts
    natoms = len(reference.labels)
    coords, forces, masses = [], [], None
    for name, expected in zip(names, geoms):
        comp = latest[name]
        if comp.software != 'qe' or comp.type != 'scf': raise ValueError("GET_VNM_FROM_FINDIFF: Displacements must be QE SCF computations")
        valid = comp.isregistered and getattr(comp,'isfinished',False) and getattr(comp,'isgood',False) and not comp.has_update
        if not valid or not Path(comp.out_path).is_file() or getattr(comp,'forces',None) is None:
            if debug > 0: print(f"GET_VNM_FROM_FINDIFF: Waiting for converged forces at {name}, run {comp.run_number}")
            return None
        if name == 'ref' and (getattr(comp,'energy',None) is None or not np.isfinite(comp.energy)):
            if debug > 0: print("GET_VNM_FROM_FINDIFF: Waiting for the reference SCF energy")
            return None
        if not Path(comp.inp_path).is_file(): raise ValueError(f"GET_VNM_FROM_FINDIFF: Missing submitted input {comp.inp_path}")
        labels, coord, atomic_masses, cell = parse_findiff_input(read_lines_file(comp.inp_path), natoms, debug=max(debug-1,0))
        coord         = np.asarray(coord, dtype=float)
        atomic_masses = np.asarray(atomic_masses, dtype=float)
        force         = np.asarray(comp.forces, dtype=float)
        if labels != list(reference.labels) or not np.allclose(coord, expected, rtol=0, atol=1e-6): raise ValueError(f"GET_VNM_FROM_FINDIFF: Geometry or atom order mismatch at {name}")
        if not np.all(np.isfinite(cell)) or not np.allclose(cell, reference.cell_vector, rtol=0, atol=1e-6): raise ValueError(f"GET_VNM_FROM_FINDIFF: Cell mismatch at {name}")
        if not np.all(np.isfinite(atomic_masses)) or np.any(atomic_masses <= 0): raise ValueError(f"GET_VNM_FROM_FINDIFF: Invalid atomic masses at {name}")
        if masses is not None and not np.array_equal(masses, atomic_masses): raise ValueError(f"GET_VNM_FROM_FINDIFF: Atomic masses changed at {name}")
        if force.shape != (natoms,3) or not np.all(np.isfinite(force)): raise ValueError(f"GET_VNM_FROM_FINDIFF: Invalid forces at {name}")
        coords.append(coord)
        forces.append(force)
        masses = atomic_masses

    # 2) Assemble columns using actual positive/negative input separations
    hessian = np.empty((3*natoms,3*natoms))
    for component in range(3*natoms):
        plus, minus = 2*component, 2*component+1
        atom, axis  = divmod(component,3)
        delta       = (coords[plus][atom,axis]-coords[minus][atom,axis])/2
        if delta <= 1e-6: raise ValueError("GET_VNM_FROM_FINDIFF: Step is unresolved at input coordinate precision")
        midpoint = (coords[plus]+coords[minus])/2
        if not np.allclose(midpoint, reference.coord, rtol=0, atol=1e-6): raise ValueError("GET_VNM_FROM_FINDIFF: Displacements are not centered on the reference geometry")
        hessian[:,component] = get_central_difference(forces[plus], forces[minus], delta, units='angstrom', force_units='ry/bohr').reshape(-1)
    return get_VNM_from_hessian(hessian, reference.labels, masses, proj_rot=proj_rot, proj_tra=proj_tra, discard_translations=proj_tra, debug=debug)
