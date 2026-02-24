####################################################
#### Functions that are Specific to Azo Species ####
####################################################

import numpy as np
import scope.constants as Constants

from scope_azo.azo_classes              import *
from copy                               import deepcopy
from scope.geometry                     import *
from scope.connectivity                 import *
from scope.operations.vecs_and_mats     import gaussian
from scope.elementdata          import ElementData
elemdatabase = ElementData()

def get_3D(smiles, debug: int=0):
    import warnings
    from scope.geometry import centercoords

    try:
        from openbabel import pybel as pb
    except Exception as exc:
        warnings.warn("Open Babel is not installed. Install it with: conda install -c conda-forge openbabel", RuntimeWarning,stacklevel=2)

    mol = pb.readstring('smiles', smiles) # first argument is a string format
    mol.addh() # add Hs for 3D
    mol.make3D()
    labels = []
    coord = []
    for at in mol:
        labels.append(elemdatabase.elementsym[at.atomicnum])
        coord.append(list([at.coords[0], at.coords[1], at.coords[2]]))
    coord = centercoords(coord, 0) 
    return labels, coord

def combine_smiles(lefts: list[str], rights: list[str], subs: list[str], systems: list['System_azo'] = None, debug: int = 0) -> list['System_azo']:
    """
    Combines left, right and substituent SMILES strings to create azo compounds. Geometries are named following the 
    order of the input lists. 
    It is hardly recommended to use this function for SMILES that do follow the following structure: c1(ccccc1)/N=N/ring2. 
    This allows a well-defined indices to set dihedral angles.

    Guide
    -----
    Selected atoms to rotate dihedral angles are  the following (example for azobenzene):

      4 --- 5                     9 --- 10
    //      \\\\                  //      \\\\
    3          0 --- 6 === 7 --- 8        11
    \\        /      N     N     \\       /
      2 === 1                     13 === 12

    where:
    - at0 = 0: is the atom with index 0, found in the left ring
    - at1 = 1: is the atom with index 1, found in the left ring
    - at2 = 6: is the Nitrogen atom with index 6, found in the Azo fragment
    - at3 = 7: is the Nitrogen atom with index 7, found in the Azo fragment
    - at4 = 8: is the atom with index 8, found in the right ring
    - at5 = 9: is the atom with index 9, found in of the right ring

    The variables at1, at2, at3 and at4 are used to define the dihedral angle of the azo fragment.

    The variables at0 and at5 are used to define the dihedral angle of the rings with
    respect to the azo fragment.

    """
    if systems is None:
        if debug != 0: print(f"AZOS.COMBINE_SMILES: systems is None, creating empty list.")
        systems = []
    else:
        if debug != 0: print(f"AZOS.COMBINE_SMILES: systems is not None, using existing list. Actual size: {len(systems)}")

    from scope_azo.azo_classes import System_azo

    existing_names = {sys.name for sys in systems}

    # SMILES combinations
    # Note: Assuming 'SUB' is only present in the left fragment.
    core_template = '(LEFT)/N=N/(RIGHT)'
    print(f"AZOS.COMBINE_SMILES: START combining {len(lefts)} left fragments, {len(rights)} right fragments, and {len(subs)} subs")
    for idx, left in enumerate(lefts):
        for jdx, right in enumerate(rights):
            for kdx, sub in enumerate(subs):
                # Combination of smiles
                left_fragment = left.replace("SUB",sub)            # Replaces SUB with i substituent
                current_smiles = core_template.replace("(LEFT)", left_fragment).replace("(RIGHT)", right)
                name = str(f"{idx}_{jdx}_{kdx}")                   # Names in (left_right_sub) order 

                try:
                    if name in existing_names:
                        print(f'Skipping {name}, already exists')
                        continue

                    new_system = System_azo(name, current_smiles)

                    new_system.create_trans(debug=debug)
                    new_system.create_cis(debug=debug)
                    new_system.create_ts(debug=debug)

                    if new_system.find_source('trans') is not None and new_system.find_source('cis') is not None:
                        systems.append(new_system)
                        existing_names.add(name)
                    else:
                        print(f"AZOS.COMBINE_SMILES: ERROR combining {name}: Trans or cis isomer not found")

                except Exception as e:
                    print(f"AZOS.COMBINE_SMILES: ERROR combining {name}: {e}")

    print(f"AZOS.COMBINE_SMILES: END. Total valid systems: {len(systems)}")
    return systems

################################
###     HELP FUNCTION
################################

def solve_dihedral(labels, coord, at0, at1, at2, at3, at4, at5, adjmat_ref, adjnum_ref, debug: int=0):
    """
    Solves the dihedral angles of a system to avoid steric hindrance. 
    It does so by rotating adjacent rings dihedral angles: 
        1 - left ring (at0-at3), 
        2 - right ring (at2-at5)

    Combinations are done in a grid of 64x64 steps, from -180 to 180 degrees.
    The first valid combination is returned. 
    """
    from itertools                  import product
    rot_steps = np.linspace(-180,180, 64).astype(int)
    rot_combinations = list(product(rot_steps, rot_steps))
    fixed_collision = False
    for angle1, angle2 in rot_combinations:
        coord_try = set_dihedral(labels, coord, angle1, at0,at1,at2,at3, adjmat=adjmat_ref, adjnum=adjnum_ref)  # Coords de tsinv                        
        coord_try = set_dihedral(labels, coord_try, angle2, at2,at3,at4,at5, adjmat=adjmat_ref, adjnum=adjnum_ref)
        _, adjmat_try, adjnum_try = get_adjmatrix(labels,coord_try)
        found_fix = np.array_equal(adjmat_try, adjmat_ref) and np.array_equal(adjnum_try, adjnum_ref)
        if found_fix:
            coord = coord_try
            if debug != 0: print(f'AZOS.SOLVE_DIHEDRAL: Found good geometry by rotating adjacent dihedrals')
            fixed_collision = True
            return fixed_collision, coord

def rotation_matrix(axis: np.ndarray, angle: float):
    """
    Returns a rotation matrix for a given axis and angle based on Rodrigues' rotation formula.

    Parameters
    ----------
    axis : array-like
        A 3D vector representing the axis of rotation.
    angle : float
        The angle of rotation in radians.
    
    Returns
    -------
    R : array-like
        A 3x3 rotation matrix.
    """
    axis = axis / np.linalg.norm(axis)  
    x, y, z = axis
    c, s = np.cos(angle), np.sin(angle)
    C = 1 - c
    R = np.array([
        [c + x*x*C,     x*y*C - z*s,   x*z*C + y*s],
        [y*x*C + z*s,   c + y*y*C,     y*z*C - x*s],
        [z*x*C - y*s,   z*y*C + x*s,   c + z*z*C]])
    return R
    
def put_atoms_on_xy(coord, atom1:int, atom2:int, atom3:int, debug=0):
    """
    Applies rotations to put atom1-atom2 bond on x-axis and atom3 on xy plane.

    Parameters
    ----------
    coord : array-like
        Array of shape (N, 3) representing the coordinates of N atoms. 
    atom1, atom2, atom3 : int
        Indices of the three atoms to be aligned.

    Returns
    -------
    rotated_coords : array-like
        Array of shape (N, 3) representing the rotated coordinates on xy plane.

    Notes
    -----
    The function first translates the coordinates so that central atom (atom2) is at the origin.
    
    """

    c1 = centercoords(coord, atom2)


    # Remove y component of atom1-atom2 bond
    x_axis = np.array([1., 0, 0])  # X-axis vector
    y_axis = np.array([0, 1., 0])  # Y-axis vector

    v1  = c1[atom1]                     # Vector from centre to atom1
    rot_axis1 = np.cross(v1, x_axis)    # Axis of rotation
    a1 = get_angle(v1, x_axis)

    # Check if the bond is already aligned with the x-axis
    if np.linalg.norm(a1) < 1e-8:
        if np.dot(v1, x_axis) < 0:
            axis = np.array([0,1,0])
            angle = np.pi
        else:
            angle = 0
    else:
        a1 = get_angle(v1, x_axis)

    R = rotation_matrix(rot_axis1, a1)  
    c2 = c1 @ R.T                               # Apply rotation to c1
    v2 = c2[atom3]                              # Vector from centre to atom2
    y_proj = np.array([0, v2[1], v2[2]])        # Projection of v2 y-axis
    a2 = get_angle(y_proj, y_axis)               # a2(y_projected, y_axis)
    if v2[2] < 0: a2 = -a2                      # Correct angle if z component is negative
    c3 = rot_in_x(a2, c2)                       # Rotate c2 around the x-axis by the angle a2
    if debug == 1: print('Angle between y_proj and y_axis',a2 * (180 / np.pi))
    return c3

################################
#
#  Functions for thermal props
#
################################

def compute_t(g_ts:float, g_iso:float, T:float=298.15):
    '''
    Computes half-life time in seconds using Eyring equation.
    
    Parameters
    ----------
    g_ts : float
        Free Gibbs energy of transition state in Hartree
    g_iso : float
        Free Gibbs energy of isomer/conformer in Hartree
    T : float
        Temperature in Kelvin
    
    Returns
    -------
    t : float
        half-life time in seconds
    k : float
        rate constant in s^-1
    '''
    k_b = Constants.boltz_J # J/K
    h = Constants.planck_Js # J·s
    R = Constants.R_J       # 8.31 J/(K·mol)

    dG = (g_ts - g_iso)* Constants.har2kJmol * 1000  # in J/mol
    k = ((k_b*T )/ h)* np.exp(-dG / (R * T))  
    t = np.log(2) / k           # Assuming a first-order reaction
    return float(t), float(k)

def correct_tripletG(isomer_state, triplet_state, T:float=298.15, overwrite = False, p_sh:float = 0.0002, debug: int=0):
    '''
    Corrects the Gtot of a triplet Molecule_azo object using the Gtot of the parent Molecule_azo object. 
    Correction is done considering the increase of energy due to surface hopping between the singlet and triplet PESs. 

    Parameters
    ----------
    triplet_specie : Molecule_azo
        The triplet Molecule_azo object to correct the Gtot of.
    T : float, optional
        The temperature in Kelvin. The default is 298.15 K.
    overwrite : bool, optional
        Whether to overwrite the existing Gtot_corr value. The default is False.
    p_sh : float, optional
        The probability of surface hopping. The default is 0.0002.
    debug : int, optional
        The debug level. The default is 0
    '''
    from scope.classes_data import Data
    k_b = Constants.boltz_J # J/K
    h = Constants.planck_Js # J·s
    R = Constants.R_J       # 8.31 J/(K·mol)

    triplet_source = triplet_state._source
    triplet_system = triplet_source._sys

    isomer_source = isomer_state._source

    exist = 'Gtot_corr' in triplet_state.results

    if not 'Gtot' in isomer_state.results or not 'Gtot' in triplet_state.results:
        print(f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: No Gtot found for {isomer_source.name} or {source_parent.name}.')
        return
    
    if not exist or overwrite:
        if debug > 0: print(f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: Found Gtot for {isomer_source.name} and {triplet_source.name}. Correcting Gtot of triplet State.')
        G_triplet = triplet_state.results['Gtot'].value
        G_iso = isomer_state.results['Gtot'].value

        if debug > 0: print(f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: {triplet_source.name} Gtot: {G_triplet} hartree, iso Gtot: {G_iso} hartree')

        dG = (G_triplet - G_iso) * Constants.har2kJmol * 1000  # in J/mol
        t, k = compute_t(G_triplet, G_iso, T)       # Adiabatic rate
        k_sh = k * p_sh                             # Surface hopping rate (non-adiabatic rate)
        deltax = - (dG + R * T * np.log((h*k_sh)/(k_b*T))) # in J/mol

        if debug > 0: print(f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: Adding deltax in kcal/mol: {deltax*0.24/1000}')
        newG = (G_triplet + deltax / (1000*Constants.har2kJmol))
        newG = float(newG)
        newG_data = Data("Gtot_corr", newG, "au", "correct_tripletG")
        
        triplet_state.add_result(newG_data, overwrite=overwrite)
        if debug > 0: print (f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: Corrected Gtot of {triplet_source.name} triplet state by {deltax*0.24/1000:.2f} Kcal/mol.')
        return newG_data

def check_triplet(input_object, overwrite:bool = False, debug:int = 0):
    """
    Checks if Triplet energies from results after computations have been corrected.
    If not, correction is applied.
    """
    if input_object.type == 'specie':
        system = input_object._sys
    elif input_object.type == 'state':
        system = input_object._source._sys
    else: 
        raise Exception('CHECK_TRIPLET: Type not recognised. Please insert a Specie or System-derived object')
    
    source_exists, trans = system.find_source('trans')
    state_exists, trans_state = trans.find_state("opt")


    triplet_indices = [i for i, source in enumerate(system.sources) if source.spin==2]
    #Check if the computation exists 

    good=0
    if source_exists and state_exists:
        for idx in triplet_indices:
            exist_state, triplet_state = system.sources[idx].find_state('opt')
            if exist_state:
                correct_tripletG(trans_state, triplet_state, overwrite=overwrite, debug=debug)
                if "Gtot_corr" in triplet_state.results.keys(): good+=1
    if good == len(triplet_indices): return True
    return False


def calculate_dG(t): 
    k_b = Constants.boltz_J # J/K
    T = 298.15              # K    
    h = Constants.planck_Js # J·s
    R = 8.31446             # J/(K·mol)

    k = np.log(2) / t          # Assuming a first-order reaction
    # dG = (g_ts - g_iso)* Constants.har2kJmol * 1000  # in J/mol
    # k = ((k_b*T )/ h)* np.exp(-dG / (R * T))  
    dG = -R * T * np.log((k * h) / (k_b * T))  # in J/mol
    dG /= 1000  # in kJ/mol
    dG *= 0.24  # in kcal/mol
    return dG
    
def show_thermal_data(systems):
    for sys in systems:
        cis = sys.find_source('cis')[1]
        trans = sys.find_source('trans')[1]
        cis_opt = cis.find_state('opt')
        trans_opt = trans.find_state('opt')
        # dG = (cis_opt[1].results['Gtot'].value - trans_opt[1].results['Gtot'].value)* har2kJmol * 0.24
        if cis_opt[0] and trans_opt[0]:
            print(f'---------------Azo: {sys.name}-----------------')
            # print(f'Delta G (cis - trans) = {dG:.2f} kcal/mol')
            print('                 CIS                         ')
            print(cis_opt[1].results['dG_cross'])
            print("t0.5 = " + format_time(cis_opt[1].results['halflife'].value))
            print(f'by {cis.mets}')
            print('                TRANS                        ')
            print(trans_opt[1].results['dG_cross'])
            print("t0.5 = " + format_time(trans_opt[1].results['halflife'].value))

            print(f'by {trans.mets}\n')
            print('\n')

def format_time(t):
    """
    Converts seconds to a more comprehensive unit.
    """
    units = [
        (1e-15, "fs", 1e15),
        (1e-12, "ps", 1e12),
        (1e-9,  "ns", 1e9),
        (1e-6,  "µs", 1e6),
        (1e-3,  "ms", 1e3),
        (1,     "s",  1),
        (60,    "min", 1/60),
        (3600,  "h",  1/3600),
        (86400, "d",  1/86400),
        (31557600, "y", 1/31557600),        # y
        (31557600000, "ky", 1/31557600000)  # millenia
    ]
    for threshold, unit, factor in units:
        if abs(t) < threshold * 100:
            return f"{t * factor:.2f} {unit}"
    return f"{t / 31557600000:.2f} ka"


#######################################
#
#         OPTICAL PROPERTIES
#             FUNCTIONS
#
#######################################

def get_PSS(cis_state, trans_state, lamp : "Lamp", phi_EZ = 0.3, phi_ZE = 0.5, t_EZ=None, t_ZE=None, debug=0):

    cis_source = cis_state._source
    trans_source = trans_state._source
    system = cis_source._sys
    sigma_, = build_spectrum(cis_state, normalize=False, units=True) # Need units to compute PSS
    sigma_E  = build_spectrum(trans_state, normalize=False, units=True) # Need units to compute PSS

    # Photon flux from lamp
    photon_flux = get_photon_flux_spectrum(lamp.wavelength, lamp.fwhm, lambda_grid, Itot=lamp.irradiance)
    
    # Half-lives in seconds
    if not 'halflife' in trans_state.results.keys():
        print('No halflife found for trans isomer. Computing...')
        trans_source.set_halflife_time(skip_triplets=False, overwrite=False)
    if not 'halflife' in cis_state.results.keys():
        print('No halflife found for cis isomer. Computing...')
        cis_source.set_halflife_time(skip_triplets=False, overwrite=True)
        raise ValueError(f'Error: No halflife found for Z or E isomers in system {system.name}')

    # If no half-life times are given (for research), halflife is taken from the corresponding state
    if t_EZ is None:
        t_EZ = trans_state.results['halflife'].value
    if t_ZE is None:
        t_ZE = cis_state.results['halflife'].value

    # Photochemical rates
    k_ph_EZ = phi_EZ * np.trapezoid(sigma_E * photon_flux, lambda_grid)
    k_ph_ZE = phi_ZE * np.trapezoid(sigma_Z * photon_flux, lambda_grid)
    # Thermal rates
    k_th_EZ = np.log(2) /  t_EZ
    k_th_ZE = np.log(2) /  t_ZE
    # Steady-state populations
    N_E = (k_th_ZE + k_ph_ZE) / (k_ph_EZ + k_ph_ZE + k_th_EZ + k_th_ZE)
    return lambda_grid, sigma_Z, sigma_E, N_E

def build_spectrum(erange, energies, fosc, sigma=0.2, normalize=False, units=True, debug: int=0):
    """
    Builds a spectrum from a list of transitions.
    
    Parameters
    ----------
    erange: array_like             Desired energy range for the spectrum, in eV.
    energies : list of floats      List of excitation energies, in eV.
    fosc : list of floats          List of oscillator strengths for each transition.
    sigma : float, optional        Standard deviation of the Gaussian, in eV. Default is 0.2 eV.
    normalize : bool, optional     Whether to normalize the Gaussian. Default is False.
    units : bool, optional         Whether to return the spectrum in units of m^2 / molecule. Default is True.
    
    Returns
    -------
    spec : array_like
        Spectrum, in m^2 / molecule if units is True, otherwise in arbitrary units.
    """
    assert len(energies) == len(fosc), "Energies and oscillator strengths must have the same length."
    
    spec = np.zeros_like(energies)
    for E0, f in zip(energies, fosc):
        spec += f * gaussian(erange, E0, sigma=sigma, normalize=normalize)

    if units: 
        # Conversion to m^2 / molecule
        K = (np.pi * Constants.planck_Js * Constants.elem_charge) / (Constants.epsilon_0 * Constants.speed_light * Constants.electron_mass)
        return spec[::-1] * K 
    return spec[::-1]

def get_photon_flux_spectrum(lam0_nm, fwhm_nm, lam_grid, Itot, power=None, debug=0):
    """
    Returns the photon flux spectrum from a given wavelength grid and intensity.
    
    Parameters
    ----------
    lam0_nm : float
        Central wavelength, in nm.
    fwhm_nm : float
        Full width at half maximum, in nm.
    lam_grid : array_like
        Wavelength grid, in nm.
    Itot : float
        Total intensity, in W/m2/nm.
    power : float, optional
        Power, in W. Default is None.
    debug : int, optional
        Debug level. Default is 0.
    
    Returns
    -------
    I_lambda : array_like
        Photon flux spectrum, in photons m-2 s-1 nm-1.
    """

    sigma = fwhm_nm / (2 * np.sqrt(2 * np.log(2)))
    profile = gaussian(lam_grid, lam0_nm, sigma=sigma)
    if power is not None:
        area = np.pi * (4.605e-3)**2  # in m2, area of a circle with diameter 0.92 cm
        I_lambda = power * profile / area   # W/m2/nm
    else:
        I_lambda = Itot * 1e-3 / 1e-6 * profile  # mW/mm2 to W/m2/nm
    # I_lambda = Itot * profile  # W/m2/nm
    lam_m = np.ones_like(lam_grid)* lam_grid * 1e-9  # in m
    photon_E = Constants.planck_Js * Constants.speed_light / lam_m  # in J
    phi = I_lambda / photon_E           # photons m-2 s-1 nm-1        
    return phi

def build_pss_spectrum(pss, sigma_Z, sigma_E, debug=0):
    """
    Builds the spectrum using 
    """
    # sigma_Z *= Constants.avogadro / (1000 * np.log(10)) * 1e4  # in M^-1 cm^-1
    # sigma_E *= Constants.avogadro / (1000 * np.log(10)) * 1e4  # in M^-1 cm^-1
    sigma_pss = pss * sigma_E + (1 - pss) * sigma_Z
    return sigma_pss

def wavelength_to_rgb(nm: float):
    """Approximate the RGB color perceived for a wavelength in nm."""
    gamma = 0.8
    # intensity = 1.0
    nm +=30
    if nm < 360 or nm > 780:
        return (0.5, 0.5, 0.5)  # gray for non-visible (e.g., 350 nm)

    if 360 <= nm < 440:
        r, g, b = -(nm - 440) / (440 - 380), 0.0, 1.0
    elif 440 <= nm < 490:
        r, g, b = 0.0, (nm - 440) / (490 - 440), 1.0
    elif 490 <= nm < 510:
        r, g, b = 0.0, 1.0, -(nm - 510) / (510 - 490)
    elif 510 <= nm < 580:
        r, g, b = (nm - 510) / (580 - 510), 1.0, 0.0
    elif 580 <= nm < 645:
        r, g, b = 1.0, -(nm - 645) / (645 - 580), 0.0
    else:  # 645–780
        r, g, b = 1.0, 0.0, 0.0

    # intensity correction at spectrum edges
    if 360 <= nm < 420:
        factor = 0.3 + 0.7 * (nm - 380) / (420 - 380)
    elif 420 <= nm <= 700:
        factor = 1.0
    else:
        factor = 0.3 + 0.7 * (780 - nm) / (780 - 700)

    r = (r * factor) ** gamma if r > 0 else 0.0
    g = (g * factor) ** gamma if g > 0 else 0.0
    b = (b * factor) ** gamma if b > 0 else 0.0

    return (r, g, b)

