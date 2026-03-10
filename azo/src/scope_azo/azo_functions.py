####################################################
#### Functions that are Specific to Azo Species ####
####################################################

import numpy as np
import scope.constants as Constants

from scope.geometry                     import centercoords
from scope.elementdata                  import ElementData
elemdatabase = ElementData()

######
def get_3D(smiles, debug: int=0):
    # Function to generate 3D coordinates from a SMILES string using Open Babel
    try:
        from openbabel import pybel as pb
    except Exception as exc:
        raise ImportError(
            "AZO_FUNCTIONS GET_3D: Open Babel (pybel) is required to generate 3D coordinates from SMILES.\n"
            "Install it with: conda install -c conda-forge openbabel\n"
        ) from exc

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

############################
#### Thermal Properties ####
############################
def eyring_halflife(g_initial: float, g_excited: float, temp: float=298.15, debug: int=0):
    '''
    Computes half-life time in seconds using Eyring equation.
    
    Parameters
    ----------
    g_initial : float       Free Gibbs energy of the ground state (in Hartree)
    g_excited : float       Free Gibbs energy of transition state (in Hartree)
    temp : float            Temperature in Kelvin
    
    Returns
    -------
    t : float               half-life time in seconds
    k : float               rate constant in s^-1
    '''
    dG = (g_excited - g_initial) * Constants.har2kJmol   # in kJ/mol
    if debug > 0: print(f'AZO.COMPUTE_T: dG: {dG} kJ/mol / {dG*Constants.kJmol2kcalmol} kcal/mol')
    dG *= 1000                                           # J/mol
    k = ((Constants.boltz_J * temp) / Constants.planck_Js) * np.exp(-dG / (Constants.R_J * temp))  
    t = np.log(2) / k  # Assuming a first-order reaction
    if debug > 0: print(f'AZO.COMPUTE_T: t05: {t:.2f} s / k: {k:.2f}')
    return float(t), float(k)

######
def calculate_dG(time: float, temp: float=298.15, debug: int=0): 
    '''
    Reverse of the Eyring equation. Gets dG (in kcal/mol) from t (in seconds)
    '''
    k = np.log(2) / time          # Assuming a first-order reaction
    dG = -Constants.R_JR * temp * np.log((k * Constants.planck_Js) / (Constants.boltz_J * temp))  # in J/mol
    dG /= 1000  # in kJ/mol
    dG *= 0.24  # in kcal/mol
    return dG

######  
def show_thermal_data(systems):                                                     ### MANEL: Esperar a definir com estan definides les funcions.
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


############################
#### Optical Properties ####
############################
 

######
def build_pss_spectrum(initial_fraction, initial_spectrum, final_spectrum, debug=0):
    """
    Builds the PSS spectrum using the fraction of trans isomer. 

    Parameters
    ----------
    initial_fraction : float    Fraction of the initial isomer.
    initial_spectrum : array    Spectrum of the initial isomer.
    final_spectrum : array      Spectrum of the final isomer.
    debug : int, optional       Verbose level.

    Returns
    -------
    array   PSS spectrum.
    """
    return initial_fraction * initial_spectrum + (1 - initial_fraction) * final_spectrum

######
def wavelength_to_rgb(nm: float):
    """Approximate the RGB color perceived for a wavelength in nm

    Parameters
    ----------
    nm : float    Wavelength in nanometers.
    
    Returns
    -------
    tuple   RGB color.
    """
    gamma = 0.8
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


###############
### For CLI ###
###############
#def combine_smiles(lefts: list[str], rights: list[str], subs: list[str], systems: list['System_azo'] = None, debug: int = 0) -> list['System_azo']:
#    """
#    Combines left, right and substituent SMILES strings to create azo compounds. Geometries are named following the 
#    order of the input lists. 
#    It is hardly recommended to use this function for SMILES that do follow the following structure: c1(ccccc1)/N=N/ring2. 
#    This allows a well-defined indices to set dihedral angles.
#
#    Guide
#    -----
#    Selected atoms to rotate dihedral angles are  the following (example for azobenzene):
#
#      4 --- 5                     9 --- 10
#    //      \\\\                  //      \\\\
#    3          0 --- 6 === 7 --- 8        11
#    \\        /      N     N     \\       /
#      2 === 1                     13 === 12
#
#    where:
#    - at0 = 0: is the atom with index 0, found in the left ring
#    - at1 = 1: is the atom with index 1, found in the left ring
#    - at2 = 6: is the Nitrogen atom with index 6, found in the Azo fragment
#    - at3 = 7: is the Nitrogen atom with index 7, found in the Azo fragment
#    - at4 = 8: is the atom with index 8, found in the right ring
#    - at5 = 9: is the atom with index 9, found in of the right ring
#
#    The variables at1, at2, at3 and at4 are used to define the dihedral angle of the azo fragment.
#
#    The variables at0 and at5 are used to define the dihedral angle of the rings with
#    respect to the azo fragment.
#
#    """
#    if systems is None:
#        if debug != 0: print(f"AZOS.COMBINE_SMILES: systems is None, creating empty list.")
#        systems = []
#    else:
#        if debug != 0: print(f"AZOS.COMBINE_SMILES: systems is not None, using existing list. Actual size: {len(systems)}")
#
#    from scope_azo.azo_classes import System_azo
#
#    existing_names = {sys.name for sys in systems}
#
#    # SMILES combinations
#    # Note: Assuming 'SUB' is only present in the left fragment.
#    core_template = '(LEFT)/N=N/(RIGHT)'
#    print(f"AZOS.COMBINE_SMILES: START combining {len(lefts)} left fragments, {len(rights)} right fragments, and {len(subs)} subs")
#    for idx, left in enumerate(lefts):
#        for jdx, right in enumerate(rights):
#            for kdx, sub in enumerate(subs):
#                # Combination of smiles
#                left_fragment = left.replace("SUB",sub)            # Replaces SUB with i substituent
#                current_smiles = core_template.replace("(LEFT)", left_fragment).replace("(RIGHT)", right)
#                name = str(f"{idx}_{jdx}_{kdx}")                   # Names in (left_right_sub) order 
#
#                try:
#                    if name in existing_names:
#                        print(f'Skipping {name}, already exists')
#                        continue
#
#                    new_system = System_azo(name, current_smiles)
#
#                    new_system.create_trans(debug=debug)
#                    new_system.create_cis(debug=debug)
#                    new_system.create_ts(debug=debug)
#
#                    if new_system.find_source('trans') is not None and new_system.find_source('cis') is not None:
#                        systems.append(new_system)
#                        existing_names.add(name)
#                    else:
#                        print(f"AZOS.COMBINE_SMILES: ERROR combining {name}: Trans or cis isomer not found")
#
#                except Exception as e:
#                    print(f"AZOS.COMBINE_SMILES: ERROR combining {name}: {e}")
#
#    print(f"AZOS.COMBINE_SMILES: END. Total valid systems: {len(systems)}")
#    return systems
