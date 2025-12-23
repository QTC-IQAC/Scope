import numpy as np
from scope import *
from scope.classes_system import *
from scope.classes_specie import *
from scope.software.gaussian.g16_parse import *
from scope.geometry import *
from scope.connectivity import *
from scope_azo.azo_functions import *

############################################
#### SYSTEM Object adapted to System_azo ###
############################################

class System_azo(System):
    def __init__(self, name: str, debug: int=0):
        System.__init__(self, name)
        self.subtype      = "system_azo"        
        self.name         = name
        self.smiles       = None
        self.fragments    = []

    ######
    def set_smiles(self, smiles):  ## En principi, el System_azo no hauria de tenir smiles associat, nomes els isomers
        self.smiles = smiles
        return self.smiles

    ######
    def get_fragments(self):
        '''Extracts the indices of the atoms of each fragment of the azo compound from the SMILES string.'''
        chars = ElementData().elementname.keys() # List of characters that define an atom. Update if needed  
        chars_lowercase = [char.lower() for char in chars]
        at_count = 1
        self.fragments = []
        block = []
        
        # New flag to track if we are inside a bracketed atom like [NH3+]
        inside_bracket = False 

        for char in self.smiles:
            if char == '/' or char == '\\': 
                self.fragments.append(block)
                block = []
            
            # CASE 1: Start of a bracketed atom (e.g. [NH3+])
            # The entire bracket group is considered as ONE atom index.
            elif char == '[':
                block.append(at_count)
                at_count += 1
                inside_bracket = True
            # CASE 2: End of bracketed atom
            elif char == ']':
                inside_bracket = False
            
            # CASE 3: Inside a bracket
            # "Symbols such as 'N', 'H', etc. are skipped because they are counted at the '['"
            elif inside_bracket:
                pass

            # CASE 4: Standard atoms outside brackets (e.g. C, N, Cl)
            elif char in chars or char in chars_lowercase:
                block.append(at_count)
                at_count += 1
            else:
                pass
        self.fragments.append(block)
        return self.fragments

    ######
    def get_dihedral_indices(self):
        '''
        Extracts the relevant indices to build molecules of azo species. It works for any azo SMILES that follows
        the following structure: c1(ccccc1)/N=N/ring2. This allows a well-defined indices to create structures.

        Parameters
        ----------
        None

        Returns
        -------
        tuple
            A tuple containing the indices of the atoms that define the dihedral angle.
                (at0, at1, at2, at3, at4, at5)

        Guide
        -----
        Indices are given in the following order (example for azobenzene):

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
        '''
        if not hasattr(self, "dihedral_indices"):
            if not hasattr(self, "fragments"):         self.get_fragments()
            fragments = self.fragments
            at0 = fragments[0][1] -1 
            at1 = fragments[0][0] -1 
            at2 = fragments[1][0] -1 
            at3 = fragments[1][1] -1 
            at4 = fragments[2][0] -1 
            at5 = fragments[2][1] -1 
            self.dihedral_indices = (at0, at1, at2, at3, at4, at5)
        return self.dihedral_indices

    ######
    def get_thermal_stability(self, target_state: str="initial", debug: int=0):
        found_cis,   cis_iso   = self.find_conformer("cis")
        found_trans, trans_iso = self.find_conformer("trans")
        if found_cis: 
            found_state_cis, state_cis = find_state(cis_iso, target_state)
            if found_state_cis:
                if hasattr(state_cis,"gs_energy"): cis_E = state_cis.gs_energy
        if found_trans:
            found_state_trans, state_trans = find_state(trans_iso, target_state)
            if found_state_trans:
                if hasattr(state_trans,"gs_energy"): trans_E = state_trans.gs_energy
        if found_cis and found_trans:
            self.dE = np.round((trans_E - cis_E) * Constants.har2kJmol/Constants.kcal2kJmol,4)
            print(self.name, self.dE, "kJ/mol. Negative means trans is more stable")

    ######
    def create_trans(self, smiles: str, debug: int = 0):
        '''
        Creates the trans structure of the azo compound from a SMILES string. 3D geometry creation is done using openbabel. 
        It sets up the trans isomer (including the creation of Specie_azo and its 'initial' State objects) and stores it as
        source of the System_azo object. 

        To ensure the geometry could be treated a posteriori to create cis or transition states geometries, the SMILES string
        should follow the template: c1(...1)/N=N/ring2. 
            
            e.g. c1(ccccc1)/N=N/c2ccccc2

        If no System_azo is provided, it returns labels, coord and the Specie_azo named 'trans'.  

        Workflow
        --------
            - Find or create trans isomer.
            - Get indices of atoms relevant to dihedral adjustment.
            - Move the azo group dihedral angle close to the target angle.
            - If there is steric hindrance, tries to solve it by moving adjacent rings.
            - If steric hindrance cannot be solved, returns None.
            - If it is solved, structure is saved to a Specie_azo object with an 'initial' State. 

        Parameters
        ----------
            self: System_azo object
                System_azo object where the trans isomer will be stored.
            smiles: str
                SMILES string of the azo compound. It should follow the template: c1(...1)/N=N/ring2.
            debug: int
                Debug level. 0: no debug, 1: verbose debug

        Returns
        -------
            labels: list
                List of atom labels.
            coord: list
                List of atomic coordinates.
            iso: Specie_azo
                Specie_azo object containing the trans isomer.
        '''

        if debug > 0: print(f"AZO.CREATE_TRANS: Creating trans isomer for {self.name}")

        if '/N=N\\' in smiles:
            if debug > 0: print(f"AZO.CREATE_TRANS: Received SMILES of CIS isomer, replacing '\\' with '/' for {self.name}")
            smiles = smiles.replace('/N=N\\', '/N=N/')

        labels,coord = get_3D_openbabel(smiles)
        coord        = centercoords(coord, 0)
        trans        = Specie_azo(labels, coord)
        trans.smiles = smiles
        if trans.check_fragmentation():
            print(f"AZO.CREATE_TRANS: Trans isomer for {self.name} is FRAGMENTED.")
            return None

        ## Aixo s'ha de revisar. La carrega s'ha de treure de l'smiles
        trans.set_total_charge(0)
        trans.set_total_spin(0)            

        trans_state = trans.add_state("initial")
        trans_state.set_geometry(labels, coord)
        self.add_source(name="trans", new_source = trans)
        return trans

    ######
    def create_cis(self, target_deg: float=40.0, max_iter: int=2000, debug: int=0) -> "Specie_azo":
        '''
        Creates the cis structure of the azo compound from a SMILES string. To avoid troubles in 3D geometry creation using openbabel, 
        the trans isomer is created using create_trans() function. It sets up the cis isomer (including the creation of Specie_azo and 
        its 'initial' State objects) and stores it as source of the System_azo object. 

        If no System_azo is provided, it returns labels, coord and the Specie_azo named 'cis'.  

        Workflow
        --------
        - Find or create trans isomer.
        - Get indices of atoms relevant to dihedral adjustment.
        - Move the azo group dihedral angle close to the target angle.
        - If there is steric hindrance, tries to solve it by moving adjacent rings.
        - If steric hindrance cannot be solved, returns None.
        - If it is solved, structure is saved to a Specie_azo object with an 'initial' State. 

        Parameters
        ----------
        target_deg: float
            Target angle of the azo dihedral in degrees. 
        max_iter: int
            Maximum number of iterations to solve steric hindrance.
        debug: int
            Debug level. 0: no debug, 1: verbose debug

        Returns
        -------
        labels: list
            List of atom labels.
        coord: list
            List of atomic coordinates.
        iso: Specie_azo
            Specie_azo object containing the cis isomer.
        '''

        # 1st-searching for trans isomer
        trans_found, trans = self.find_source('trans')
        if not trans_found: raise Exception(f"AZO.CREATE_CIS: [ERROR] Trans isomer not found. Create it first with self.create_trans()")
        
        # Get trans isomer geometry as reference
        labels = trans.labels
        coord  = trans.coord

        # Get indices for the azo group (at1 - at2 = at3 - at4) and neighbours (at0, at5)
        at0, at1, at2, at3, at4, at5 = self.get_dihedral_indices()
        if debug > 0: print(f"AZO.CREATE_CIS: dihedral indices: {at0}, {at1}, {at2}, {at3}, {at4}, {at5}")

        # Get the adjacency matrix for reference when rotating the dihedral angle of the azo group
        adjmat_ref, adjnum_ref = trans.get_adjmatrix()

        # Initial dihedral angle, for info 
        current_rad = get_dihedral(coord[at1],coord[at2],coord[at3],coord[at4]) # Initial dihedral angle
        current_deg = np.degrees(current_rad)
        if debug > 0: print(f"AZO.CREATE_CIS: Initial dihedral angle: {current_deg} degrees")

        # Pre-Rotation. Getting dihedral close to target
        if abs(current_deg - target_deg) >= 5:
            if debug > 0: print(f"AZO.CREATE_CIS: current dihedral {abs(current_deg)} is further than {target_deg+5} degrees from target. Adjusting...")
            jump_target = target_deg+1 if current_deg > 0 else -target_deg-1 
            if debug > 0: print(f"AZO.CREATE_CIS: Jumping to {jump_target} degrees")
            if debug > 0: print(f"AZO.CREATE_CIS: Selected atoms: {at1}, {at2}, {at3}, {at4}")
            coord_next = set_dihedral(labels, coord, jump_target, at1,at2,at3,at4,adjmat=adjmat_ref, adjnum=adjnum_ref)
            angle_next = np.degrees(get_dihedral(coord_next[at1],coord_next[at2],coord_next[at3],coord_next[at4]))
            if debug > 0:  print(f"AZO.CREATE_CIS: Changed dihedral in {self.name} from {current_deg} to {angle_next}, it should be near +-{target_deg}")
            _, adjmat_cis, adjnum_cis = get_adjmatrix(labels,coord_next)
        else:
            _, adjmat_cis, adjnum_cis = get_adjmatrix(labels,coord)


        found_geometry = False 

        matrices_match = np.array_equal(adjmat_cis, adjmat_ref) and np.array_equal(adjnum_cis, adjnum_ref)
 
        if debug > 0: print(f"AZO.CREATE_CIS: Matrices match: {matrices_match}")
        if matrices_match:
            coord = coord_next
            current_deg = angle_next
            if debug > 0: print(f"AZO.CREATE_CIS: Angle {current_deg:.2f} (OK)")
            if abs(current_deg) <= (target_deg + 5):
                found_geometry = True
        else:
            fixed_collision, coord = solve_dihedral(labels, coord_next, at0, at1, at2, at3, at4, at5, adjmat_ref=adjmat_ref, adjnum_ref=adjnum_ref,debug=debug)
            current_deg = np.degrees(get_dihedral(coord[at1],coord[at2],coord[at3],coord[at4]))
            if debug > 0: print(f"AZO.CREATE_CIS: Fixed collision: {fixed_collision}")
            if abs(current_deg) <= (target_deg + 5) and fixed_collision:
                if debug > 0: print(f"AZO.CREATE_CIS: Angle {current_deg:.2f} (OK)")
                found_geometry = True
                if debug > 0: print(f'AZO.CREATE_CIS: Found good geometry for {self.name} by rotating adjacent dihedrals')
            else: 
                if debug > 0: print(f'AZO.CREATE_CIS: Failed to find good geometry for {self.name} by rotating adjacent dihedrals')
        
        if found_geometry: 
            coord      = centercoords(coord, at1)
            cis        = Specie_azo(labels, coord) 
            cis.smiles = trans.smiles.replace('/N=N/','/N=N\\')

            # Aixo s'ha de revisar
            cis.set_total_charge(0)
            cis.set_total_spin(0)

            self.add_source(name="cis", new_source = cis)

            if cis.check_fragmentation():
                print(f"AZO.CREATE_CIS: Cis isomer for {self.name} is FRAGMENTED.")
                return None, None, None
            
            ini_state = cis.add_state("initial")
            ini_state.set_geometry(labels, coord)
            return cis
        else:
            raise Exception(f'AZO.CREATE_CIS: Target dihedral for {self.name} could not be reached. Reached max. iterations: {max_iter}.')

    ######
    def create_ts(self, ts_list:list = ['TSrot', 'TSinv_l', 'TSinv_r', 'triplet'], debug: int=0):
        """
        Creates a set of TS for a given System_azo. Users can select which TS to create from a list of options, 
        it must be in ts_list (['TSrot', 'TSinv', 'triplet']). By default, TSrot_A, TSrot_B, TSinv_L, TSinv_R are created, 
        including triplet states TSrot_A_T, TSrot_B_T. 
        
        If 'triplet' is selected on ts_list, rotation TS in triplet state are created.

        Once the TS are created, they are added to the System_azo object.

        TS labels and coordinates are returned as a dictionary, containing the labels (created_ts['key'][0]) and coordinates (created_ts['key'][1]) of each TS. 
        They can be accessed using the keys 'TSrot_A', 'TSrot_B', 'TSinv_L', 'TSinv_R', 'TSrot_A_T', 'TSrot_B_T'.

        Important notes
        ---------------

        -------
        TSrot
        -------
        TSrot is created by rotating the trans isomer by +/-90 degrees around the azo dihedral. 
        If 'triplet' is selected on ts_list, rotation TS in triplet state are created.

        Nomenclature: 
            - TSrot_A_S: Rotation TS in singlet state. Dihedral angle is set as +90º. Total spin is set as 1.
            - TSrot_A_T: Rotation TS in triplet state. Dihedral angle is set as +90º. Total spin is set as 3. MULTIPLICITY
            - TSrot_B_S: Rotation TS in singlet state. Dihedral angle is set as -90º. Total spin is set as 1.
            - TSrot_B_T: Rotation TS in triplet state. Dihedral angle is set as -90º. Total spin is set as 3.

        -------
        TSinv
        -------
        TSinv is created by inverting the trans isomer around the azo dihedral (setting dihedral angle to 180º).
        By default, total spin is set as 1. It can be changed using the function for Species as Specie.set_total_spin(value).
        Two versions are created:
            - TSinv_L: Inversion TS involving inversion of left ring (at0 - at1 = at2). Total spin is set as 1.
            - TSinv_R: Inversion TS involving inversion of right ring (at1 = at2 - at3). Total spin is set as 1.

        They are added as sources of the azosystem. They can be accessed by using azosystem.find_source('TSrot_A_T') or azosystem.find_source('TSinv_R')
        
        """

        # 1st-searching for trans isomer
        trans_found, trans = self.find_source('trans')
        if not trans_found: raise Exception(f"AZO.CREATE_TS: [ERROR] Trans isomer not found. Create it first with self.create_trans()")

        # Get trans isomer geometry as reference
        labels = trans.labels
        coord = trans.coord

        # Get indices for the azo group (at1 - at2 = at3 - at4) and neighbours (at0, at5)
        at0, at1, at2, at3, at4, at5 = self.get_dihedral_indices()
        if debug > 0: print(f"AZO.CREATE_CIS: dihedral indices: {at0}, {at1}, {at2}, {at3}, {at4}, {at5}")

        # Get the adjacency matrix for reference when rotating the dihedral angle of the azo group
        adjmat_ref, adjnum_ref = trans.get_adjmatrix()

        if 'TSrot' in ts_list:
            ## TSrot_A ##
            coord = set_dihedral(labels, trans.coord, 90, at1,at2,at3,at4, adjmat=adjmat_ref, adjnum=adjnum_ref)
            coord = set_dihedral(labels, coord, 0, at2,at3,at4,at5, adjmat=adjmat_ref, adjnum=adjnum_ref)
            _, adjmat, adjnum = get_adjmatrix(labels,coord)
            is_equal = np.array_equal(adjmat, adjmat_ref) and np.array_equal(adjnum, adjnum_ref)
            if is_equal:
                found = True
                coord = centercoords(coord, at0)
            else:
                found, coord = solve_dihedral(labels, coord, at0, at1, at2, at3, at4, at5, adjmat_ref=adjmat_ref, adjnum_ref=adjnum_ref, debug=debug)
            
            if found:
                ts = Specie_azo(labels, coord)
                isFragmented = ts.check_fragmentation()  # Check if the TSrot is fragmented
                if not isFragmented:
                    state = ts.add_state("initial")
                    state.set_geometry(labels, coord)
                    # Aixo s'ha de revisar
                    ts.set_total_charge(0)
                    ts.set_total_spin(0)
                    self.add_source('TSrot_A_S', ts)

                    if 'triplet' in ts_list:
                        ts_triplet = Specie_azo(labels, coord)
                        # Aixo s'ha de revisar
                        ts_triplet.set_total_charge(0)
                        ts_triplet.set_total_spin(2)
                        triplet_state = ts_triplet.add_state("initial")
                        triplet_state.set_geometry(labels, coord)
                        self.add_source('TSrot_A_T', ts_triplet)
                        if debug > 0: print(f'AZOS.CREATE_TS.TSROT_A: TSrot_A_T Specie_Azo successfully created for {self.name}')
                    if debug > 0: print(f'AZOS.CREATE_TS.TSROT_A: TSrot_A_S Specie_Azo successfully created for {self.name}')
                else:
                    raise Exception(f'AZOS.CREATE_TS.TSROT_A: [ERROR] TSrot_A fragmented for {self.name}')

            ## TSrot_B ##                
            if debug > 0:
                dg_deg = np.degrees(get_dihedral(trans.coord[at1], trans.coord[at2], trans.coord[at3], trans.coord[at4]))
                print(f'AZOS.CREATE_TS.TSROT_B: Dihedral angle for reference geometry: {dg_deg} degrees')
            coord = set_dihedral(labels, trans.coord, -90, at1,at2,at3,at4, adjmat=adjmat_ref, adjnum=adjnum_ref)            # Coords de tsrot
            _, adjmat, adjnum = get_adjmatrix(labels,coord)
            is_equal = np.array_equal(adjmat, adjmat_ref) and np.array_equal(adjnum, adjnum_ref)
            if is_equal: found = True
            else:        found, coord = solve_dihedral(labels, coord, at0, at1, at2, at3, at4, at5, adjmat_ref=adjmat_ref, adjnum_ref=adjnum_ref, debug=debug)
    
            if found:
                coord = centercoords(coord, at0)
                ts = Specie_azo(labels, coord)
                isFragmented = ts.check_fragmentation()  # Check if the TSrot is fragmented
                if not isFragmented:
                    # Aixo s'ha de revisar
                    ts.set_total_charge(0)
                    ts.set_total_spin(0)
                    state = ts.add_state("initial")
                    state.set_geometry(labels, coord)
                    self.add_source('TSrot_B_S', ts) # tsrot created from E-isomer 
                    if 'triplet' in ts_list:
                        ts_triplet = Specie_azo(labels, coord, name="TSrot_B_T")
                        # Aixo s'ha de revisar
                        ts_triplet.set_total_charge(0)
                        ts_triplet.set_total_spin(2)
                        triplet_state = ts_triplet.add_state("initial")
                        triplet_state.set_geometry(labels, coord)
                        self.add_source('TSrot_B_T', ts_triplet)
                        if debug > 0: print(f'AZOS.CREATE_TS.TSROT_B: TSrot_B_T Specie_Azo successfully created for {self.name}')
                    if debug > 0: print(f'AZOS.CREATE_TS.TSROT_B: TSrot_B_S Specie_Azo successfully created for {self.name}')
                else:
                    raise Exception(f'WARNING: TSrot_B fragmented for {self.name}')
        
        ## TSinv Left ##
        if 'TSinv_l' in ts_list:
            coord = set_angle(labels, trans.coord, 179.9, at1,at2,at3)
            if debug > 0: 
                angle_deg = np.degrees(get_angle(coord[at1]-coord[at2], coord[at3]-coord[at2]))
                print(f'AZOS.CREATE_TS.TSINV_L: Angle between {at1}, {at2} and {at3} set to {angle_deg} degrees.')
            angles = np.concatenate(([0],[val for i in range(1, 12) for val in (15 * i, -15 * i)]))
            for a0 in angles:
                coord = set_dihedral(labels, coord, a0, at3,at2,at1,at0, adjmat=adjmat_ref, adjnum=adjnum_ref)
                _, adjmat, adjnum = get_adjmatrix(labels,coord)
                is_equal = np.array_equal(adjmat, adjmat_ref) and np.array_equal(adjnum, adjnum_ref)
                if is_equal:
                    ts = Specie_azo(labels, coord)
                    ts_isFragmented = ts.check_fragmentation()  # Check if the TSinv_l is fragmented
                    if not ts_isFragmented:
                        # Aixo s'ha de revisar
                        ts.set_total_charge(0)
                        ts.set_total_spin(0)
                        ts_state = ts.add_state("initial")
                        ts_state.set_geometry(labels, coord)
                        self.add_source('TSinv_l', ts)
                        if debug > 0: print(f'AZOS.CREATE_TS.TSINV_L: TSinv_l Specie_Azo successfully created for {self.name}')
                        break
                    else:
                        raise Exception(f'AZOS.CREATE_TS.TSINV_L: [ERROR] TSinv_l fragmented for {self.name}')
        
        ## TSinv Right ##
        if 'TSinv_r' in ts_list:
            coord = set_angle(labels, coord, 179.9, at2,at3,at4) ## Angle value of 179.9 instead of 180 to avoid numerical issues
            if debug > 0: 
                angle_deg = np.degrees(get_angle(coord[at3]-coord[at2], coord[at4]-coord[at2]))
                print(f'AZOS.CREATE_TS.TSINV_R: Angle between {at2}, {at3} and {at4} set to {angle_deg} degrees.')
            angles = np.concatenate(([0],[val for i in range(1, 12) for val in (15 * i, -15 * i)]))
            for a0 in angles:
                coord = set_dihedral(labels, coord, a0, at2,at3,at4,at5, adjmat=adjmat_ref, adjnum=adjnum_ref)
                _, adjmat, adjnum = get_adjmatrix(labels,coord)
                is_equal = np.array_equal(adjmat, adjmat_ref) and np.array_equal(adjnum, adjnum_ref)
                if is_equal:
                    ts = Specie_azo(labels, coord)
                    ts_isFragmented = ts.check_fragmentation()  # Check if the TSinv_l is fragmented
                    if not ts_isFragmented:
                        # Aixo s'ha de revisar
                        ts.set_total_charge(0)
                        ts.set_total_spin(0)
                        ts_state = ts.add_state("initial")
                        ts_state.set_geometry(labels, coord)
                        self.add_source('TSinv_r', ts)
                        if debug > 0: print(f'AZOS.CREATE_TS.TSINV_R: TSinv_r Specie_Azo successfully created for {self.name}')
                        break
                    else:
                        raise Exception(f'AZOS.CREATE_TS.TSINV_R: [ERROR] TSinv_r fragmented for {self.name}')
        return 

    def get_PSS(self, lamp : "Lamp", phi_EZ = 0.3, phi_ZE = 0.5, t_EZ=None, t_ZE=None, debug=0):

        lambda_grid, sigma_Z, sigma_E = self.get_abs_spectrum(normalize=False, units=True, get_PSS=True)
        Z_exists, cis = self.find_conformer('cis')
        E_exists, trans = self.find_conformer('Trans')
        if not Z_exists or not E_exists:
            print(f'Z_exists: {Z_exists}, E_exists: {E_exists}')
            raise ValueError(f'Error: No cis or trans conformers found for system {self.name}')
        cis_state = cis.find_state('opt')[1]
        trans_state = trans.find_state('opt')[1]
        # Photon flux from lamp
        photon_flux = get_photon_flux_spectrum(lamp.wavelength, lamp.fwhm, lambda_grid, Itot=lamp.irradiance)
        # Half-lives in seconds
        if not hasattr(trans_state, 'halflife'):
            print('No halflife found for trans isomer. Computing...')
            self.set_iso_halftime('trans', skip_triplets=True, overwrite=False)
        if not hasattr(cis_state, 'halflife'):
            print('No halflife found for cis isomer. Computing...')
            self.set_iso_halftime('cis', skip_triplets=True, overwrite=False)
            raise ValueError(f'Error: No halflife found for Z or E isomers in system {self.name}')
        if t_EZ is None:
            t_EZ = trans_state.results['halflife'].value
        if t_ZE is None:
            t_ZE = cis_state.results['halflife'].value
        # Photochemical rates
        k_ph_EZ = phi_EZ * np.trapz(sigma_E * photon_flux, lambda_grid)
        k_ph_ZE = phi_ZE * np.trapz(sigma_Z * photon_flux, lambda_grid)
        # Thermal rates
        k_th_EZ = np.log(2) /  t_EZ
        k_th_ZE = np.log(2) /  t_ZE
        # Steady-state populations
        N_E = (k_th_ZE + k_ph_ZE) / (k_ph_EZ + k_ph_ZE + k_th_EZ + k_th_ZE)
        return lambda_grid, sigma_Z, sigma_E, N_E
    
    def get_abs_spectrum(self, normalize: bool = False, units: bool = False, get_PSS: bool = False, custom_cis: str = None, custom_trans: str = None):
        if custom_cis is not None:
            name = str(custom_cis)
            Z_exists, cis = self.find_conformer(custom_cis)
        else:
            Z_exists, cis = self.find_conformer('cis')
        if custom_trans is not None:
            name = str(custom_trans)
            E_exists, trans = self.find_conformer(custom_trans)
        else:
            E_exists, trans = self.find_conformer('Trans')
        if not Z_exists or not E_exists:
            print(f'Z_exists: {Z_exists}, E_exists: {E_exists}')
            raise ValueError(f'Error: No cis or trans conformers found for system {self.name}')
        opt_Z_exists, cis_state = find_state(cis, 'opt')
        opt_E_exists, trans_state = find_state(trans, 'opt')
        if not opt_Z_exists or not opt_E_exists:
            print(f'opt_Z_exists: {opt_Z_exists}, opt_E_exists: {opt_E_exists}')
            raise ValueError(f'Error: No opt state found for Z or E isomers')

        # Add checks for thermal stability
        if not hasattr(cis_state, 'es_list') or not hasattr(trans_state, 'es_list'):
            print('WARNING: No TDDFT data found for cis or trans isomers')
            return None, None, None, None

        Z_e = [es.energy for es in cis_state.es_list]
        Z_f = [es.fosc for es in cis_state.es_list]
        E_e = [es.energy for es in trans_state.es_list]
        E_f = [es.fosc for es in trans_state.es_list]
        Emin = min(min(Z_e), min(E_e)) - 1
        Emax = max(max(Z_e), max(E_e)) + 1
        x = np.linspace(Emin, Emax, 5000)
        sigma_Z = build_sigma(zip(Z_e, Z_f), x, normalize=False,units=False)     # Absolute
        sigma_E = build_sigma(zip(E_e, E_f), x, normalize=False,units=False)
        if get_PSS:
            sigma_Z = build_sigma(zip(Z_e, Z_f), x, normalize=False,units=True)     # Absolute
            sigma_E = build_sigma(zip(E_e, E_f), x, normalize=False,units=True)
        return 1240 / x[::-1], sigma_Z, sigma_E  

    def __repr__(self):
        to_print = ""
        to_print += f'------------- SCOPE Azo System --------------\n'
        to_print += f' Name:              {self.name}\n'
        if hasattr(self,"dE"):      to_print += f' Thermal Stability = {self.dE} kJ/mol (- means trans is more stable)'
        to_print += '---------------------------------------------\n'
        to_print += '                                              \n'
        return to_print

###################################
##### SPECIE Object Adapted to Specie_azo Class #####
###################################

class Specie_azo(Specie):

    def __init__(self, labels, coord, name: str):
        Specie.__init__(self, labels, coord)
        self.subtype  = "specie_azo"

    def correct_tripletG(self, triplet_specie, T:float=298.15, overwrite = False, p_sh:float = 0.0002, debug: int=0):
        '''
        Corrects the Gtot of a triplet Specie_azo object using the Gtot of the parent specie. 
        Correction is done considering the increase of energy due to surface hopping between the singlet and triplet PESs. 

        Parameters
        ----------
        triplet_specie : Specie_azo
            The triplet specie object to correct the Gtot of.
        T : float, optional
            The temperature in Kelvin. The default is 298.15 K.
        overwrite : bool, optional
            Whether to overwrite the existing Gtot_corr value. The default is False.
        p_sh : float, optional
            The probability of surface hopping. The default is 0.0002.
        debug : int, optional
            The debug level. The default is 0
        '''
        k_b = Constants.boltz_J # J/K
        h = Constants.planck_Js # J·s
        R = Constants.R_J       # 8.31 J/(K·mol)

        found_iso_opt, iso_opt = self.find_state("opt")
        found_triplet_opt, triplet_opt = triplet_specie.find_state("opt")

        if not found_iso_opt or not found_triplet_opt:
            print(f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: No opt state found for {self.name} or {triplet_specie.name}.')
            return

        parent = self._sys
        exist = 'Gtot_corr' in triplet_opt.results

        if not 'Gtot' in iso_opt.results or not 'Gtot' in triplet_opt.results:
            print(f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: No Gtot found for {self.name} or {triplet_specie.name}.')
            return
        
        if not exist or overwrite:
            if debug > 0: print(f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: Found Gtot for {self.name} and {triplet_specie.name}. Correcting Gtot of {parent.name} triplet state.')
            G_triplet = triplet_opt.results['Gtot'].value
            G_iso = iso_opt.results['Gtot'].value

            if debug > 0: print(f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: {parent.name} triplet G: {G_triplet} hartree, iso G: {G_iso} hartree')

            dG = (G_triplet - G_iso) * Constants.har2kJmol * 1000  # in J/mol
            t, k = compute_t(G_triplet, G_iso, T)
            k_sh = k * p_sh *p_sh 
            deltax = - (dG + R * T * np.log((h*k_sh)/(k_b*T))) # in J/mol

            if debug > 0: print(f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: Adding deltax in kcal/mol: {deltax*0.24/1000}')
            newG = (G_triplet + deltax / (1000*Constants.har2kJmol))
            newG = float(newG)
            newdata = Data("Gtot_corr", newG, "au", "correct_triplet_G")
            
            triplet_opt.add_result(newdata)
            if debug > 0: print(f'AZO.SPECIE_AZO.CORRECT_TRIPLETG: Corrected Gtot of {parent.name} triplet state by {deltax*0.24/1000:.2f} Kcal/mol.')
            return newG

    def set_iso_halftime(self, skip_triplets : bool = True, overwrite = False):
        '''
        Computes t0.5 in seconds for a given conformer/isomer e.g. cis or trans using the Eyring equation.
        Saves the result as a data object with the key 'halftime' in the conformer object.
        Minimum energy Transition State can be accessed using the key 'mets' in the conformer object. E. g. cis.mets

        Parameters
        ----------
        self : Specie_azo
            The specie_azo object to compute the halftime for.
        skip_triplets : bool
            Skip triplet conformers in halftime calculation.
        overwrite : bool
            Overwrite existing t0.5 values.

        Notes
        -----
        The function will only consider Specie_azo objects that have an opt state with a Gtot value.
        The function will only consider TSs that have an opt state with a Gtot value, or Gtot_corr value if skip_triplets is False.
        
        '''
        ts_values = []
        ts_names = []

        found_iso_opt, iso_state = self.find_state("opt")

        if not found_iso_opt:
            raise Exception(f'AZO.SPECIE_AZO.SET_HALFTIMEOptimization state not found for {self.name}.')
             
        if 'Gtot' in iso_state.results.keys(): g_iso = iso_state.results['Gtot'].value
        else: raise ValueError('AZO.SPECIE_AZO.SET_HALFTIME: Gtot not found for isomer')

        parent = self._sys
        candidates = [source for source in parent.sources if not source.name.lower().startswith('cis') or not source.name.lower().startswith('trans')]
        candidates_names = [source.name for source in candidates]

        if 'TSrot_A_T' in candidates_names or 'TSrot_B_T' in candidates_names:
            triplets = [source for source in parent.sources if source.name.lower().startswith('TSrot_A_T') or source.name.lower().startswith('TSrot_B_T')]
            for triplet in triplets:
                correct_tripletG(self, triplet)

        for ts in candidates:
            found_ts_state, ts_state = ts.find_state("opt")
            if not found_ts_state or not 'Gtot' in ts_state.results.keys():
                print(f'AZO.SPECIE_AZO.SET_HALFTIME: [WARNING] Optimization state or Gtot not found for {ts.name}.')
                continue
            # Use only TSs with Gtot or Gtot_corr

            if ts.spin == 3:     # Use corrected Gtot for triplets
                if not skip_triplets:
                    if 'Gtot_corr' in ts_state.results.keys(): energy = ts_state.results['Gtot_corr'].value
                    else: raise ValueError(f'AZO.SPECIE_AZO.SET_HALFTIME: Corrected Gtot for {ts.name} Specie_azo not found for Triplet TS, altough it was corrected with correct_tripletG() function.')
                else:
                    continue
            else:   energy = ts_state.results['Gtot'].value

            name = ts.name 
            # if type(energy) is not float:
            #     energy = 0.
            #     print("Warning: energy is not a float for", name, 'Setting to 0.')
            ts_names.append(name)
            ts_values.append(energy)
        
        if not hasattr(self, 'halflife') or overwrite:
            if debug > 0: print(rf'AZO.SPECIE_AZO.SET_HALFTIME: Collected {len(ts_values)} TSs for {self.name} : {ts_names} with energies {ts_values}.')
            if debug > 0: print(f'AZO.SPECIE_AZO.SET_HALFTIME:Doing halftime for {parent.name} {self.name}')
            # Choosing Minimum Energy TS (mets)
            min_idx = int(np.argmin(ts_values))
            mets = ts_names[min_idx]                        
            g_cross = ts_values[min_idx]
            dG_cross = (float(g_cross)- float(g_iso)) * Constants.har2kJmol * 0.24 # in Kcal/mol
            # Compute and store halftime
            t,k = compute_t(float(g_cross), float(g_iso))           
            if debug > 0: print(t, ' s for isomer ', self.name)
            new_time = Data('halflife', float(t), 's', 'compute_t')
            self.halflife = float(t)
            self.mets = mets
            self.add_result(new_time, overwrite)
            
            if debug > 0: print(dG_cross, 'for isomer ', self.name)
            self.dG_cross = dG_cross
            newdata = Data('dG_cross', float(dG_cross), 'kcal/mol', 'set_iso_halftime')
            self.add_result(newdata, overwrite)
        else: 
            print(f'State not found in {self.name}.')
        print('Done! Note that missing energy values have been set to 0.')

    def link_tda_to_state(self, state: object, filepath: str, overwrite: bool=False):
        '''
        Parses excited states from a TDDFT Gaussian16 .log file to a ExcitedState object.
        '''
        if not hasattr(state,"es_list") or overwrite: 
            if os.path.exists(filepath):
                state.tda_filepath = filepath
                state.es_list = []
                lines = read_lines_file(filepath)
                state.gs_energy = parse_energy_from_step(lines)
                if parse_status_finished(lines):
                    for st_num in range(1,11):
                        if st_num < 10:    line_nums, found = search_string(f"Excited State   {st_num}:", lines, typ='first')
                        elif st_num == 10: line_nums, found = search_string(f"Excited State  {st_num}:", lines, typ='first')
                        if found:
                            dummy, dummy, idx, dummy, energy, dummy, wavelength, dummy, fosc, s2 = lines[line_nums].split() 
                            idx  = int(idx.replace(':',''))
                            s2   = float(s2.replace('<S**2>=',''))
                            fosc = float(fosc.replace('f=',''))
                            wavelength = float(wavelength)
                            energy = float(energy)
                            if st_num == idx:
                                new_es = ExcitedState(self, st_num, energy, wavelength, fosc, s2)
                                state.es_list.append(new_es)
                        else: print(f"Excited State {st_num} not found in {filepath}.") 
            else: print(f"File {filepath} does not exists.") 

    def link_opt_to_state(self, state:object, filepath:str, overwrite: bool=False):
        '''
        Parses a Gaussian16 .log file to save geometry and free Gibbs energy to a state object.
        '''
        if not hasattr(state, 'opt') or overwrite:
            if os.path.exists(filepath):

                state.opt_filepath = filepath
                lines = read_lines_file(filepath)
                opt_finished = parse_opt_status(lines)
                if not opt_finished: raise ValueError('Optimization not finished')
                state.energy = parse_energy(lines)
                state.gtot = parse_free_energy(lines)

                labels, coord = parse_last_geometry(lines, debug=1)
                state.set_geometry(labels, coord)
                newG = data("Gtot", state.gtot, "au", "parse_free_energy")
                state.add_result(newG)
        else:
            print(f"File {filepath} does not exist.")

    def get_abs_spectrum(self, normalize: bool = False, units: bool = False, custom_cis: str = None, custom_trans: str = None):
        if custom_cis is not None:
            name = str(custom_cis)
            Z_exists, cis = self.find_conformer(custom_cis)
        if custom_trans:
            name = str(custom_trans)
            E_exists, trans = self.find_conformer(custom_trans)
        if not Z_exists or not E_exists:
            print(f'Z_exists: {Z_exists}, E_exists: {E_exists}')
            raise ValueError(f'Error: No cis or trans conformers found for system {self.name}')
        opt_Z_exists, cis_state = find_state(cis, 'opt')
        opt_E_exists, trans_state = find_state(trans, 'opt')
        if not opt_Z_exists or not opt_E_exists:
            print(f'opt_Z_exists: {opt_Z_exists}, opt_E_exists: {opt_E_exists}')
            raise ValueError(f'Error: No opt state found for Z or E isomers')

        # Add checks for thermal stability
        if not hasattr(cis_state, 'es_list') or not hasattr(trans_state, 'es_list'):
            print('WARNING: No TDDFT data found for cis or trans isomers')
            return None, None, None, None

        Z_e = [es.energy for es in cis_state.es_list]
        Z_f = [es.fosc for es in cis_state.es_list]
        E_e = [es.energy for es in trans_state.es_list]
        E_f = [es.fosc for es in trans_state.es_list]
        Emin = min(min(Z_e), min(E_e)) - 1
        Emax = max(max(Z_e), max(E_e)) + 1
        x = np.linspace(Emin, Emax, 5000)
        sigma_Z = build_sigma(zip(Z_e, Z_f), x, normalize=False,units=False)     # Absolute
        sigma_E = build_sigma(zip(E_e, E_f), x, normalize=False,units=False)
        return 1240 / x[::-1], sigma_Z, sigma_E 

    def __repr__(self):
        to_print = ""
        to_print += f'------ Specie_azo custom SPECIE object ------\n'
        to_print += Specie.__repr__(self, indirect=True)
        to_print += '-----------------------------------------------\n'
        return to_print

#######################################
#####     General Lamp Class      #####
#######################################
class Lamp:
    def __init__(self, name : str, wavelength : float, fwhm_nm : float = None, power : float = None, shift_nm : float = 0):
        '''
        Parameters
        ----------
            name: string 
                Name of the lamp.
            wavelength: float
                Irradiation wavelength in nm. 
            fwhm_nm: float
                Full Width at Half Maximum in nm (if None, FWHM=1nm).
            power: float, 
                Power of the lamp in W. If None, power is taken from literature data.

        
        '''
        power_data = {      # Power data (W) from https://doi.org/10.1016/j.bcp.2025.117065
            365: 0.10,
            385: 0.16,
            405: 0.10,
            435: 0.025,
            460: 0.125,
            470: 0.35,
            500: 0.025
            }        
        irr_data = {        # Irradiance at 50% intensity (mW/mm2)
            365: 1.04,
            385: 2.6,
            405: 2.10,
            435: 0.72,
            460: 2.17,
            470: 1.02,
            490: 0.95,
            500: 0.3,
            525: 0.36,
            550: 1.57
            }
        
        fwhm_list = {            # Fwhm data (nm) from COOLLED light source
            365: 12.54,
            385: 11.38,
            405: 16.08,
            435: 13.02,
            460: 17.61,
            470: 23.26,
            500: 24.72,
            525: 28.32,
            550: 82.79
        }
        
        self.name = name
        self.wavelength = float(wavelength) # in nm

        if fwhm_nm is None:
            if wavelength not in fwhm_list.keys():
                raise ValueError(f'Wavelength {wavelength} nm not available for lamp {name}. Choose from {list(fwhm_list.keys())} nm.')
            else:
                self.fwhm = fwhm_list.get(int(round(wavelength)))
        else:
            self.fwhm = float(fwhm_nm)

        if wavelength not in irr_data.keys():   raise ValueError(f'Wavelength {wavelength} nm not available for lamp {name}. Choose from {list(irr_data.keys())} nm.')
        else:   self.irradiance = irr_data.get(int(wavelength)) # in mW/mm2   

        if shift_nm != float(0.0):  self.eff_wavelength = self.wavelength - shift_nm
        else:   self.eff_wavelength = self.wavelength

    def __repr__(self):
        print('--------------LAMP OBJECT--------------')
        print(f'Name:           {self.name}')
        print(f'Wavelength:                 {self.wavelength} nm')
        if hasattr(self, "eff_wavelength"): 
            print(f'Wavelength (after shift):   {self.eff_wavelength} nm')
            
        print(f'FWHM:                       {self.fwhm} nm')
        if hasattr(self, "power"): 
            print(f'Power:                      {self.power} W')
        print('---------------------------------------')
        print('NOTE: If no shift is applied, both wavelengths are equal.')
        return ''
