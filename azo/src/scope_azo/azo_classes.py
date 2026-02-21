import numpy as np
from    scope                               import *
from    scope.classes_system                import *
from    scope.classes_specie                import *
from    scope.software.gaussian.g16_parse   import *
from    scope.geometry                      import *
from    scope.connectivity                  import *
from    scope_azo.azo_functions             import *
from    scope.parse_general                 import search_string, read_lines_file
from    scope.software.gaussian.g16_parse   import * 
from    scope.software.gaussian.g16_output  import * 
from    scope.classes_data                  import *
from    scope.classes_qc                    import *

############################################
#### SYSTEM Object adapted to System_azo ###
############################################

class System_azo(System):
    def __init__(self, name: str, smiles: str, debug: int=0):
        System.__init__(self, name)
        self.subtype          = "system_azo"        
        self.name             = name
        self.smiles           = smiles
        self.dihedral_indices = self.get_dihedral_indices()

    ######
    def get_dihedral_indices(self):
        '''
        Extracts the indices of the atoms involved in the dihedral angle describing the isomerization in azo species. 
        It works for any azo SMILES with the following structure: ring1/N=N/ring2. 
        It is based on the assumption that when the 3D geometry is created from the smiles, the H atoms are added at the end.

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
        - at0 would be the atom with index 0, found in the left ring 
        - at1 would be the atom with index 1, found in the left ring 
        - at2 would be the Nitrogen atom with index 6, found in the Azo fragment
        - at3 would be the Nitrogen atom with index 7, found in the Azo fragment
        - at4 would be the atom with index 8, found in the right ring 
        - at5 would be the atom with index 9, found in of the right ring

        Atoms at1, at2, at3 and at4 can be used to define the dihedral angle of the azo fragment.
        Atoms at0 and at5 can be used to define the dihedral angle of the adjacent rings with respect to the azo fragment.
        '''

        # 1) First, we identify the two rings from the smiles string
        chars           = ElementData().elementname.keys() # List of characters that define an atom. Update if needed  
        chars_lowercase = [char.lower() for char in chars]
        at_count        = 1
        fragments       = []
        block           = []
        
        # New flag to track if we are inside a bracketed atom like [NH3+]
        inside_bracket = False 

        for char in self.smiles:
            if char == '/' or char == '\\': 
                fragments.append(block)
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
        fragments.append(block)

        #2) Second, we get the actual atom indices
        at0 = fragments[0][1] -1 
        at1 = fragments[0][0] -1 
        at2 = fragments[1][0] -1 
        at3 = fragments[1][1] -1 
        at4 = fragments[2][0] -1 
        at5 = fragments[2][1] -1 
        return list([at0, at1, at2, at3, at4, at5])

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
    def create_trans(self, overwrite: bool=False, debug: int = 0):
        '''
        Creates the trans structure of the azo compound from a SMILES string. 3D geometry creation is done using openbabel. 
        It sets up the trans isomer (including the creation of Molecule_azo and its 'initial' State objects) and stores it as
        source of the System_azo object. 

        To ensure the geometry could be treated a posteriori to create cis or transition states geometries, the SMILES string
        should follow the template: c1(...1)/N=N/ring2. 
            
            e.g. c1(ccccc1)/N=N/c2ccccc2

        If no System_azo is provided, it returns labels, coord and the Molecule_azo named 'trans'.  

        Workflow
        --------
            - Find or create trans isomer.
            - Get indices of atoms relevant to dihedral adjustment.
            - Move the azo group dihedral angle close to the target angle.
            - If there is steric hindrance, tries to solve it by moving adjacent rings.
            - If steric hindrance cannot be solved, returns None.
            - If it is solved, structure is saved to a Molecule_azo object with an 'initial' State. 

        Parameters
        ----------
            self: System_azo object
                System_azo object where the trans isomer will be stored.
            overwrite: bool
                In case the source already has a "trans" isomer, whether it should overwrite it. 
            debug: int
                Debug level. 0: no debug, 1: verbose debug

        Returns
        -------
            trans: Molecule_azo
                Molecule_azo object containing the trans isomer.
        '''

        if debug > 0: print(f"AZO.CREATE_TRANS: Creating trans isomer for {self.name}")
        smiles = self.smiles

        if '/N=N\\' in smiles:
            if debug > 0: print(f"AZO.CREATE_TRANS: Received SMILES of CIS isomer, replacing '\\' with '/' for {self.name}")
            smiles = smiles.replace('/N=N\\', '/N=N/')

        labels, coord          = get_3D(smiles)
        coord                  = centercoords(coord, 0)
        trans                  = Molecule_azo(labels, coord)
        trans.smiles           = smiles
        trans.dihedral_indices = self.dihedral_indices
        if trans.check_fragmentation():
            print(f"AZO.CREATE_TRANS: Trans isomer for {self.name} is FRAGMENTED.")
            return None

        trans.set_total_charge(0)
        trans.set_total_spin(0)            

        trans_state = trans.add_state("initial")
        trans_state.set_geometry(labels, coord)
        self.add_source(name="trans", new_source = trans, overwrite=overwrite)
        return trans

    ######
    def create_cis(self, target_deg: float=40.0, max_iter: int=2000, overwrite: bool=False, debug: int=0) -> "Molecule_azo":
        '''
        Creates the cis structure of the azo compound from a SMILES string. To avoid troubles in 3D geometry creation using openbabel, 
        the trans isomer is created using create_trans() function. It sets up the cis isomer (including the creation of Molecule_azo and 
        its 'initial' State objects) and stores it as source of the System_azo object. 

        If no System_azo is provided, it returns labels, coord and the Molecule_azo named 'cis'.  

        Workflow
        --------
        - Find or create trans isomer.
        - Get indices of atoms relevant to dihedral adjustment.
        - Move the azo group dihedral angle close to the target angle.
        - If there is steric hindrance, tries to solve it by moving adjacent rings.
        - If steric hindrance cannot be solved, returns None.
        - If it is solved, structure is saved to a Molecule_azo object with an 'initial' State. 

        Parameters
        ----------
        target_deg: float
            Target angle of the azo dihedral in degrees. 
        max_iter: int
            Maximum number of iterations to solve steric hindrance.
        overwrite: bool
            In case the source already has a "cis" isomer, whether it should overwrite it. 
        debug: int
            Debug level. 0: no debug, 1: verbose debug

        Returns
        -------
        iso: Molecule_azo
            Molecule_azo object containing the cis isomer.
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
            coord                = centercoords(coord, at1)
            cis                  = Molecule_azo(labels, coord) 
            cis.smiles           = trans.smiles.replace('/N=N/','/N=N\\')
            cis.dihedral_indices = self.dihedral_indices

            cis.set_total_charge(0)
            cis.set_total_spin(0)

            self.add_source(name="cis", new_source = cis, overwrite=overwrite)

            if cis.check_fragmentation():
                print(f"AZO.CREATE_CIS: Cis isomer for {self.name} is FRAGMENTED.")
                return None
            
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
        By default, total spin is set as 1. It can be changed using the function for Molecule objects as Molecule.set_total_spin(value).
        Two versions are created:
            - TSinv_L: Inversion TS involving inversion of left ring (at0 - at1 = at2). Total spin is set as 1.
            - TSinv_R: Inversion TS involving inversion of right ring (at1 = at2 - at3). Total spin is set as 1.

        They are added as sources of the azosystem. They can be accessed by using azosystem.find_source('TSrot_A_T') or azosystem.find_source('TSinv_R')
        
        """
        created_ts = []
        # 1st-searching for trans isomer
        trans_found, trans = self.find_source('trans')
        if not trans_found: raise Exception(f"AZO.CREATE_TS: [ERROR] Trans isomer not found. Create it first with self.create_trans()")

        # Get trans isomer geometry as reference
        labels = trans.labels
        coord = trans.coord

        # Get indices for the azo group (at1 - at2 = at3 - at4) and neighbours (at0, at5)
        at0, at1, at2, at3, at4, at5 = self.dihedral_indices
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
                ts = Molecule_azo(labels, coord)
                isFragmented = ts.check_fragmentation()  # Check if the TSrot is fragmented
                if not isFragmented:
                    state = ts.add_state("initial")
                    state.set_geometry(labels, coord)

                    ts.set_total_charge(0)
                    ts.set_total_spin(0)

                    ts.dihedral_indices = self.dihedral_indices
                    self.add_source('TSrot_A_S', ts)
                    created_ts.append(ts)

                    if 'triplet' in ts_list:
                        ts_triplet = Molecule_azo(labels, coord)

                        ts_triplet.set_total_charge(0)
                        ts_triplet.set_total_spin(2)

                        ts_triplet.dihedral_indices = self.dihedral_indices
                        triplet_state = ts_triplet.add_state("initial")
                        triplet_state.set_geometry(labels, coord)
                        self.add_source('TSrot_A_T', ts_triplet)
                        created_ts.append(ts_triplet)
                        if debug > 0:   print(f'AZOS.CREATE_TS.TSROT_A: TSrot_A_T Molecule_azo successfully created for {self.name}')
                    if debug > 0:       print(f'AZOS.CREATE_TS.TSROT_A: TSrot_A_S Molecule_azo successfully created for {self.name}')
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
                ts = Molecule_azo(labels, coord)
                isFragmented = ts.check_fragmentation()  # Check if the TSrot is fragmented
                if not isFragmented:

                    ts.set_total_charge(0)
                    ts.set_total_spin(0)

                    ts.dihedral_indices = self.dihedral_indices
                    state = ts.add_state("initial")
                    state.set_geometry(labels, coord)
                    self.add_source('TSrot_B_S', ts)
                    created_ts.append(ts)

                    if 'triplet' in ts_list:
                        ts_triplet = Molecule_azo(labels, coord)

                        ts_triplet.set_total_charge(0)
                        ts_triplet.set_total_spin(2)

                        ts_triplet.dihedral_indices = self.dihedral_indices
                        triplet_state = ts_triplet.add_state("initial")
                        triplet_state.set_geometry(labels, coord)
                        self.add_source('TSrot_B_T', ts_triplet)
                        created_ts.append(ts_triplet)
                        if debug > 0:   print(f'AZOS.CREATE_TS.TSROT_B: TSrot_B_T Molecule_Azo successfully created for {self.name}')
                    if debug > 0:       print(f'AZOS.CREATE_TS.TSROT_B: TSrot_B_S Molecule_Azo successfully created for {self.name}')
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
                    ts = Molecule_azo(labels, coord)
                    ts_isFragmented = ts.check_fragmentation()  # Check if the TSinv_l is fragmented
                    if not ts_isFragmented:

                        ts.set_total_charge(0)
                        ts.set_total_spin(0)

                        ts.dihedral_indices = self.dihedral_indices
                        ts_state = ts.add_state("initial")
                        ts_state.set_geometry(labels, coord)
                        self.add_source('TSinv_l', ts)
                        created_ts.append(ts)
                        if debug > 0: print(f'AZOS.CREATE_TS.TSINV_L: TSinv_l Molecule_azo successfully created for {self.name}')
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
                    ts = Molecule_azo(labels, coord)
                    ts_isFragmented = ts.check_fragmentation()  # Check if the TSinv_l is fragmented
                    if not ts_isFragmented:

                        ts.set_total_charge(0)
                        ts.set_total_spin(0)

                        ts.dihedral_indices = self.dihedral_indices
                        ts_state = ts.add_state("initial")
                        ts_state.set_geometry(labels, coord)
                        self.add_source('TSinv_r', ts)
                        created_ts.append(ts)
                        if debug > 0: print(f'AZOS.CREATE_TS.TSINV_R: TSinv_r Molecule_azo successfully created for {self.name}')
                        break
                    else:
                        raise Exception(f'AZOS.CREATE_TS.TSINV_R: [ERROR] TSinv_r fragmented for {self.name}')
        return created_ts

    def plot_energies(self):
        import matplotlib.pyplot as plt
        import numpy as np

        plt.rcParams.update({
            'font.size': 12,
            'font.family': 'sans-serif',
            'axes.linewidth': 1.5,
            'xtick.major.width': 1.5,
            'ytick.major.width': 1.5
        })

        # --- SOURCE DETECTION ---
        # Find Transition States
        ts_candidates = [i for i, source in enumerate(self.sources) 
            if source.name.lower().startswith('ts') 
            and "Gtot" in source.find_state('opt')[1].results]
        
        # Find Isomers
        isomer_candidates = [i for i, source in enumerate(self.sources) 
            if (source.name.lower().startswith('cis') or source.name.lower().startswith('trans')) 
            and "Gtot" in source.find_state('opt')[1].results]
        candidates_idx = ts_candidates + isomer_candidates
        if not candidates_idx:
            print("No valid isomers or TS structures with 'Gtot' found.")
            return

        energies = []
        labels = []

        # --- DATA EXTRACTION ---
        for i in candidates_idx:
            source = self.sources[i]
            opt_results = source.find_state('opt')[1].results
            
            if source.spin == 2 and "Gtot_corr" in opt_results:
                e_val = opt_results["Gtot_corr"].value
            else:
                e_val = opt_results["Gtot"].value
                print('AZO.PLOT_ENERGIES: WARNING! A source in the triplet state does not have its Free Energy corrected.')
                
            energies.append(e_val)
            labels.append(source.name)

        energies = np.array(energies)
        labels = np.array(labels)

        # --- CLASSIFICATION & COORDINATES ---
        is_trans = np.char.startswith(np.char.lower(labels), "trans")
        is_cis = np.char.startswith(np.char.lower(labels), "cis")

        # X-axis: TS default to 2.0
        x_coords = np.full_like(energies, 2.0)
        x_coords[is_trans] = 1.0  # E isomer
        x_coords[is_cis] = 3.0    # Z isomer

        # --- NORMALIZATION ---
        if np.any(is_trans):
            # Use the first trans isomer as the reference zero
            ref_idx = np.where(is_trans)[0][0]
            ref_energy = energies[ref_idx]
            
            energies = (energies - ref_energy) * Constants.har2kJmol * Constants.kJmol2kcal  # Hartree to kcal/mol
            ylabel = r"$\Delta G$ (kcal mol$^{-1}$)"
        else:
            ylabel = "Gibbs Energy (Hartree)"

        # --- HIGHLIGHT LOGIC ---
        ts_indices = (x_coords == 2.0)
        # Added a check to ensure we don't pass an empty array to np.min
        min_ts_energy = np.min(energies[ts_indices]) if np.any(ts_indices) else None

        # --- PLOTTING ---
        fig, ax = plt.subplots(figsize=(7, 5))

        for x, e, name in zip(x_coords, energies, labels):
            
            color = 'black'
            font_weight = 'normal'
            
            # Lowest Energy TS is highlighted
            if x == 2.0 and min_ts_energy is not None and np.isclose(e, min_ts_energy):
                color = '#D55E00'  
                font_weight = 'bold'

            ax.plot(x, e, marker='_', markersize=50, markeredgewidth=3, color=color)
            
            clean_name = name.replace("trans", "").replace("cis", "").strip("-_ ")
            ax.text(x + 0.4, e - 0.5, clean_name, 
                    ha='center', va='bottom', 
                    fontsize=10, color=color, weight=font_weight)

        # --- AXIS FORMATTING ---
        ax.set_ylabel(ylabel)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(['E', 'TS', 'Z'])
        ax.set_xlim(0.5, 3.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        plt.show()
        
        print("----SUMMARY OF ENERGY BARRIERS----")
        for name, energy in zip(labels,energies):
            print(name, ' ', energy)


    def get_PSS(self, lamp : "Lamp", phi_EZ = 0.3, phi_ZE = 0.5, t_EZ=None, t_ZE=None, debug=0):

        # Search for cis and trans sources.
        found_cis_opt, cis_state = self.find_source('cis')[1].find_state('opt')
        found_trans_opt, trans_state = self.find_source('trans')[1].find_state('opt')

        if not (found_cis_opt and found_trans_opt):
            raise Exception('SYSTEM_AZO.GET_PSS: One of the optimized states were not found.')

        print(f'SYSTEM_AZO.GET_PSS: TRANS: {trans_state.results}')
        print(f'SYSTEM_AZO.GET_PSS: CIS: {cis_state.results}')
        return get_PSS(cis_state, trans_state, lamp = lamp, phi_EZ = phi_EZ, phi_ZE = phi_ZE, t_EZ=t_EZ, t_ZE=t_ZE, debug=debug)
    
    def get_abs_spectrum(self, normalize: bool = False, units: bool = False, get_PSS: bool = False,
                         custom_cis= None, custom_trans= None):

        if custom_cis is None:
            found_cis, cis = self.find_source('cis')
        if custom_trans is None:
            found_trans, trans = self.find_source('trans')

        if not found_cis or not found_trans:
            print(f'Z_exists: {found_cis}, E_exists: {found_trans}')
            raise ValueError(f'Error: No cis or trans conformers found for system {self.name}')

        opt_cis_exists, cis_state = cis.find_state('opt')
        opt_trans_exists, trans_state = trans.find_state('opt')

        if not opt_cis_exists or not opt_trans_exists:
            print(f'opt_Z_exists: {opt_Z_exists}, opt_E_exists: {opt_E_exists}')
            raise ValueError(f'Error: No opt state found for Z or E isomers')

        return get_abs_spectrum(cis_state, trans_state, normalize= normalize, 
            units= units, get_PSS= get_PSS)

    def compute_pss(self, wl_list, state_trans=None, state_cis=None, phi_EZ=0.3, phi_ZE=0.5, t_EZ=None, t_ZE=None, shift_nm=None, debug =0): 
        """
        Calculates the photostationary state for a A <--> B photo-interconversion.
        """
       
        if state_trans is None: 
            state_trans = self.find_source('trans')[1].find_state('opt')[1]
        if state_cis is None: 
            state_cis = self.find_source('cis')[1].find_state('opt')[1]

        if debug>0: print(f'SYSTEM_AZO.COMPUTE_PSS: Trans results: {state_trans.results}')
        if debug>0: print(f'SYSTEM_AZO.COMPUTE_PSS: Cis results: {state_cis.results}')
        name_trans = state_trans._source.name # Trans
        name_cis = state_cis._source.name # Cis

        frac_list = []
        pss_list = []

        # 1. Initial dark condition
        lamp = Lamp(name='DARK', wavelength=365)
        
        lambda_grid, sigma_cis, sigma_trans, pss_B = self.get_PSS(lamp, phi_EZ=0, phi_ZE=0, t_EZ=t_EZ, t_ZE=t_ZE, debug=debug) 
        
        sigma_cis *= Constants.avogadro / (1000 * np.log(10)) * 1e4
        sigma_trans *= Constants.avogadro / (1000 * np.log(10)) * 1e4

        sigma_pss = build_pss_spectrum(pss_B, sigma_cis, sigma_trans, debug = debug) 
        frac_list.append(pss_B)     
        pss_list.append(sigma_pss)  

        # 2. Iterate through irradiations
        for irr_wl in wl_list:
            lamp = Lamp(name='COOLLED', wavelength=irr_wl, shift_nm=shift_nm)
            _, _, _, pss_B = self.get_PSS(lamp, phi_EZ=phi_EZ, phi_ZE=phi_ZE, t_EZ=t_EZ, t_ZE=t_ZE, debug=debug) 
            sigma_pss = build_pss_spectrum(pss_B, sigma_cis, sigma_trans)
            frac_list.append(pss_B)    
            pss_list.append(sigma_pss)  

        # 3. Store the extracted names AND the extracted half-lives!
        # This saves the plotting function from having to look them up again.
        self.pss_data = {
            "name_trans": name_trans,
            "name_cis": name_cis,
            "lambda_grid": lambda_grid,
            "wl_list": wl_list,
            "sigma_cis": sigma_cis,
            "sigma_trans": sigma_trans,
            "frac_list": frac_list,      
            "pss_list": pss_list,
            # Safely grab the halflife from the State object's results dictionary
            "t_EZ": t_EZ if t_EZ is not None else state_trans.results['halflife'].value, 
            "t_ZE": t_ZE if t_ZE is not None else state_cis.results['halflife'].value   
        }
        print(f"PSS data for {name_trans}/{name_cis} successfully computed!")

    def plot_pss(self, solvent="Unknown", xlim=(250, 600), skip_spectra=False, only_dark=False, savetoimg=False, plots_dir="./plots/"):
        import matplotlib.pyplot as plt
        import os

        if not hasattr(self, 'pss_data'):
            raise AttributeError(f"No PSS data found for {self.name}. Please run .compute_pss() first.")

        data = self.pss_data
        name_trans = data["name_trans"]
        name_cis = data["name_cis"]
        
        # The half-lives were already safely extracted and stored by compute_pss!
        t_EZ_str = format_time(data["t_EZ"])
        t_ZE_str = format_time(data["t_ZE"])
        
        fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
        
        if not skip_spectra:
            ax.plot(data["lambda_grid"], data["sigma_trans"], color='black', linestyle='dashed', linewidth=1, label=f'{name_trans}  t = {t_EZ_str}')
            ax.plot(data["lambda_grid"], data["sigma_cis"], color='blue', linestyle='dashed', linewidth=1, label=f'{name_cis}  t = {t_ZE_str}')

        for i, pss_wl in enumerate(data["pss_list"]):
            frac_A = data["frac_list"][i]
            if i == 0:
                ax.plot(data["lambda_grid"], pss_wl, color='black', linewidth=2, label=f'DARK ({frac_A:.2f} {name_trans})')
                if only_dark: break
            else:
                color = wavelength_to_rgb(data["wl_list"][i-1])
                ax.plot(data["lambda_grid"], pss_wl, linewidth=0.75, c=color, label=f'{data["wl_list"][i-1]} nm ({frac_A:.2f} {name_trans})')

        ax.set_title(f'{self.name} in {solvent}')
        ax.legend(fontsize=9)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel(r'$\epsilon$ (M$^{-1}$ cm$^{-1}$)')
        ax.set_xlim(xlim[0], xlim[1]) 
        plt.tight_layout()

        if savetoimg:
            os.makedirs(plots_dir, exist_ok=True)
            plt.savefig(os.path.join(plots_dir, f"{self.name}_{name_trans}_{name_cis}_pss.png"), dpi=300)
            plt.close(fig)
        else:
            plt.show()

    def __repr__(self):
        to_print = ""
        to_print += f'------------- SCOPE System_azo Object -------\n'
        to_print += f' Name                  = {self.name}\n'
        to_print += f' Atom Indices for Dihedral = {self.dihedral_indices}\n'
        to_print += System.__repr__(self, indirect=True)
        if hasattr(self,"dE"): to_print += f' Thermal Stability     = {self.dE} kJ/mol (- means trans is more stable)\n'
        to_print += '---------------------------------------------\n'
        to_print += '\n'
        return to_print

##################################################################
##### MOLECULE - SPECIE Object Adapted to Molecule_azo Class #####
##################################################################

class Molecule_azo(Molecule):

    def __init__(self, labels, coord):
        Molecule.__init__(self, labels, coord)
        self.subtype  = "molecule_azo"

    def set_halflife_time(self, skip_triplets : bool = True, overwrite = False, debug: int = 0):
        '''
        Computes t0.5 in seconds for a given conformer/isomer stored in a Molecule_azo object e.g. cis or trans using the Eyring equation.
        Saves the result as a data object with the key 'halflife' in the Molecule_azo object.
        Minimum energy Transition State can be accessed using the key 'mets' in the Molecule_azo object. E. g. cis.mets
        The argument skip_triplets is used to whether take into account triplet states, since their optimized geometry could 
        

        Parameters
        ----------
        self : Molecule_azo
            The Molecule_azo object to compute the halflife for.
        skip_triplets : bool
            Skip triplet conformers in halflife calculation.
        overwrite : bool
            Overwrite existing t0.5 values.

        Notes
        -----
        The function will only consider Molecule_azo objects that have an opt state with a Gtot value.
        The function will only consider TSs that have an opt state with a Gtot value, or Gtot_corr value if skip_triplets is False.
        
        '''
        ts_values = []
        ts_names = []

        # Check if triplets are corrected.
        iscorrect = check_triplet(self, overwrite=overwrite, debug=debug) 

        if not iscorrect:
            raise Exception("MOLECULE_AZO.SET_HALFLIFE_TIME: There has been an error with correcting triplet Gtot. Check if computations have finished ")

        found_iso_opt, iso_state = self.find_state("opt")
        if not found_iso_opt:
            raise Exception(f'MOLECULE_AZO.SET_HALFLIFE_TIME: Optimization state not found for {self.name}.')

        # Find Gtot from selected isomer, stored as Molecule_azo object
        if 'Gtot' in iso_state.results.keys(): g_iso = iso_state.results['Gtot'].value
        else: raise ValueError('MOLECULE_AZO.SET_HALFLIFE_TIME: Gtot not found for isomer')

        parent = self._sys
        candidates = [source for source in parent.sources if source.name.lower().startswith('ts')]
        candidates_names = [source.name for source in candidates]

        for ts in candidates:
            found_ts_state, ts_state = ts.find_state("opt")
            if not found_ts_state or not 'Gtot' in ts_state.results.keys():
                print(f'MOLECULE_AZO.SET_HALFLIFE_TIME: [WARNING] Optimization state or Gtot not found for {ts.name}.')
                continue
            # Use only TSs with Gtot or Gtot_corr

            if ts.spin == 2:     # Use corrected Gtot for triplets
                if not skip_triplets:
                    if 'Gtot_corr' in ts_state.results.keys(): energy = ts_state.results['Gtot_corr'].value
                    else: raise ValueError(f'MOLECULE_AZO.SET_HALFLIFE_TIME: Corrected Gtot for {ts.name} Molecule_azo not found for Triplet TS, altough it was corrected with correct_tripletG() function.')
                else:
                    continue
            else:   energy = ts_state.results['Gtot'].value

            name = ts.name 
            ts_names.append(name)
            ts_values.append(energy)
        
        if not hasattr(self, 'halflife') or overwrite:
            if debug > 0: print(rf'MOLECULE_AZO.SET_HALFLIFE_TIME: Collected {len(ts_values)} TSs for {self.name} : {ts_names} with energies {ts_values}.')
            if debug > 0: print(f'MOLECULE_AZO.SET_HALFLIFE_TIME:Doing halflife for {parent.name} {self.name}')
            
            # Choosing Minimum Energy TS (mets)
            min_idx = int(np.argmin(ts_values))
            mets = ts_names[min_idx]        # Name of the METS 
            g_cross = ts_values[min_idx]
            dG_cross = (float(g_cross)- float(g_iso)) * Constants.har2kJmol * 0.24 # in Kcal/mol
            
            # Compute and store halflife
            t,k = compute_t(float(g_cross), float(g_iso))           
            if debug > 0: print(t, ' s for isomer ', self.name)
            new_time = Data('halflife', float(t), 's', 'compute_t')
            self.halflife = float(t)
            self.mets = mets
            iso_state.add_result(new_time, overwrite=overwrite)
            
            if debug > 0: print(dG_cross, 'for isomer ', self.name)
            self.dG_cross = dG_cross        # Stored in kcal/mol
            newdata = Data('dG_cross', float(dG_cross), 'kcal/mol', 'set_halflife_time')
            iso_state.add_result(newdata, overwrite=overwrite)
        else: 
            print(f'MOLECULE_AZO.SET_HALFLIFE_TIME: State not found in {self.name}.')
            
        print(f'MOLECULE_AZO.SET_HALFLIFE_TIME: Half-life time has been computed for {self.name} and was stored as a result to the corresponding State!')
        print(f'MOLECULE_AZO.SET_HALFLIFE_TIME: ---- Attributes such as dG_cross and mets can be accessed by self.dG_cross and self.mets)')
        

    def link_tda_to_state(self, state: object, filepath: str, overwrite: bool=False):
        '''
        Parses excited states from a TDDFT Gaussian16 .log file to a ExcitedState object.
        '''
        if not hasattr(state,"es_list") or overwrite: 
            if os.path.exists(filepath):
                state.tda_filepath = filepath
                state.es_list = []
                lines = read_lines_file(filepath)
                # state.gs_energy = parse_energy_from_step(lines)
                if parse_status_finished(lines):
                    for st_num in range(1,11):
                        if st_num < 10:    
                            line_nums, found = search_string(f"Excited State   {st_num}:", lines, typ='first')
                        elif st_num >= 10 and st_num < 100: 
                            line_nums, found = search_string(f"Excited State  {st_num}:", lines, typ='first')
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
                output= G16_output(lines)
                opt_finished = output.get_optimization_finished()
                if opt_finished:
                    state.energy = output.get_last_energy()
                    state.gtot = output.get_free_energy()

                    labels, coord = output.get_last_geometry(lines)
                    state.set_geometry(labels, coord)
                    newG = Data("Gtot", state.gtot, "au", "get_free_energy")
                    newH = Data("energy", state.gtot, "au", "get_last_energy")
                    state.add_result(newG)
                else:
                    print('LINK_OPT_TO_STATE: WARNING: An optimization did not finished')
                    print(f'LINK_OPT_TO_STATE: Optimization file: {filepath}')


        else:
            print(f"File {filepath} does not exist.")


    def get_azo_substituents(self, debug: int=0):
        self.azo_substituents = []
        azo_idx   = self.dihedral_indices[2:4]
        rest_idx  = list(idx for idx in self.indices if idx not in azo_idx)
        rest_indices = extract_from_list(rest_idx, self.indices, dimension=1)
        rest_labels  = extract_from_list(rest_idx, self.labels, dimension=1)
        rest_coord   = extract_from_list(rest_idx, self.coord, dimension=1)
        rest_radii   = extract_from_list(rest_idx, self.radii, dimension=1)
        rest_atoms   = extract_from_list(rest_idx, self.atoms, dimension=1)
        blocklist = split_species(rest_labels, rest_coord)
        for b in blocklist:
            if debug > 0: print(f"GET_AZO_SUBSTITUENTS. PREPARING BLOCK: {b}")
            sub_indices      = extract_from_list(b, rest_indices, dimension=1)
            sub_labels       = extract_from_list(b, rest_labels, dimension=1)
            sub_coord        = extract_from_list(b, rest_coord, dimension=1)
            sub_radii        = extract_from_list(b, rest_radii, dimension=1)
            sub_atoms        = extract_from_list(b, rest_atoms, dimension=1)
            new_substituent  = Molecule(sub_labels, sub_coord, radii=sub_radii)
            new_substituent.origin = "get_azo_substituents"
            new_substituent.add_parent(self, indices=sub_indices)
            new_substituent.set_atoms(atomlist=sub_atoms)
            self.azo_substituents.append(new_substituent)
        return self.azo_substituents

    def __repr__(self):
        to_print = ""
        to_print += f'---------- SCOPE Molecule_azo Object ------------\n'
        to_print += Molecule.__repr__(self, indirect=True)
        # if hasattr(self,"halflife"): to_print += f'     Half-life (t1/2):{self.halflife}                         \n'
        to_print +=  '-------------------------------------------------\n'
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
        
        fwhm_list = {       # Fwhm data (nm) from COOLLED light source
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

        if name=='DARK':
            self.wavelength = 365
            self.irradiance = 0 
            self.power = 0


    def __repr__(self):
        to_print = f'--------- Scope Lamp Object -----------\n'
        to_print += f' Name                     = {self.name}\n'
        to_print += f' Wavelength               = {self.wavelength} nm\n'
        to_print += f' FWHM                     = {self.fwhm} nm\n'
        if hasattr(self, "eff_wavelength"): to_print += f' Wavelength (after shift) = {self.eff_wavelength} nm\n'
        if hasattr(self, "power"):          to_print += f' Power                    = {self.power} W\n'
        to_print += f'---------------------------------------\n'
        return to_print
