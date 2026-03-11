import numpy as np
from scope import constants                    
from scope.classes_data                  import Data, Collection
from scope.classes_system                import System
from scope.classes_state                 import State, find_state
from scope.classes_specie                import Molecule
from scope.elementdata                   import ElementData
from scope.geometry                      import *
from scope_azo.azo_functions import * 

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
        self.results          = dict()

    #################
    #### Results ####
    #################
    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  
            self.results[result.key] = result

    def remove_result(self, key: str):
        return self.results.pop(key, None)

    ######
    def add_source(self, name: str, new_source: object, overwrite: bool=False, debug: int=0):
        '''
        Same function than vanilla System class. It is copied here to make sure that set_initial_state() calls it.
        '''

        ## Source names are de-capitalized and spaces replaced by underscores
        name = name.lower()
        name = name.replace(" ","_")
        ## Sources Must Have a Name
        if not hasattr(new_source,"name"): new_source.name = name
        ## Search if source with the same name already exists
        found, old_source = self.find_source(name, debug=debug)
        ## If not, it is added
        if not found: 
            new_source._sys = self ## Links the system to the source 
            new_source.set_initial_state(debug=debug)
            self.sources.append(new_source)
        ## If it exists, it is overwritten if specified 
        elif found and overwrite: 
            new_source._sys = self ## Links the system to the source 
            new_source.set_initial_state(debug=debug)
            self.sources = [s for s in self.sources if s.name.lower() != name.lower()]
            self.sources.append(new_source)
        else: 
            print(f"SYSTEM_AZO.ADD_SOURCE: Source with name '{new_source.name}' already exists in system '{self.name}'") 
            print(f"SYSTEM_AZO.ADD_SOURCE: If you would like to Overwrite, specify overwrite=True")
        return self.sources

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

    ############################
    ## Creation of Structures ##
    ############################
    def create_trans(self, charge: int=0, overwrite: bool=False, debug: int=0):
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
            self: System_azo        System_azo object where the trans isomer will be stored.
            overwrite: bool         In case the source already has a "trans" isomer, whether it should overwrite it. 
            debug: int              Debug level. 0: no debug, 1: verbose debug

        Returns
        -------
            trans: Molecule_azo     Molecule_azo object containing the trans isomer.
        '''
        from scope_azo.azo_functions import get_3D

        if debug > 0: print(f"AZO.CREATE_TRANS: Creating trans isomer for {self.name}")
        smiles = self.smiles

        if '/N=N\\' in smiles:
            if debug > 0: print(f"AZO.CREATE_TRANS: Received SMILES of CIS isomer, replacing '\\' with '/' for {self.name}")
            smiles = smiles.replace('/N=N\\', '/N=N/')

        labels, coord          = get_3D(smiles)
        coord                  = centercoords(coord, 0)
        trans                  = Molecule_azo(labels, coord)

        # Checks fragmentation to verify that the 3D structure is correct (or, at least, not fragmented)
        if trans.check_fragmentation():
            print(f"AZO.CREATE_TRANS: Trans isomer for {self.name} is FRAGMENTED.")
            return None

        # Sets charge and spin
        trans.set_total_charge(charge)
        #trans.set_total_spin(0)                    ## Not needed, it should be the default      

        # Sets other attributes
        trans.smiles           = smiles
        trans.dihedral_indices = self.dihedral_indices

        # Adds to System as source. Initial State is created automatically when sourcing
        self.add_source("trans",trans,overwrite=overwrite) 
        return trans

    ######
    def create_cis(self, target_deg: float=40.0, overwrite: bool=False, debug: int=0):
        '''
        Creates the cis structure of the azo compound from a SMILES string. To avoid troubles in 3D geometry creation using openbabel, 
        The trans isomer is first created using create_trans() function. Then, it sets up the cis isomer and stores it as source of the System_azo object. 

        Workflow
        --------
        - Find or create trans isomer.
        - Get indices of atoms relevant to dihedral adjustment.
        - Move the azo dihedral angle to the target angle.
        - If there is steric hindrance, tries to solve it by moving adjacent rings.
        - If steric hindrance cannot be solved, returns None.
        - If it is solved, structure is saved to a Molecule_azo object with an 'initial' State. 

        Parameters
        ----------
        target_deg: float       Target dihedral angle (in degrees). Default is 40º 
        overwrite: bool         In case the System already has a "cis" source, whether it should overwrite it. 
        debug: int              Debug level. 0: no, 1: verbose 

        Returns
        -------
        cis: Molecule_azo       Molecule_azo object containing the cis isomer.

        '''
        from scope.connectivity import get_adjmatrix

        # 1st-searching for trans isomer
        found, trans = self.find_source('trans')
        if not found: trans = self.create_trans(debug=debug)
        
        # Get trans isomer geometry as reference
        labels = trans.labels
        coord  = trans.coord

        # Get indices for the azo group (at1 - at2 = at3 - at4) and neighbours (at0, at5)
        at0, at1, at2, at3, at4, at5 = self.get_dihedral_indices()
        if debug > 0: print(f"AZO.CREATE_CIS: dihedral indices: {at0}, {at1}, {at2}, {at3}, {at4}, {at5}")

        # Get the adjacency matrix. Used as a reference when rotating the dihedral angle
        adjmat_ref, adjnum_ref = trans.get_adjmatrix()

        # Initial dihedral angle, for info 
        current_deg = np.degrees(get_dihedral(coord[at1],coord[at2],coord[at3],coord[at4])) # Initial dihedral angle
        if debug > 0: print(f"AZO.CREATE_CIS: Initial dihedral angle: {current_deg} degrees")

        # Pre-Rotation. Setting dihedral to target
        if debug > 0: print(f"AZO.CREATE_CIS: Jumping to {target_deg} degrees. Selected atoms: {at1}, {at2}, {at3}, {at4}")
        coord_next = set_dihedral(labels, coord, target_deg, at1,at2,at3,at4,adjmat=adjmat_ref, adjnum=adjnum_ref)
        angle_next = np.degrees(get_dihedral(coord_next[at1],coord_next[at2],coord_next[at3],coord_next[at4]))
        if debug > 0:  print(f"AZO.CREATE_CIS: Changed dihedral in {self.name} from {current_deg} to {angle_next}, it should be near +-{target_deg}")

        # Get the adjacency matrix after rotating dihedral for comparison.
        _, adjmat_cis, adjnum_cis = get_adjmatrix(labels,coord_next)

        # Ensure atomic connectivity remains intact after dihedral modification. 
        matrices_match = np.array_equal(adjmat_cis, adjmat_ref) and np.array_equal(adjnum_cis, adjnum_ref)
        if debug > 0: print(f"AZO.CREATE_CIS: Matrices match: {matrices_match}")
        if matrices_match:
            coord = coord_next
            current_deg = angle_next
            if debug > 0: print(f"AZO.CREATE_CIS: Angle {current_deg:.2f} (OK). Saving it as source")

        # Connectivity change indicates steric clashes. Attempt to resolve by rotating adjacent rings.
        else:
            fixed_collision, coord = solve_dihedral(labels, coord_next, at0, at1, at2, at3, at4, at5, adjmat_ref=adjmat_ref, adjnum_ref=adjnum_ref,debug=debug)
            current_deg = np.degrees(get_dihedral(coord[at1],coord[at2],coord[at3],coord[at4]))
            if debug > 0: print(f"AZO.CREATE_CIS: Fixed collision: {fixed_collision}")
            # If the collision is fixed by rotating adjacent rings and the target dihedral remains almost the same 
            if abs(current_deg) <= (target_deg + 5) and fixed_collision:
                if debug > 0: print(f'AZO.CREATE_CIS: Found good geometry for {self.name} with angle {current_deg:.2f} by rotating adjacent dihedrals')
            else: 
                if debug > 0: print(f'AZO.CREATE_CIS: Failed to find good geometry for {self.name}.')
                return None

        # Here, it means it found a good geometry for the cis isomer, without clashes
        coord = centercoords(coord, at1)
        cis   = Molecule_azo(labels, coord) 

        # Checks fragmentation to verify that the 3D structure is correct (or, at least, not fragmented)
        if cis.check_fragmentation(debug=debug):
            print(f"AZO.CREATE_CIS: Cis isomer for {self.name} is FRAGMENTED.")
            return None

        # Sets other attributes
        cis.smiles           = trans.smiles.replace('/N=N/','/N=N\\')
        cis.dihedral_indices = self.dihedral_indices
        cis.set_total_charge(trans.charge)
        cis.set_total_spin(trans.spin)

        # Adds to System as source. Initial State is created automatically when sourcing
        self.add_source("cis",cis,overwrite=overwrite)

        return cis

    ######
    def create_ts(self, ts_list:list = ['TSrot', 'TSinv', 'triplet'], debug: int=0):
        """
        Creates a set of TS for a System_azo. Users can select which TS must be created from a list of options (['TSrot', 'TSinv', 'triplet']). 
        Once the TS are created, they are added to the System_azo object as sources.

        -------
        TSrot
        -------
        TSrot is created by rotating the trans isomer by +/-90 degrees around the azo dihedral. 

        Nomenclature: 
            - TSrot_A_S: Rotation TS in singlet state. Dihedral angle is set as +90º. Total spin is set to 0.
            - TSrot_A_T: Rotation TS in triplet state. Dihedral angle is set as +90º. Total spin is set as 2.
            - TSrot_B_S: Rotation TS in singlet state. Dihedral angle is set as -90º. Total spin is set as 0.
            - TSrot_B_T: Rotation TS in triplet state. Dihedral angle is set as -90º. Total spin is set as 2.

        -------
        TSinv
        -------
        TSinv are created by rotating the azo N=N-ring angles to 180º.
        By default, total spin is set as 0. It can be changed using the function for Molecule objects as Molecule.set_total_spin(value).
        2 versions are created:
            - tsinv_l: Inversion TS involving inversion of left ring (at0 - at1 = at2).
            - tsinv_r: Inversion TS involving inversion of right ring (at1 = at2 - at3).

        -------
        Triplets
        -------
        Similar to the TSrot, but setting the molecule in a triplet state. Technically, the final structures won't be a TS, but a triplet minimum in the PES. 
        However, it becomes a sort of 'TS' during the cis-trans thermal relaxation involving the Triplet manifold, hence the nomenclature.

        """
        from scope.connectivity import get_adjmatrix

        created_ts = []
        # 1st-searching for trans isomer
        trans_found, trans = self.find_source('trans')
        if not trans_found: raise Exception(f"AZO.CREATE_TS: [ERROR] Trans isomer not found. Create it first with self.create_trans()")

        # Get trans isomer geometry as reference
        labels = trans.labels
        coord = trans.coord

        # Get indices for the azo group (at1 - at2 = at3 - at4) and neighbours (at0, at5)
        at0, at1, at2, at3, at4, at5 = self.dihedral_indices
        if debug > 0: print(f"SYSTEM_AZO.CREATE_CIS: dihedral indices: {at0}, {at1}, {at2}, {at3}, {at4}, {at5}")

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
                    ts.dihedral_indices = self.dihedral_indices
                    self.add_source('TSrot_A_S', ts)
                    created_ts.append(ts)
                    if debug > 0:       print(f'SYSTEM_AZO.CREATE_TS: TSrot_A_S Molecule_azo successfully created for {self.name}')
                    if 'triplet' in ts_list:
                        ts_triplet = Molecule_azo(labels, coord)
                        ts_triplet.set_total_spin(2)
                        ts_triplet.dihedral_indices = self.dihedral_indices
                        self.add_source('TSrot_A_T', ts_triplet)
                        created_ts.append(ts_triplet)
                        if debug > 0:   print(f'SYSTEM_AZO.CREATE_TS: TSrot_A_T Molecule_azo successfully created for {self.name}')
                else:
                    raise Exception(f'SYSTEM_AZO.CREATE_TS: [ERROR] TSrot_A fragmented for {self.name}')

            ## TSrot_B ##
            if debug > 0:
                dg_deg = np.degrees(get_dihedral(trans.coord[at1], trans.coord[at2], trans.coord[at3], trans.coord[at4]))
                print(f'SYSTEM_AZO.CREATE_TS: Dihedral angle for reference geometry: {dg_deg} degrees')
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
                    ts.dihedral_indices = self.dihedral_indices
                    self.add_source('TSrot_B_S', ts)
                    created_ts.append(ts)
                    if debug > 0:       print(f'SYSTEM_AZO.CREATE_TS: TSrot_B_S Molecule_Azo successfully created for {self.name}')
                    if 'triplet' in ts_list:
                        ts_triplet = Molecule_azo(labels, coord)
                        ts_triplet.set_total_spin(2)
                        ts_triplet.dihedral_indices = self.dihedral_indices
                        self.add_source('TSrot_B_T', ts_triplet)
                        created_ts.append(ts_triplet)
                        if debug > 0:   print(f'SYSTEM_AZO.CREATE_TS: TSrot_B_T Molecule_Azo successfully created for {self.name}')
                else:
                    raise Exception(f'SYSTEM_AZO.CREATE_TS: WARNING: TSrot_B fragmented for {self.name}')
        
        ## TSinv Left (TSinv_l) ##
        if 'TSinv' in ts_list:
            init_coord = set_angle(labels, trans.coord, 179.9, at1,at2,at3)
            if debug > 0: 
                angle_deg = np.degrees(get_angle_points(coord[at1], coord[at2], coord[at3]))
                print(f'SYSTEM_AZO.CREATE_TS: Angle between {at1}, {at2} and {at3} set to {angle_deg} degrees.')
            
            ## TSinv_l_a (Left adjacent dihedral angle starting at 0º) 
            angles = np.concatenate(([0],[val for i in range(1, 12) for val in (15 * i, -15 * i)]))
            for a0 in angles:
                coord = set_dihedral(labels, init_coord, a0, at3,at2,at1,at0, adjmat=adjmat_ref, adjnum=adjnum_ref)
                _, adjmat, adjnum = get_adjmatrix(labels,coord)
                is_equal = np.array_equal(adjmat, adjmat_ref) and np.array_equal(adjnum, adjnum_ref)
                if is_equal:
                    ts = Molecule_azo(labels, coord)
                    ts_isFragmented = ts.check_fragmentation()  # Check if the TSinv_l is fragmented
                    if not ts_isFragmented:
                        ts.set_total_charge(0)
                        ts.dihedral_indices = self.dihedral_indices
                        self.add_source('TSinv_l_a', ts)
                        created_ts.append(ts)
                        if debug > 0: print(f'SYSTEM_AZO.CREATE_TS: TSinv_l_a Molecule_azo successfully created for {self.name}')
                        break
                    else:
                        raise Exception(f'SYSTEM_AZO.CREATE_TS: [ERROR] TSinv_l_a fragmented for {self.name}')

            ## TSinv_l_b (Left adjacent dihedral angle starting at 180º) 
            angles = np.insert(angles[::-1],0, 180.)
            if debug >0: print(f'SYSTEM_AZO.CREATE_TS: Angles for TSinv_l_b: {angles}')
            for a0 in angles:
                coord = set_dihedral(labels, init_coord, a0, at3,at2,at1,at0, adjmat=adjmat_ref, adjnum=adjnum_ref)
                _, adjmat, adjnum = get_adjmatrix(labels,coord)
                is_equal = np.array_equal(adjmat, adjmat_ref) and np.array_equal(adjnum, adjnum_ref)
                if is_equal:
                    ts = Molecule_azo(labels, coord)
                    ts_isFragmented = ts.check_fragmentation()  # Check if the TSinv_l is fragmented
                    if not ts_isFragmented:
                        ts.set_total_charge(0)
                        ts.dihedral_indices = self.dihedral_indices
                        self.add_source('TSinv_l_b', ts)
                        created_ts.append(ts)
                        if debug > 0: print(f'SYSTEM_AZO.CREATE_TS: TSinv_l_b Molecule_azo successfully created for {self.name}')
                        break
                    else:
                        raise Exception(f'SYSTEM_AZO.CREATE_TS: [ERROR] TSinv_l_b fragmented for {self.name}')
        
        ## TSinv Right (TSinv_r) ##
            coord = set_angle(labels, trans.coord, 179.9, at2,at3,at4) ## Angle value of 179.9 instead of 180 to avoid numerical issues
            if debug > 0: 
                angle_deg = np.degrees(get_angle_points(coord[at2], coord[at3], coord[at4]))
                print(f'SYSTEM_AZO.CREATE_TS: Angle between {at2}, {at3} and {at4} set to {angle_deg} degrees.')
            
            ## TSinv_r_a (Right adjacent dihedral angle starting at 0º)
            angles = np.concatenate(([0],[val for i in range(1, 12) for val in (15 * i, -15 * i)]))
            if debug >0: print(f'SYSTEM_AZO.CREATE_TS: Angles for TSinv_r_a: {angles}')
            for a0 in angles:
                coord = set_dihedral(labels, coord, a0, at2,at3,at4,at5, adjmat=adjmat_ref, adjnum=adjnum_ref)
                _, adjmat, adjnum = get_adjmatrix(labels,coord)
                is_equal = np.array_equal(adjmat, adjmat_ref) and np.array_equal(adjnum, adjnum_ref)
                if is_equal:
                    ts = Molecule_azo(labels, coord)
                    ts_isFragmented = ts.check_fragmentation()  # Check if the TSinv_l is fragmented
                    if not ts_isFragmented:
                        ts.dihedral_indices = self.dihedral_indices
                        self.add_source('TSinv_r_a', ts)
                        created_ts.append(ts)
                        if debug > 0: print(f'SYSTEM_AZO.CREATE_TS.TSINV_R: TSinv_r Molecule_azo successfully created for {self.name}')
                        break
                    else:
                        raise Exception(f'SYSTEM_AZO.CREATE_TS.TSINV_R: [ERROR] TSinv_r fragmented for {self.name}')
            
            # TSinv_r_b (Right adjacent dihedral angle starting at 180º)
            angles = np.insert(angles[::-1],0,180.)
            if debug >0: print(f'SYSTEM_AZO.CREATE_TS: Angles for TSinv_r_b: {angles}')
            for a0 in angles:
                coord = set_dihedral(labels, coord, a0, at2,at3,at4,at5, adjmat=adjmat_ref, adjnum=adjnum_ref)
                _, adjmat, adjnum = get_adjmatrix(labels,coord)
                is_equal = np.array_equal(adjmat, adjmat_ref) and np.array_equal(adjnum, adjnum_ref)
                if is_equal:
                    ts = Molecule_azo(labels, coord)
                    ts_isFragmented = ts.check_fragmentation()  # Check if the TSinv_l is fragmented
                    if not ts_isFragmented:
                        ts.dihedral_indices = self.dihedral_indices
                        self.add_source('TSinv_r_b', ts)
                        created_ts.append(ts)
                        if debug > 0: print(f'SYSTEM_AZO.CREATE_TS.TSINV_R: TSinv_r_b Molecule_azo successfully created for {self.name}')
                        break
                    else:
                        raise Exception(f'SYSTEM_AZO.CREATE_TS.TSINV_R: [ERROR] TSinv_r_b fragmented for {self.name}')
        return created_ts

    #######################
    ## Thermal Stability ##
    #######################
    def get_mets(self, temp: float=298.15, target_state: str='opt', skip_triplets: bool=False, force: bool=False, p_sh: float=0.0002, debug: int=0):
        '''
        Find the minimum-energy transition structure (METS) for thermal isomerization.

        This method scans the available sources in the system and builds a list of
        candidate "barrier" structures at the requested temperature:
        - Singlet (`spin == 0`): source is considered only if the selected
          `target_state` is identified as a transition state (`state.is_ts`).
        - Triplet (`spin == 2`): source is considered only if the selected
          `target_state` is a minimum on the triplet PES (`state.isminimum`) and
          `skip_triplets` is `False`.

        Energies are compared in atomic units (au):
        - Singlets use `Gtot(T)`.
        - Triplets use the corrected free energy `Gtot_eff(T)`.

        The source with the lowest collected energy is stored in `self.mets` and
        returned.

        Parameters
        ----------
        temp : float, default=298.15           Temperature in Kelvin used to retrieve/compute thermal free energies.
        target_state : str, default='opt'      Name of the state to inspect in each source (for example, `'opt'`).
        skip_triplets : bool, default=False    If `True`, triplet minima are ignored and only singlet TS candidates are used.
        force : bool, default=False            If `False` and `self.mets` already exists, the cached value is returned.
                                               If `True`, candidates are recomputed from current sources/results.         
        debug : int, default=0                 Verbosity level for diagnostic prints.

        Returns
        -------
        Molecule_azo                           The source selected as minimum-energy transition structure.

        Raises
        ------
        Exception
            If a required thermal free energy (`Gtot` / `Gtot_eff`) cannot be obtained for a candidate at the requested temperature.
        ValueError
            If no valid candidates are found (for example, empty candidate list).
        '''
        if hasattr(self,"mets") and not force: return self.mets

        ts_names    = []
        ts_energies = []
        for sou in self.sources:
            found, state = sou.find_state(target_state)
            if not found: 
                if debug > 0: print(f'SYSTEM_AZO.GET_METS: Target state not found in source={sou.name}, spin={sou.spin}')
                continue

            # Collects Energy of Sources that behave as TS in the thermal Cis-to-trans relaxation. 
            # This includes Triplet minima. So, technically, these are not actual TS. 
            if sou.spin == 2 and not skip_triplets:     # Use corrected Gtot for triplets
                if hasattr(state,"isminimum"): 
                    if not state.isminimum:
                        if debug > 0: print(f'SYSTEM_AZO.GET_METS: The target state in source {sou.name} is NOT a minimum in the Triplet PES. Skipping')
                        continue
                else: 
                    if debug > 0: print(f'SYSTEM_AZO.GET_METS: The target state in source {sou.name} has no VNMs. Skipping')
                    continue
            elif sou.spin == 0:   
                ## Skips Sources for which we didn't find a TS.
                if hasattr(state,"is_ts"): 
                    if not state.is_ts:
                        if debug > 0: print(f'SYSTEM_AZO.GET_METS: The target state in source {sou.name} is not a TS. Skipping')
                        continue
                else: 
                    if debug > 0: print(f'SYSTEM_AZO.GET_METS: The target state in source {sou.name} has no VNMs. Skipping')
                    continue

            # Searches Gtot entry
            energy = state.get_gtot_eff(temp, p_sh=p_sh, debug=debug).convert_to_units('au').value
            if energy is None: 
                raise Exception (f'SYSTEM_AZO.GET_METS: Could not fint Gtot_eff({temp}) for {sou.name}.')
            ts_names.append(sou.name)
            ts_energies.append(energy)

        #if debug == 1: print(f'SYSTEM_AZO.GET_METS: Collected {len(ts_energies)} TSs for {self.name}')
        if debug >= 1: 
            print(f'SYSTEM_AZO.GET_METS: Collected {len(ts_energies)} TSs for {self.name}. Printing entries:')
            print(f'\tName, Energy:')
            for idx in range(len(ts_energies)):
                print(f'\t{ts_names[idx]} {ts_energies[idx]}')        
            
        # Choosing Minimum Energy TS (mets)
        min_idx     = int(np.argmin(ts_energies)) 
        mets        = ts_names[min_idx]        
        self.mets   = self.find_source(mets)[1]
        if debug > 0: print(f'SYSTEM_AZO.GET_METS: Selected {self.mets.name}')
        return self.mets

    ######
    def get_thermal_stability(self, target_state: str, temp: float=298.15, debug: int=0):
        # Finds Trans and Target State
        found, trans = self.find_source('trans')
        if not found:       raise Exception('SYSTEM_AZO.get_thermal_stability: Trans source not found.')
        found_state, trans_state = trans.find_state(target_state)
        if not found_state: raise Exception(f'SYSTEM_AZO.get_thermal_stability: Target state: {target_state} not found.')

        # Finds Cis and Target State
        found, cis = self.find_source('cis')
        if not found:       raise Exception('SYSTEM_AZO.get_thermal_stability: cis source not found.')
        found_state, cis_state = cis.find_state(target_state)
        if not found_state: raise Exception(f'SYSTEM_AZO.get_thermal_stability: Target state: {target_state} not found.')

        # Getting the energy of the Trans and Cis isomers
        gtot_trans = trans_state.get_gtot_eff(temp=temp, debug=debug).convert_to_units('au').value
        gtot_cis   = cis_state.get_gtot_eff(temp=temp, debug=debug).convert_to_units('au').value

        dE = (gtot_cis - gtot_trans) * constants.har2kJmol/constants.kcal2kJmol ## Returns value in Kcal/mol
        dE_data = Data("dG_cis-trans",dE,'kcal/mol',"system_azo.get_thermal_stability()")
        self.add_result(dE_data, overwrite=True)
        return dE_data

    ######
    def get_trans_halflife_time(self, temp: float=298.15, target_state: str='opt', skip_triplets: bool=False, p_sh: float=0.0002, debug: int=0):
        from scope.thermodynamics import eyring_equation
        '''
        
        '''
        found, source = self.find_source('trans')
        if not found:       raise Exception('SYSTEM_AZO.GET_TRANS_HALFLIFE_TIME: Trans source not found.')
        found_state, ground_state = source.find_state(target_state)
        if not found_state: raise Exception(f'SYSTEM_AZO.GET_TRANS_HALFLIFE_TIME: Target state: {target_state} not found.')

        # Getting the energy of the Trans and Cis isomers
        gtot_trans = ground_state.get_gtot_eff(temp=temp, debug=debug).convert_to_units('au').value

        # Getting the energy of the METS (Minimum Energy Transition State)
        mets = self.get_mets(temp=temp, target_state=target_state, skip_triplets=skip_triplets, p_sh=p_sh, debug=debug)
        found, mets_state = mets.find_state(target_state)
        gtot_mets = mets_state.get_gtot_eff(temp=temp, debug=debug).convert_to_units('au').value
        
        #Computing halflife time (t05) and kinetic constant (k)
        t05, k_th = eyring_equation(gtot_trans, gtot_mets, temp)
        if debug > 0: print(f'SYSTEM_AZO.GET_TRANS_T05: Obtained t05: {t05}  k_th: {k_th}') 
        dG_from_trans = (float(gtot_mets)- float(gtot_trans)) * constants.har2kJmol * constants.kJmol2kcal # in Kcal/mol
        if debug > 0: print(f'SYSTEM_AZO.GET_TRANS_T05: dG_from_trans = {dG_from_trans}') 

        units_time = 's'
        units_kcal = 'kcal/mol'

        # Creates Data-class objects to store results, and sets temperature as property
        k_data   = Data("rate_thermal_trans2cis", k_th, units_time, "system_azo.get_trans_halflife_time()")
        k_data.add_property("temperature", temp)
        dG_data  = Data("dG_from_trans",dG_from_trans,units_kcal,"system_azo.get_trans_halflife_time()")
        dG_data.add_property("temperature", temp)
        t05_data = Data("trans_halflife_time",t05,units_time,"system_azo.get_trans_halflife_time()")
        t05_data.add_property("temperature", temp)
        self.add_result(k_data, overwrite=True)
        self.add_result(dG_data, overwrite=True)
        self.add_result(t05_data, overwrite=True)
        return self.results['trans_halflife_time']
        
    ######
    def get_cis_halflife_time(self, temp: float=298.15, target_state: str='opt', skip_triplets: bool=False, p_sh: float=0.0002, debug: int=0):
        from scope.thermodynamics import eyring_equation
        '''
        
        '''
        found, source = self.find_source('cis')
        if not found:       raise Exception('SYSTEM_AZO.GET_CIS_HALFLIFE_TIME: cis source not found.')
        found_state, ground_state = source.find_state(target_state)
        if not found_state: raise Exception(f'SYSTEM_AZO.GET_CIS_HALFLIFE_TIME: Target state: {target_state} not found.')

        # Getting the energy of the Trans and Cis isomers
        gtot_cis = ground_state.get_gtot_eff(temp=temp, debug=debug).convert_to_units('au').value

        # Getting the energy of the METS (Minimum Energy Transition State)
        mets = self.get_mets(temp=temp, target_state=target_state, skip_triplets=skip_triplets, p_sh=p_sh, debug=debug)
        found, mets_state = mets.find_state(target_state)
        gtot_mets = mets_state.get_gtot_eff(temp=temp, debug=debug).convert_to_units('au').value
        
        #Computing halflife time (t05) and kinetic constant (k)
        t05, k_th = eyring_equation(gtot_cis, gtot_mets, temp)
        if debug > 0: print(f'SYSTEM_AZO.GET_CIS_T05: Obtained t05: {t05}  k_th: {k_th}') 
        dG_from_cis = (float(gtot_mets)- float(gtot_cis)) * constants.har2kJmol * constants.kJmol2kcal # in Kcal/mol
        if debug > 0: print(f'SYSTEM_AZO.GET_CIS_T05: dG_from_cis = {dG_from_cis}') 

        units_time = 's'
        units_kcal = 'kcal/mol'

        # Creates Data-class objects to store results, and sets temperature as property
        k_data   = Data("rate_thermal_cis2trans", k_th, units_time, "system_azo.get_cis_halflife_time()")
        k_data.add_property("temperature", temp)
        dG_data  = Data("dG_from_cis",dG_from_cis,units_kcal,"system_azo.get_cis_halflife_time()")
        dG_data.add_property("temperature", temp)
        t05_data = Data("cis_halflife_time",t05,units_time,"system_azo.get_cis_halflife_time()")
        t05_data.add_property("temperature", temp)
        self.add_result(k_data, overwrite=True)
        self.add_result(dG_data, overwrite=True)
        self.add_result(t05_data, overwrite=True)
        return self.results['cis_halflife_time']

    ########################
    ## Optical Properties ##
    ########################
    def get_PSS(self, lamp_name: str="default", target_state: str = 'opt', temp: float=298.15, phi_EZ: float=0.3, phi_ZE: float=0.5, pw_int: float=1.0, trans_k=None, cis_k=None, lmin: float=200, lmax: float=1000, debug=0):
        """
        Creates and initializes a PSS object for this azo system using thermal and
        photochemical rates.

        The method extracts cis/trans absorption spectra from the selected electronic
        state, resolves thermal rates (provided or computed), configures the lamp,
        and evaluates the photostationary-state spectrum at each wavelength available
        in that lamp profile.

        Workflow
        --------
        - Finds cis/trans sources and the requested target state in each source.
        - Computes absorption spectra for both isomers in the same wavelength range.
        - Resolves thermal rates:
          - Uses existing 'trans_k'/'cis_k' rates if provided and consistent with 'temp'.
          - Otherwise computes/retrieves them from system results.
        - Builds a `PSS` object and attaches the selected lamp profile.
        - Applies lamp power/intensity scaling with `pw_int`.
        - Computes and stores PSS spectra for all lamp wavelengths.

        Parameters
        ----------
        lamp_name    : str, optional     Name of the lamp profile (e.g. `"COOLLED"`, `"DARK"`). Default is `"default"`.
        target_state : str, optional     State label to use in cis/trans sources (e.g. `"opt"`). Default is `"opt"`.
        temp         : float, optional   Temperature in Kelvin used for thermal-rate handling. Default is 298.15 K.
        phi_EZ       : float, optional   Quantum yield for trans -> cis photoisomerization. Default is 0.3.
        phi_ZE       : float, optional   Quantum yield for cis -> trans photoisomerization. Default is 0.5.
        pw_int       : float, optional   Lamp power/intensity scaling factor passed to `lamp.set_power_intensity()`. Default is 1.0.        
        trans_k      : float, optional   Thermal rate constant for trans -> cis (s^-1). If `None`, it is computed or retrieved from stored results.
        cis_k        : float, optional   Thermal rate constant for cis -> trans (s^-1). If `None`, it is computed or retrieved from stored results.
        lmin         : float, optional   Minimum wavelength (nm) for absorption spectra. Default is 200.
        lmax         : float, optional   Maximum wavelength (nm) for absorption spectra. Default is 1000.
        debug        : int, optional     Debug level. 0: silent, >0: verbose.

        Returns
        -------
        PSS             Initialized PSS object containing spectra, lamp settings, thermal rates, and computed PSS results across lamp wavelengths.
        """

        # Search for cis and trans sources and states.
        found_cis_state, cis_state = self.find_source('cis')[1].find_state(target_state)
        if not found_cis_state:   raise Exception('SYSTEM_AZO.GET_PSS: The target state for the CIS isomer was not found.')
        found_trans_state, trans_state = self.find_source('trans')[1].find_state(target_state)
        if not found_trans_state: raise Exception('SYSTEM_AZO.GET_PSS: The target state for the TRANS isomer was not found.')

        # Extract absorption spectra from two isomers. 
        trans_x, trans_y = trans_state.get_abs_spectrum(lmin=lmin, lmax=lmax, debug=debug) # Need units to compute PSS
        cis_x, cis_y     = cis_state.get_abs_spectrum(lmin=lmin, lmax=lmax, debug=debug)   # Need units to compute PSS
        assert all(trans_x == cis_x)
        if debug > 0: print(f"SYSTEM_AZO.GET_PSS: Spectra of cis and trans computed") 

        # Thermal rates are normally extracted, but can be provided as well
        if trans_k is None: 
            if not 'rate_thermal_trans2cis' in self.results: self.get_trans_halflife_time(temp=temp, target_state=target_state, debug=debug) # Computed if not exists 
            trans_k = self.results['rate_thermal_trans2cis'].value
        else:
            if self.results['rate_thermal_trans2cis'].temperature == temp:
                trans_k = float(trans_k)
            else:
                trans_k = self.get_trans_halflife_time(temp=temp, target_state=target_state, debug=debug).value  # Computed if exists, but wrong temperature

        if cis_k is None: 
            if not 'rate_thermal_cis2trans'   in self.results: self.get_cis_halflife_time(temp=temp, target_state=target_state, debug=debug) # Computed if not exists 
            cis_k   = self.results['rate_thermal_cis2trans'].value
        else:
            if self.results['rate_thermal_cis2trans'].temperature == temp:
                cis_k = float(cis_k)
            else:
                cis_k = self.get_cis_halflife_time(temp=temp, target_state=target_state, debug=debug).value # Computed if exists, but wrong temperature 

        if debug > 0: print(f"SYSTEM_AZO.GET_PSS: Obtained Cis and Trans halflife times")
        if debug > 0: print(f"SYSTEM_AZO.GET_PSS: Creating PSS-Class object")

        # Initializes PSS and computes it for all wl
        self.PSS = PSS(self, trans_x, trans_y, cis_y)
        self.PSS.set_lamp(lamp_name)
        self.PSS.lamp.set_power_intensity(pw_int=pw_int)
        self.PSS.set_thermal_rates(trans_k, cis_k)
        for wl in self.PSS.lamp.wavelengths:
            if debug > 0: print(f"SYSTEM_AZO.GET_PSS: Evaluating PSS for wavelength: {wl}") 
            self.PSS.get_pss_spectrum(wl, phi_EZ=phi_EZ, phi_ZE=phi_ZE, debug=debug) # get_pss_ratio computes the photo_rates already with the same parameters (wl and phi's)
        return self.PSS

    ######
    def __repr__(self):
        to_print = ""
        to_print += '-------------------------------------\n'
        to_print += '   >>> SCOPE System_azo Object >>>   \n'
        to_print += '-------------------------------------\n'
        to_print += f' Name                  = {self.name}\n'
        to_print += f' Dihedral Indices:     = {self.dihedral_indices}\n'
        to_print += System.__repr__(self, indirect=True)
        to_print += '---------------------------------------------\n'
        to_print += '\n'
        return to_print

################
##### PSS Object 
#################
class PSS(object):
    def __init__(self, system, wl_range, trans_spectrum, cis_spectrum): 
        self.system         = system
        self.wl_range       = wl_range
        self.trans_spectrum = trans_spectrum
        self.cis_spectrum   = cis_spectrum
        self.pss_results    = dict()

    def set_thermal_rates(self, trans_t05: float, cis_t05: float, debug: int=0):
        # Thermal rates
        self.rate_thermal_trans2cis = trans_t05
        self.rate_thermal_cis2trans = cis_t05

    def set_lamp(self, name: str="dark", debug: int=0):
        self.lamp = Lamp(name)
        return self.lamp
        
    def set_diameter(self, diameter: float=10):  
        # Relevant to calculate the illuminated area. Provided in mm.
        self.diameter = diameter
    
    def set_area(self, area=None):                          
        # This is to set the illuminated area, relevant to know the amount of light received by the samples
        # Typically, it is a circular beam with diameter=d, so we use the area of a circle
        if area is None: 
            if not hasattr(self,"diameter"): self.set_diameter()
            self.area = np.pi * (self.diameter/2)**2        # diameter in mm 
        else:
            self.area = area                                # value in square milimeters (mm2)
        return self.area

    def get_irradiance(self):
        # Evaluates irradiance (I) as: I = self.lamp.power / self.area
        if not hasattr(self,"area"):        self.set_area()
        self.irradiance = self.lamp.power * self.lamp.power_intensity / self.area 
        return self.irradiance

    def get_photo_rates(self, wl: float, phi_EZ: float=0.3, phi_ZE: float=0.5, debug: int=0):
        # phi_EZ is the quantum yield for the CIS-to-TRANS photoisomerization
        # phi_ZE is the quantum yield for the TRANS-to-CIS photoisomerization
        photon_flux          = self.get_photon_flux_spectrum(wl, debug=debug)

        # Apply proper units to both spectra 
        K = (np.pi * constants.planck_Js * constants.elem_charge) / (constants.epsilon_0 * constants.speed_light * constants.electron_mass)
        rate_photo_trans2cis = phi_EZ * np.trapezoid((self.trans_spectrum * K) * photon_flux, self.wl_range)
        rate_photo_cis2trans = phi_ZE * np.trapezoid((self.cis_spectrum * K) * photon_flux, self.wl_range) 
        return rate_photo_trans2cis, rate_photo_cis2trans 

    def get_pss_ratio(self, wl: float, phi_EZ: float=0.3, phi_ZE: float=0.5, overwrite: bool=True, debug: int=0):
        # Returns relative population of trans isomer
        if not hasattr(self,"rate_thermal_trans2cis"): raise ValueError("PSS.GET_RATIO: Thermal reaction constants missing, set them first with self.set_thermal_rates() before proceeding") 
        if wl not in self.pss_results.keys() or overwrite: 
            # Computes photo rates at this wl and quantum yields (phi_EZ, phi_ZE)
            rate_photo_trans2cis, rate_photo_cis2trans = self.get_photo_rates(wl, phi_EZ, phi_ZE, debug=debug) 
            pss = (self.rate_thermal_cis2trans + rate_photo_cis2trans) / (rate_photo_trans2cis + rate_photo_cis2trans + self.rate_thermal_trans2cis + self.rate_thermal_cis2trans)
            assert pss >= 0.0 and pss <= 1.0

            # Creates a dictionary with results
            results = {
                "wavelength": wl, 
                "pss_trans_ratio": pss,
                "rate_photo_trans2cis": rate_photo_trans2cis, 
                "rate_photo_cis2trans": rate_photo_cis2trans,
                "rate_thermal_trans2cis": self.rate_thermal_trans2cis,
                "rate_thermal_cis2trans": self.rate_thermal_cis2trans, 
                "QY_trans2cis": phi_EZ,
                "QY_cis2trans": phi_ZE,
            } 
            self.pss_results.update({wl: results})
        return pss

    def get_pss_spectrum(self, wl: float, phi_EZ: float=0.3, phi_ZE: float=0.5, overwrite: bool=True, debug: int=0):
        if wl not in self.pss_results.keys() or overwrite: 
            self.get_pss_ratio(wl, phi_EZ, phi_ZE, debug=debug)
        result = self.pss_results[wl]
        pss_spectrum = result['pss_trans_ratio'] * self.trans_spectrum + (1 - result['pss_trans_ratio']) * self.cis_spectrum
        result.update({'pss_spectrum': pss_spectrum})
        return pss_spectrum

    def build_pss_spectrum(self, pss_trans_ratio: float, debug: int=0):
        assert pss_trans_ratio <= 1.0 and pss_trans_ratio >= 0.0
        pss_spectrum = pss_trans_ratio * self.trans_spectrum + (1 - pss_trans_ratio) * self.cis_spectrum
        return pss_spectrum

    def get_photon_flux_spectrum(self, wl: float, debug: int=0):
        '''
        Returns the photon flux spectrum at a given wavelength grid and intensity.
        '''
        from scope.operations.dicts_and_lists import where_in_array
        from scope.operations.vecs_and_mats   import gaussian

        if not hasattr(self,"area"):        self.set_area()
        if not hasattr(self,"lamp"):        self.set_lamp(debug=debug)
        if not hasattr(self.lamp, "power"): self.lamp.set_power()
        if not wl in self.lamp.wavelengths: raise ValueError(f"PSS.GET_PHOTON_FLUX: requested {wl=} is not available in lamp") 

        idx      = where_in_array(self.lamp.wavelengths, wl)[0]
        fwhm     = self.lamp.fwhm[idx]
        # Power Intensity Controls the intensity of the power source. In practice, it acts as a multiplier of self.power
        power    = self.lamp.power[idx] * self.lamp.power_intensity 
        sigma    = fwhm / (2 * np.sqrt(2 * np.log(2)))     # Conversion from FWHM to sigma 
        profile  = gaussian(self.wl_range, wl, sigma=sigma)
        ## Find Intensity along the wavelength space. 
        ## Irradiances are also computed as self.get_irradiance(), but in discrete points, not the spectrum as we do here.  
        irr_spec = power * profile / self.area   # in W/m2/nm

        ## Convert wavelength to meters
        wl_meters = self.wl_range * 1e-9
        # Compute photonic energy as E = h*c/lambda
        photonic_energy = constants.planck_Js * constants.speed_light / wl_meters   # in J   
        return irr_spec / photonic_energy                               # phi: photons * m-2 s-1 nm-1

    #############################
    ## Plots and Visualization ##
    #############################
    def plot_spectra(self, typ: str='all'):
        import matplotlib.pyplot as plt

        for wl, result in self.pss_results.items():
            color = wavelength_to_rgb(float(wl))
            pss_spectrum = result['pss_spectrum']
            pss_ratio = result['pss_trans_ratio']

            # Plots all spectra, irrespectively of the pss_trans_ratio
            plot = False
            if typ == 'all':
                plot = True
            elif typ == 'mixed':  
                if pss_ratio >= 0.01 and pss_ratio <= 0.99:
                    plot = True 

            if plot == True: 
                plt.plot(self.wl_range,pss_spectrum,color=color,label=f"{wl} nm pss: {100 * pss_ratio:.1f}% E")

        plt.plot(self.wl_range, self.trans_spectrum, color='black', label='trans')
        plt.plot(self.wl_range, self.cis_spectrum, color='black',  linestyle='dashed', label='cis')
        plt.legend()
        plt.show()

    def __repr__(self):
        to_print = ""
        to_print += '------------------------------\n'
        to_print += '   >>> SCOPE PSS Object >>>   \n'
        to_print += '------------------------------\n'
        to_print += f' System             = {self.system.name}\n'
        to_print += f' Lamp Name          = {self.lamp.name}\n'
        to_print += f' Lamp Wavelengths   = {self.lamp.wavelengths} (nm)\n'
        to_print += f' Lamp FWHM          = {self.lamp.fwhm} (nm)\n'
        to_print += f' Lamp Powers        = {self.lamp.power} (mW)\n'
        if hasattr(self,"area"): to_print += f' Illuminated Area   = {self.area:8.6f} (mm2)\n' 
        if hasattr(self,"rate_thermal_trans2cis"): to_print += f' Rate Thermal E->Z  = {self.rate_thermal_trans2cis:8.6e} (s)\n' 
        if hasattr(self,"rate_thermal_cis2trans"): to_print += f' Rate Thermal Z->E  = {self.rate_thermal_cis2trans:8.6e} (s)\n' 
        if len(self.pss_results) > 0:              to_print += f' Num of Results     = {len(self.pss_results)}\n'
        to_print += '\n'
        return to_print

############################################
##### MOLECULE Object Adapted to Azo's #####
############################################
class Molecule_azo(Molecule):
    def __init__(self, labels, coord):
        Molecule.__init__(self, labels, coord)
        self.subtype  = "molecule_azo"

    def set_dihedral_indices(self, dih: list):
        self.dihedral_indices = dih
        return self.dihedral_indices

    ##############
    ## Geometry ##
    ##############
    def get_NN_distance(self):
        if not hasattr(self, "dihedral_indices"): 
            raise ValueError(f"MOLECULE_AZO: Please set the dihedral indices with self.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self.dihedral_indices
        return get_dist(self.coord[at2], self.coord[at3])

    def get_inv_angle_left(self):
        if not hasattr(self, "dihedral_indices"): 
            raise ValueError(f"MOLECULE_AZO: Please set the dihedral indices with self.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self.dihedral_indices
        return get_angle_points(self.coord[at1], self.coord[at2], self.coord[at3]) 

    def get_inv_angle_right(self):
        if not hasattr(self, "dihedral_indices"): 
            raise ValueError(f"MOLECULE_AZO: Please set the dihedral indices with self.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self.dihedral_indices
        return get_angle_points(self.coord[at2], self.coord[at3], self.coord[at4]) 

    def get_main_dihedral(self):
        if not hasattr(self, "dihedral_indices"): 
            raise ValueError(f"MOLECULE_AZO: Please set the dihedral indices with self.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self.dihedral_indices
        return get_dihedral(self.coord[at1], self.coord[at2], self.coord[at3], self.coord[at4]) 

    def get_adj_dihedral_left(self):
        if not hasattr(self, "dihedral_indices"): 
            raise ValueError(f"MOLECULE_AZO: Please set the dihedral indices with self.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self.dihedral_indices
        return get_dihedral(self.coord[at0], self.coord[at1], self.coord[at2], self.coord[at3]) 

    def get_adj_dihedral_right(self):
        if not hasattr(self, "dihedral_indices"): 
            raise ValueError(f"MOLECULE_AZO: Please set the dihedral indices with self.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self.dihedral_indices
        return get_dihedral(self.coord[at2], self.coord[at3], self.coord[at4], self.coord[at5]) 
        
    def get_geometry_summary(self, do_print: bool=True):
        if not hasattr(self, "dihedral_indices"): 
            raise ValueError(f"MOLECULE_AZO: Please set the dihedral indices with self.set_dihedral_indices()") 
        to_print =  '-------------------------------------------------\n'
        to_print += ' Main Geometry Features of the Azo group:        \n'
        to_print += '-------------------------------------------------\n'
        to_print += f" NN distance (dNN):         {np.round(self.get_NN_distance(),4)} \n"
        to_print += f" Left Inv. Angle (CNN):     {np.round(np.degrees(self.get_inv_angle_left()),4)} \n"
        to_print += f" Right Inv. Angle (NNC):    {np.round(np.degrees(self.get_inv_angle_right()),4)} \n"
        to_print += f" Main Dihedral (CNNC):      {np.round(np.degrees(self.get_main_dihedral()),4)} \n"
        to_print += f" Adj Dihedral Left (CCNN):  {np.round(np.degrees(self.get_adj_dihedral_left()),4)} \n"
        to_print += f" Adj Dihedral Right (NNCC): {np.round(np.degrees(self.get_adj_dihedral_right()),4)} \n"
        if do_print: print(to_print); return None
        else:        return to_print

    #########################################
    ### Functions to Interact with States ###
    #########################################
    def set_initial_state(self, name: str='initial', debug: int=0):
        ## Same functions than System, but adapted to create State_azo instances
        """Creates the initial state of the specie, with only the geometry"""
        ini_state = self.add_state(name)
        ini_state.set_geometry(self.labels, self.coord)
        return ini_state

    def add_state(self, name: str, debug: int=0):
        ## Same functions than System, but adapted to create State_azo instances
        if not hasattr(self,"states"): setattr(self,"states",list([]))
        exists, new_state = self.find_state(name)
        if exists:  
            if debug > 0: print(f"MOLECULE_AZO.ADD_STATE. State with same {name=} found, returning it")
            return new_state
        else:
            if debug > 0: print("MOLECULE_AZO.ADD_STATE. Creating new state, returning it")
            new_state = State_azo(self, name, debug=debug)  ## Here's the difference with respect to vanilla State class
            self.states.append(new_state)
        return new_state

    def __repr__(self):
        to_print = ""
        to_print += f'---------- SCOPE Molecule_azo Object ------------\n'
        to_print += Molecule.__repr__(self, indirect=True)
        to_print += f"{self.get_geometry_summary(do_print=False)}" 
        return to_print
  
#########################################
##### STATE Object Adapted to Azo's #####
#########################################
class State_azo(State):
    def __init__(self, _source: object, name: str, debug: int=0):
        State.__init__(self, _source, name, debug=debug)
        self.subtype  = "state_azo" 

    ##############
    ## Geometry ##  # Same as in Molecule_azo, but adapted to read dihedral indices from source. Coordinates are read from self
    ##############
    def get_NN_distance(self):
        if not hasattr(self._source, "dihedral_indices"): 
            raise ValueError(f"STATE_AZO: Please set the dihedral indices of the source with source.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self._source.dihedral_indices
        return get_dist(self.coord[at2], self.coord[at3])

    def get_inv_angle_left(self):
        if not hasattr(self._source, "dihedral_indices"): 
            raise ValueError(f"STATE_AZO: Please set the dihedral indices of the source with source.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self._source.dihedral_indices
        return get_angle_points(self.coord[at1], self.coord[at2], self.coord[at3]) 

    def get_inv_angle_right(self):
        if not hasattr(self._source, "dihedral_indices"): 
            raise ValueError(f"STATE_AZO: Please set the dihedral indices of the source with source.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self._source.dihedral_indices
        return get_angle_points(self.coord[at2], self.coord[at3], self.coord[at4]) 

    def get_main_dihedral(self):
        if not hasattr(self._source, "dihedral_indices"): 
            raise ValueError(f"STATE_AZO: Please set the dihedral indices of the source with source.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self._source.dihedral_indices
        return get_dihedral(self.coord[at1], self.coord[at2], self.coord[at3], self.coord[at4]) 

    def get_adj_dihedral_left(self):
        if not hasattr(self._source, "dihedral_indices"): 
            raise ValueError(f"STATE_AZO: Please set the dihedral indices of the source with source.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self._source.dihedral_indices
        return get_dihedral(self.coord[at0], self.coord[at1], self.coord[at2], self.coord[at3]) 

    def get_adj_dihedral_right(self):
        if not hasattr(self._source, "dihedral_indices"): 
            raise ValueError(f"STATE_AZO: Please set the dihedral indices of the source with source.set_dihedral_indices()") 
        at0, at1, at2, at3, at4, at5 = self._source.dihedral_indices
        return get_dihedral(self.coord[at2], self.coord[at3], self.coord[at4], self.coord[at5]) 

    def get_geometry_summary(self, do_print: bool=True):
        if not hasattr(self._source, "dihedral_indices"): 
            raise ValueError(f"STATE_AZO: Please set the dihedral indices of the source with source.set_dihedral_indices()") 
        to_print =  '-------------------------------------------------\n'
        to_print += ' Main Geometry Features of the Azo group:        \n'
        to_print += '-------------------------------------------------\n'
        to_print += f" NN distance (dNN):         {np.round(self.get_NN_distance(),4)} \n"
        to_print += f" Left Inv. Angle (CNN):     {np.round(np.degrees(self.get_inv_angle_left()),4)} \n"
        to_print += f" Right Inv. Angle (NNC):    {np.round(np.degrees(self.get_inv_angle_right()),4)} \n"
        to_print += f" Main Dihedral (CNNC):      {np.round(np.degrees(self.get_main_dihedral()),4)} \n"
        to_print += f" Adj Dihedral Left (CCNN):  {np.round(np.degrees(self.get_adj_dihedral_left()),4)} \n"
        to_print += f" Adj Dihedral Right (NNCC): {np.round(np.degrees(self.get_adj_dihedral_right()),4)} \n"
        if do_print: print(to_print); return None
        else:        return to_print

    ################################
    ## Thermal Properties for Azo ##
    ################################
    def get_gtot_eff(self, temp: float=298.15, p_sh: float=0.0002, debug: int=0):
        from math import log as ln
        '''
        Corrects the Gtot of a triplet Molecule_azo object using the Gtot of the parent Molecule_azo object. 
        Correction is done considering the increase of energy due to surface hopping between the singlet and triplet PESs. 

        Parameters
        ----------
        temp : float, optional         The temperature in Kelvin. The default is 298.15 K.
        overwrite : bool, optional     Whether to overwrite the existing Gtot_corr value. The default is False.
        p_sh : float, optional         The probability of surface hopping. The default is 0.0002.
        debug : int, optional          The debug level. The default is 0
        '''
        # Searches Gtot_eff entry,  in case it already exists
        if 'Gtot_eff' in self.results:
            gtot_eff_col = self.results['Gtot_eff']
            if isinstance(gtot_eff_col, Collection):
                gtot_eff_data = gtot_eff_col.find_value_with_property("temperature", temp)
                if gtot_eff_data is not None: return gtot_eff_data
            else:
                gtot_eff_col = Collection("Gtot_eff", "temperature")
        else:
            gtot_eff_col = Collection("Gtot_eff", "temperature")

        if not 'Gtot' in self.results:
            self.get_thermal_data(temp=temp)
        # Searches specific Gtot(T) entry
        gtot = self.results['Gtot'].find_value_with_property("temperature", temp)
        if gtot is None: 
            self.get_thermal_data(temp=temp)
        # Now it should exists for sure
        gtot = self.results['Gtot'].find_value_with_property("temperature", temp)
        if gtot is None: 
            raise Exception (f'STATE_AZO.GET_Gtot_EFF: Gtot for State {self.name} at {temp=} not found.')
        # Gets Value (energy) is a Data class
        gtot = self.results['Gtot'].find_value_with_property("temperature", temp).convert_to_units("au")

        if self._source.spin == 2:
            # Compute the Penalty, based on Surface Hopping Probability (p_sh)
            dx = - constants.R_J * temp * ln(p_sh)     # J/mol
            dx /= 1000                                 # kJ/mol
            if debug>0: print(f'STATE_AZO.GET_Gtot_EFF: Adding penalty of {dx*constants.kJmol2kcal:.2f} kcal/mol.')
            dx *= constants.kJmol2har
        else: 
            dx = 0

        # Create a Collection, with a single data entry, that of Gtot_eff at the requested temperature
#        gtot_eff_col  = Collection("Gtot_eff", "temperature")
        gtot_eff_data = Data("Gtot_eff", float(gtot.value + dx), "au", "state_azo.get_gtot_eff()") 
        gtot_eff_data.add_property("temperature", temp, overwrite=True)
        gtot_eff_col.add_data(gtot_eff_data)

        # Add the Gtot_eff Collection as a result
        self.add_result(gtot_eff_col, overwrite=True)

        if debug > 0: print (f'STATE_AZO.GET_Gtot_EFF: Corrected Gtot. Data stored in state: {gtot_eff_data.value:12.8f}')
        if debug > 0: print (f'STATE_AZO.GET_Gtot_EFF: Before correction: {gtot.value:12.8f} au. After correction: {gtot_eff_data.value:12.8f}')
        return gtot_eff_data

    ######
    def __repr__(self):
        to_print = ""
        to_print += f'---------- SCOPE State_azo Object ------------\n'
        to_print += State.__repr__(self, indirect=True)
        to_print += f"{self.get_geometry_summary(do_print=False)}" 
        return to_print

################
## Lamp Class ##
################
class Lamp:
    from pathlib import Path
    _LAMP_DATA_DIR = Path(__file__).resolve().parent / "lamp_data" ## Path of where lamp-data files are stored
    def __init__(self, name: str, debug: int=0):
        self.name   = str(name).lower()
        self.data   = self._load_lamp_data(name)
        self.read_wl(debug=debug)
        self.read_fwhm(debug=debug)
        self.read_power(debug=debug)
        self.read_power_intensity(debug=debug)

    ####################
    ## Read from JSON ##
    ####################
    def read_wl(self, debug: int=0):
        if "fwhm_data_nm" in self.data: self.wavelengths   = np.array(list(self._coerce_wavelength_dict(self.data.get("fwhm_data_nm")).keys()))
        else:                           raise ValueError(f"LAMP: No 'fwhm_data_nm' found for lamp '{self.name}'. I can't read the wavelengths")
        return self.wavelengths

    def read_fwhm(self, debug: int=0):
        ## Loads data from JSON.
        if "fwhm_data_nm" in self.data:
            self.fwhm          = np.array(list(self._coerce_wavelength_dict(self.data.get("fwhm_data_nm")).values()))
        elif "default_fwhm_nm" in self.data: 
            self.fwhm          = np.full(len(self.wavelengths), float(self.data.get("default_fwhm_nm"))) # in nm
        else:
            raise ValueError(f"LAMP: No 'fwhm_data_nm' nor default value found for lamp '{self.name}'")
        return self.fwhm

    def read_power(self, debug: int=0):
        if "power_data_mW" in self.data:
            self.power         = np.array(list(self._coerce_wavelength_dict(self.data.get("power_data_mW")).values()))
        elif "default_power_mW" in self.data: 
            self.power         = np.full(len(self.wavelengths), float(self.data.get("default_power_mW"))) # in mW
        else:
            raise ValueError(f"LAMP: Neither 'power_data_mW' nor default value found for lamp '{self.name}'")
        return self.power

    def read_power_intensity(self, debug: int=0):
        if "default_power_intensity" in self.data: 
            self.power_intensity = float(self.data.get("default_power_intensity")) # in mW
            assert self.power_intensity >= 0.0 and self.power_intensity <= 1.0 

    def _load_lamp_data(self, name: str) -> dict:
        from scope.read_write import load_json
        lamp_name = str(name).strip().lower()
        file_path = self._LAMP_DATA_DIR / f"{lamp_name}.json"
        if not file_path.exists():
            available = ", ".join(self.available_lamps())
            raise FileNotFoundError(f"LAMP: No data file found for lamp '{name}'.\n Expected '{file_path.name}' file in '{self._LAMP_DATA_DIR}'. \n Available lamps are: {available}")
        return load_json(file_path)

    def available_lamps(self):
        if not self._LAMP_DATA_DIR.exists():
            return []
        return sorted(path.stem for path in self._LAMP_DATA_DIR.glob("*.json"))

    ###########
    ## Other ##
    ###########
    def set_power_intensity(self, pw_int: float, debug: int=0):
        # Controls the intensity of the power source. In practice, it acts as a multiplier of self.power
        assert pw_int >= 0.0 and pw_int <= 1.0 
        self.power_intensity = pw_int
        return self.power_intensity

    @staticmethod
    def _coerce_wavelength_dict(raw_data):
        if raw_data is None: return {}
        converted = {}
        for key, value in raw_data.items():
            converted[int(round(float(key)))] = float(value)
        return converted

    ##########
    ## repr ##
    ##########
    def __repr__(self):
        to_print = f'--------- Scope Lamp Object -----------\n'
        to_print += f' Name             = {self.name}\n'
        to_print += f' Wavelengths      = {self.wavelengths} nm\n'
        to_print += f' FWHM             = {self.fwhm} nm\n'
        to_print += f' Power            = {self.power} mW\n'
        to_print += f' Power Intensity  = {self.power_intensity}\n'
        to_print += f'---------------------------------------\n'
        return to_print
