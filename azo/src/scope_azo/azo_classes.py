import numpy as np
from    math                                import log as ln
import  scope.constants                     as Constants
from    scope.classes_data                  import Data
from    scope.classes_system                import System
from    scope.classes_state                 import State, find_state
from    scope.classes_specie                import Molecule
from    scope.elementdata                   import ElementData
from    scope.geometry                      import centercoords, set_dihedral, get_dihedral, solve_dihedral, set_angle, get_angle
from    scope.connectivity                  import get_adjmatrix, split_species
from    scope.operations.dicts_and_lists    import extract_from_list
from    scope_azo.azo_functions             import get_3D, compute_t, build_spectrum, get_photon_flux_spectrum, build_pss_spectrum

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
        Creates a set of TS for a System_azo. Users can select which TS to create from a list of options, 
        it must be in ts_list (['TSrot', 'TSinv', 'triplet']). 

        Once the TS are created, they are added to the System_azo object as sources.

        -------
        TSrot
        -------
        TSrot is created by rotating the trans isomer by +/-90 degrees around the azo dihedral. 
        If 'triplet' is selected on ts_list, rotation TS in triplet state are created.

        Nomenclature: 
            - TSrot_A_S: Rotation TS in singlet state. Dihedral angle is set as +90º. Total spin is set to 0.
            - TSrot_A_T: Rotation TS in triplet state. Dihedral angle is set as +90º. Total spin is set as 2.
            - TSrot_B_S: Rotation TS in singlet state. Dihedral angle is set as -90º. Total spin is set as 0.
            - TSrot_B_T: Rotation TS in triplet state. Dihedral angle is set as -90º. Total spin is set as 2.

        -------
        TSinv
        -------
        TSinv is created by inverting the trans isomer around the azo dihedral (setting N=N-ring angle to 180º).
        By default, total spin is set as 0. It can be changed using the function for Molecule objects as Molecule.set_total_spin(value).
        4 versions are created:
            - tsinv_l: Inversion TS involving inversion of left ring (at0 - at1 = at2).
            - tsinv_r: Inversion TS involving inversion of right ring (at1 = at2 - at3).

        They are added as sources of the System_azo. They can be accessed by using System_azo.find_source('TSrot_A_T') or system_azo.find_source('TSinv_R_a'). 
        
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
                angle_deg = np.degrees(get_angle(coord[at1]-coord[at2], coord[at3]-coord[at2]))
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
                angle_deg = np.degrees(get_angle(coord[at3]-coord[at2], coord[at4]-coord[at2]))
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
    def get_mets(self, temp: float=298, target_state: str='opt', skip_triplets: bool=False, debug: int=0):

        ts_names    = []
        ts_energies = []
        for sou in self.sources:
            found, state = sou.find_state(target_state)
            if not found: 
                if debug > 0: print(f'SYSTEM_AZO.GET_METS: Target state not found in source={sou.name}, spin={sou.spin}.')
                continue
            if not hasattr(state,'is_ts'):
                if debug>0: print(f'SYSTEM_AZO.GET_METS: The target state in source {sou.name} is not a TS.')
                continue
            if debug > 0: print(f'SYSTEM_AZO.GET_METS: Found TS: {sou.name}, spin: {sou.spin}, is_ts: {state.is_ts}')

            if state.is_ts:
                # Collects Energy of TSs
                if sou.spin == 2:     # Use corrected Gtot for triplets
                    if not skip_triplets:
                        # Searches Gtot entry
                        if not 'Gtot_eff' in state.results.keys(): 
                            print(f'SYSTEM_AZO.GET_METS: [WARNING] Gtot not found for State={state.name} of source={sou.name}. Computing')
                            state.correct_tripletG(temp=temp, debug=debug).convert_to_units("au").value  
                        # Searches specific Gtot(T) entry
                        energy = state.results['Gtot_eff'].value
                        if energy is None: 
                            state.correct_tripletG(temp=temp, debug=debug)
                        # Now it should exists for sure
                        energy = state.results['Gtot_eff'].value
                        if energy is None: 
                            raise Exception (f'SYSTEM_AZO.GET_METS: Could not fint Gtot_eff({temp}) for {sou.name}.')
                        # Gets Value (energy) is a Data class
                        energy = state.results['Gtot_eff'].convert_to_units("au").value                              

                elif sou.spin == 0:   
                    # Searches Gtot entry
                    if not 'Gtot' in state.results.keys(): 
                        state.get_thermal_data(Trange=range(temp-5,temp+5,1))     
                    # Searches specific Gtot(T) entry
                    energy = state.results['Gtot'].find_value_with_property("temperature", temp)             
                    if energy is None: 
                        state.get_thermal_data(Trange=range(temp-5,temp+5,1))     
                    # Now it should exists for sure
                    energy = state.results['Gtot'].find_value_with_property("temperature", temp)             
                    if energy is None: 
                        raise Exception (f'SYSTEM_AZO.GET_METS: Could not fint Gtot({temp}) for {sou.name}.')
                    # Gets Value (energy) is a Data class
                    energy = state.results['Gtot'].find_value_with_property("temperature", temp).convert_to_units("au").value                              

                ts_names.append(sou.name)
                ts_energies.append(energy)
        if debug == 1: print(f'SYSTEM_AZO.GET_METS: Collected {len(ts_energies)} TSs for {self.name}')
        if debug >= 2: 
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

    def get_thermal_stability(self, target_state: str, debug: int=0):
        # Requires having the energies of Cis and Trans
        found_cis,   cis_iso   = self.find_conformer("cis")
        found_trans, trans_iso = self.find_conformer("trans")
        if found_cis: 
            found_state_cis, state_cis = find_state(cis_iso, target_state)
            if found_state_cis:
                if "energy" in state_cis.results: cis_E = state_cis.results["gs_energy"]
        if found_trans:
            found_state_trans, state_trans = find_state(trans_iso, target_state)
            if found_state_trans:
                if "energy" in state_trans.results: trans_E = state_trans.results["gs_energy"]
        if found_cis and found_trans:
            self.dE = np.round((trans_E - cis_E) * Constants.har2kJmol/Constants.kcal2kJmol,4)
            print(self.name, self.dE, "kJ/mol. Negative means trans is more stable")
        return self.dE

    def get_trans_halflife_time(self, target_state: str='opt', skip_triplets: bool=False, temp: float=298, debug: int=0):
        '''
        
        '''
        found, source = self.find_source('trans')
        if not found:       raise Exception('SYSTEM_AZO.GET_TRANS_HALFLIFE_TIME: Trans source not found.')
        found_state, ground_state = source.find_state(target_state)
        if not found_state: raise Exception(f'SYSTEM_AZO.GET_TRANS_HALFLIFE_TIME: Target state: {target_state} not found.')

        # Getting the energy of the Trans isomer
        if not 'Gtot' in ground_state.results.keys(): 
            ground_state.get_thermal_data(Trange=range(temp-5,temp+5,1)) ## If not available, computes it
        # Searches specific Gtot(T) entry
        Gtot_T = ground_state.results['Gtot'].find_value_with_property("temperature", temp)                     ## Gtot is a collection. We get the value for temp=T 
        if Gtot_T is None: 
            ground_state.get_thermal_data(Trange=range(temp-5,temp+5,1)) ## If not available, computes it
        # Now it should exists for sure
        Gtot_T = ground_state.results['Gtot'].find_value_with_property("temperature", temp)                     ## Gtot is a collection. We get the value for temp=T 
        if Gtot_T is None: 
            raise Exception (f'SYSTEM_AZO.GET_TRANS_HALFLIFE_TIME: Gtot for the Ground State not found.')
        # Gets Value (energy) is a Data class
        gtot_ground = ground_state.results['Gtot'].find_value_with_property("temperature", temp).convert_to_units("au").value

        # Getting the energy of the METS (Minimum Energy Transition State)
        mets = self.get_mets(temp=temp, target_state=target_state, skip_triplets=skip_triplets, debug=debug)
        found, mets_state = mets.find_state(target_state)

        if mets.spin == 0:     # Singlets
            gtot_mets = mets_state.results['Gtot'].find_value_with_property("temperature", temp).convert_to_units("au").value     ## Gtot is a collection. We get the value for temp=T 
        elif mets.spin == 2:   # Triplets 
            gtot_mets = mets_state.results['Gtot_eff'].find_value_with_property("temperature", temp).convert_to_units("au").value ## Gtot is a collection. We get the value for temp=T 
        else:                  # Other 
            raise ValueError(f"SYSTEM_AZO.GET_TRANS_T05: Gtot_eff only implemented for triplet states. Current State has spin={mets_state.spin}")
        
        #Computing halflife time (t05) and kinetic constant (k)
        t05, k_th = compute_t(gtot_ground, gtot_mets, temp)
        if debug > 0: print(f'SYSTEM_AZO.GET_TRANS_T05: Obtained t05: {t05}  k_th: {k_th}') 
        dG_from_trans = (float(gtot_mets)- float(gtot_ground)) * Constants.har2kJmol * Constants.kJmol2kcal # in Kcal/mol
        if debug > 0: print(f'SYSTEM_AZO.GET_TRANS_T05: dG_from_trans = {dG_from_trans}') 

        units_time = 's'
        units_kcal = 'kcal/mol'

        self.add_result(Data("dG_from_trans",dG_from_trans,units_kcal,"system_azo.get_trans_halflife_time()"), overwrite=True)
        self.add_result(Data("trans_halflife_time",t05,units_time,"system_azo.get_trans_halflife_time()"), overwrite=True)
        return self.results['trans_halflife_time']
        
    ######
    def get_cis_halflife_time(self, target_state: str='opt', skip_triplets: bool=False, temp: float=298, debug: int=0):
        '''
        
        '''
        found, source = self.find_source('cis')
        if not found:       raise Exception('SYSTEM_AZO.GET_CIS_HALFLIFE_TIME: cis source not found.')
        found_state, ground_state = source.find_state(target_state)
        if not found_state: raise Exception(f'SYSTEM_AZO.GET_CIS_HALFLIFE_TIME: Target state: {target_state} not found.')

        # Getting the energy of the cis isomer
        if not 'Gtot' in ground_state.results.keys(): 
            ground_state.get_thermal_data(Trange=range(temp-5,temp+5,1)) ## If not available, computes it
        # Searches specific Gtot(T) entry
        Gtot_T = ground_state.results['Gtot'].find_value_with_property("temperature", temp)                     ## Gtot is a collection. We get the value for temp=T 
        if Gtot_T is None: 
            ground_state.get_thermal_data(Trange=range(temp-5,temp+5,1)) ## If not available, computes it
        # Now it should exists for sure
        Gtot_T = ground_state.results['Gtot'].find_value_with_property("temperature", temp)                     ## Gtot is a collection. We get the value for temp=T 
        if Gtot_T is None: 
            raise Exception (f'SYSTEM_AZO.GET_CIS_HALFLIFE_TIME: Gtot for the Ground State not found.')
        # Gets Value (energy) is a Data class
        gtot_ground = ground_state.results['Gtot'].find_value_with_property("temperature", temp).convert_to_units("au").value

        # Getting the energy of the METS (Minimum Energy Transition State)
        mets = self.get_mets(temp=temp, target_state=target_state, skip_triplets=skip_triplets, debug=debug)
        found, mets_state = mets.find_state(target_state)

        if mets.spin == 0:     # Singlets
            gtot_mets = mets_state.results['Gtot'].find_value_with_property("temperature", temp).convert_to_units("au").value     ## Gtot is a collection. We get the value for temp=T 
        elif mets.spin == 2:   # Triplets 
            gtot_mets = mets_state.results['Gtot_eff'].find_value_with_property("temperature", temp).convert_to_units("au").value ## Gtot is a collection. We get the value for temp=T 
        else:                  # Other 
            raise ValueError(f"SYSTEM_AZO.GET_CIS_T05: Gtot_eff only implemented for triplet states. Current State has spin={mets_state.spin}")
        
        #Computing halflife time (t05) and kinetic constant (k)
        t05, k_th = compute_t(gtot_ground, gtot_mets, temp)
        if debug > 0: print(f'SYSTEM_AZO.GET_CIS_T05: Obtained t05: {t05}  k_th: {k_th}') 
        dG_from_cis = (float(gtot_mets)- float(gtot_ground)) * Constants.har2kJmol * Constants.kJmol2kcal # in Kcal/mol
        if debug > 0: print(f'SYSTEM_AZO.GET_CIS_T05: dG_from_cis = {dG_from_cis}') 

        units_time = 's'
        units_kcal = 'kcal/mol'

        self.add_result(Data("dG_from_cis",dG_from_cis,units_kcal,"system_azo.get_cis_halflife_time()"), overwrite=True)
        self.add_result(Data("cis_halflife_time",t05,units_time,"system_azo.get_cis_halflife_time()"), overwrite=True)
        return self.results['cis_halflife_time']

    ########################
    ## Optical Properties ##
    ########################
    def get_PSS(self, target_state: str, lamp : "Lamp", phi_EZ = 0.3, phi_ZE = 0.5, t_EZ=None, t_ZE=None, debug=0):
        """
        Function to calculate the photostationary state (PSS) for a given System_azo, based on the photochemical and thermal rates.
        Returns the PSS value, which is the fraction of the Trans isomer at the PSS.
        """

        # Search for cis and trans sources and states.
        found_cis_state, cis_state = self.find_source('cis')[1].find_state(target_state)
        if not found_cis_state:   raise Exception('SYSTEM_AZO.GET_PSS: The target state for the CIS isomer was not found.')
        found_trans_state, trans_state = self.find_source('trans')[1].find_state(target_state)
        if not found_trans_state: raise Exception('SYSTEM_AZO.GET_PSS: The target state for the TRANS isomer was not found.')

        trans_abs_spectrum = trans_state.get_abs_spectrum(normalize=False, units=True, debug=debug) # Need units to compute PSS
        cis_abs_spectrum   = cis_state.get_abs_spectrum(normalize=False, units=True, debug=debug)   # Need units to compute PSS

        # Photon flux from lamp
        photon_flux = get_photon_flux_spectrum(lamp.wavelength, lamp.fwhm, lambda_grid, Itot=lamp.irradiance)

        # Gets Half-life times (in seconds). If not provided, it computes them using the corresponding functions.
        if t_EZ is None:
            if not 'trans_halflife_time' in self.results.keys():
                print('SYSTEM_AZO.get_PSS: No halflife time found for trans isomer. Computing...')
                self.get_trans_halflife_time(skip_triplets=False, overwrite=False)               
            assert self.results['trans_halflife_time'].units == 'seconds'
            t_EZ = self.results['trans_halflife_time'].value
        else: 
            assert t_EZ > 0, "SYSTEM_AZO.get_PSS: t_EZ must be a positive value in seconds."
        
        if t_ZE is None:
            if not 'cis_halflife_time' in self.results.keys():
                print('SYSTEM_AZO.get_PSS: No halflife time found for cis isomer. Computing...')
                self.get_cis_halflife_time(skip_triplets=False, overwrite=False)                  
            assert self.results['cis_halflife_time'].units == 'seconds'
            t_ZE = self.results['cis_halflife_time'].value
        else: 
            assert t_ZE > 0, "SYSTEM_AZO.get_PSS: t_ZE must be a positive value in seconds."

        # Photochemical rates
        k_ph_EZ = phi_EZ * np.trapezoid(trans_abs_spectrum * photon_flux, lambda_grid)      ## MANEL: ojo, ara trans_abs_spectrum es un array de 2 columnes, (x, y). Adapta-ho com calgui
        k_ph_ZE = phi_ZE * np.trapezoid(cis_abs_spectrum   * photon_flux, lambda_grid)      ##      per altra banda, lambda grid no se d'on surt. Pot ser que sigui trans_abs_spectrum[:,0]? 
        # Thermal rates
        k_th_EZ = np.log(2) / t_EZ
        k_th_ZE = np.log(2) / t_ZE
        # Steady-state populations
        self.PSS = (k_th_ZE + k_ph_ZE) / (k_ph_EZ + k_ph_ZE + k_th_EZ + k_th_ZE)        # Fraction of Trans isomer aka N_E. 
        return self.PSS 

    ######
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

        sigma_pss = build_pss_spectrum(pss_B, sigma_trans, sigma_cis, debug=debug) 
        frac_list.append(pss_B)     
        pss_list.append(sigma_pss)  

        # 2. Iterate through irradiations
        for irr_wl in wl_list:
            lamp = Lamp(name='COOLLED', wavelength=irr_wl, shift_nm=shift_nm)
            _, _, _, pss_B = self.get_PSS(lamp, phi_EZ=phi_EZ, phi_ZE=phi_ZE, t_EZ=t_EZ, t_ZE=t_ZE, debug=debug) 
            sigma_pss = build_pss_spectrum(pss_B, sigma_trans, sigma_cis, debug=debug)
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

    ######
    def __repr__(self):
        to_print = ""
        to_print += f'------------- SCOPE System_azo Object -------\n'
        to_print += f' Name                  = {self.name}\n'
        to_print += f' Dihedral Indices:     = {self.dihedral_indices}\n'
        to_print += System.__repr__(self, indirect=True)
        if hasattr(self,"dE"): to_print += f' Thermal Stability     = {self.dE} kJ/mol (- means trans is more stable)\n'
        to_print += '---------------------------------------------\n'
        to_print += '\n'
        return to_print


#########################################
##### STATE Object Adapted to Azo's #####
#########################################
class State_azo(State):
    def __init__(self, _source: object, name: str, debug: int=0):
        State.__init__(self, _source, name, debug=debug)
        self.subtype  = "state_azo" 

    ######
    def correct_tripletG(self, temp: float=298.15, p_sh: float=0.0002, debug: int=0):
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
        # Searches Gtot entry
        if not 'Gtot' in self.results:
            self.get_thermal_data(Trange=range(temp-5,temp+5,1))                                         
        # Searches specific Gtot(T) entry
        gtot = self.results['Gtot'].find_value_with_property("temperature", temp)
        if gtot is None: 
            self.get_thermal_data(Trange=range(temp-5,temp+5,1))                                         
        # Now it should exists for sure
        gtot = self.results['Gtot'].find_value_with_property("temperature", temp)
        if gtot is None: 
            raise Exception (f'STATE_AZO.CORRECT_TRIPLETG: Gtot for State {self.name} at {temp=} not found.')
        # Gets Value (energy) is a Data class
        gtot = self.results['Gtot'].find_value_with_property("temperature", temp).convert_to_units("au").value   

        dx = - Constants.R_J * temp * ln(p_sh)     # J/mol
        dx /= 1000                                 # kJ/mol
        if debug>0: print(f'STATE_AZO.CORRECT_TRIPLETG: Adding penalty of {dx*Constants.kJmol2kcal:.2f} kcal/mol.')
        dx *= Constants.kJmol2har
        gtot_eff = Data("Gtot_eff", float(gtot + dx), "au", "state_azo.correct_tripletg()")
        gtot_eff.add_property("temperature", temp, overwrite=True)
        self.add_result(gtot_eff, overwrite=True)
        if debug > 0: print (f'STATE_AZO.CORRECT_TRIPLETG: Corrected Gtot. Data stored in state: {gtot_eff}')
        if debug > 0: print (f'STATE_AZO.CORRECT_TRIPLETG: Before correction: {gtot} au. After correction: {gtot_eff}')
        return gtot_eff

    ######
    def get_abs_spectrum(self, normalize: bool = False, units: bool = False, lmin: float=200, lmax: float=1000, debug: int=0):
        # Check if TDDFT data exists.
        if not hasattr(self, 'exc_states'): raise ValueError('AZO.GET_ABS_SPECTRUM: [WARNING] No TDDFT data found in this state')

        # Collects Values
        energies = [es.energy for es in self.exc_states]
        fosc     = [es.fosc for es in self.exc_states]
        erange   = np.linspace(lmin, lmax, lmax-lmin)

        # Builds the spectrum from discrete values, using Gaussian broadening
        self.abs_spectrum = build_spectrum(erange, energies, fosc, normalize=normalize, units=units)
        return self.abs_spectrum


############################################
##### MOLECULE Object Adapted to Azo's #####
############################################
class Molecule_azo(Molecule):
    def __init__(self, labels, coord):
        Molecule.__init__(self, labels, coord)
        self.subtype  = "molecule_azo"

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

    ######        
    ## MANEL: es fa servir aquesta funcio?

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
