##############################
##### General CELL Class #####
##############################

import numpy as np
from scope.connectivity     import *
from scope.geometry         import cellparam_2_cellvec, cellvec_2_cellparam, cart2frac, get_unit_cell_volume
from scope.elementdata      import ElementData
from scope.classes_specie   import Molecule
elemdatabase = ElementData()

##############
#### CELL ####
##############
class Cell(object):
    """
    Represent a periodic crystal cell and its molecular decomposition.

    Attributes:
        object_type (str):              Object category (`"cell"`).
        name (str):                     Cell name.
        labels (list):                  Atomic symbols.
        coord (list):                   Cartesian coordinates.
        cell_vector (list):             Lattice vectors.
        cell_param (list):              Lattice parameters.
        volume (float):                 Unit-cell volume.

    Methods:
        get_molecules():                Split the cell into molecular fragments.
        get_atoms():                    Collect atom objects from molecules.
        get_z():                        Compute the stoichiometric multiplicity.
        set_spin_metals():              Assign spins to metal centers.
        associate_cif():                Link the cell to a CIF object.
    """

    def __init__(self, name: str, labels: list, coord: list, cell_vector: list=None, cell_param: list=None) -> None:
        self.version              = "1.0"
        self.object_type          = "cell"
        self.object_subtype       = "cell"
        self.origin               = "created"
        self.name                 = name
        self.labels               = labels 
        self.coord                = coord
        self.formula              = labels2formula(labels)
        self.natoms               = len(labels)

        ## Gets Cell Parameters, Vectors and Fractional Coordinates
        if   cell_vector is None and cell_param is None:
            raise ValueError("CELL: Either cell_vector or cell_param must be provided to create a cell object")
        elif cell_vector is None and cell_param is not None:
            self.cell_param       = cell_param
            self.cell_vector      = cellparam_2_cellvec(cell_param)
        elif cell_vector is not None and cell_param is None:
            self.cell_vector      = cell_vector
            self.cell_param       = cellvec_2_cellparam(cell_vector)
        else:
            self.cell_vector      = cell_vector
            self.cell_param       = cell_param
        self.frac_coord           = cart2frac(self.coord, self.cell_vector)

        ## Gets the volume
        self.volume               = get_unit_cell_volume(*self.cell_param) 

    #####################
    ## Charge and Spin ##
    #####################
    @property
    def charge(self):
        if not hasattr(self,"atoms"): self.get_atoms()
        if not all(hasattr(at,"charge") for at in self.atoms): return None
        return sum(at.charge for at in self.atoms)

    @property
    def atomic_charges(self):
        if not hasattr(self,"atoms"): self.get_atoms()
        if not all(hasattr(at,"charge") for at in self.atoms): return None
        return [at.charge for at in self.atoms]

    @property
    def spin(self):
        if not hasattr(self,"atoms"): self.get_atoms()
        if not all(hasattr(at,"spin") for at in self.atoms): return None
        return sum(at.spin for at in self.atoms)

    @property
    def atomic_spins(self):
        if not hasattr(self,"atoms"): self.get_atoms()
        if not all(hasattr(at,"spin") for at in self.atoms): return None
        return [at.spin for at in self.atoms]

    @property
    def ismagnetic(self):
        for at_s in self.atomic_spins:
            if at_s != 0: return True
        return False
   
    @property
    def spin_multiplicity(self):
        return int(self.spin + 1) 

    ######################
    ## Other Properties ##
    ######################
    @property
    def Z(self):
        """Return `z` as a convenience alias."""
        return self.z

    @property
    def z(self):
        """Return the cached stoichiometric multiplicity, computing it if needed."""
        if not hasattr(self, "_z"):
            return self.get_z()
        return self._z

    #### 
    def reset_charge(self, debug: int=0):
        """Reset charges for all molecules stored in the cell."""
        if not hasattr(self,"molecules"): return None
        for mol in self.molecules:
            mol.reset_charge()

    def reset_spin(self, debug: int=0):
        """Reset spins for all molecules stored in the cell."""
        if not hasattr(self,"molecules"): return None
        for mol in self.molecules:
            mol.reset_spin()

    ####
    def set_spin_metals(self, spins: list | int, debug: int=0):
        """Assign spin values to metal centers in transition-metal complexes.

        Parameters:
            spins (list | int):          One spin per complex, or one value for all.
            debug (int):                 Verbosity level.

        Returns:
            int | None: Total cell spin, or `None` if no complexes are found.
        """
        ## Function to simplify setting the spin for the Transition Metal Complexes of this cell
        if not hasattr(self,"molecules"): self.get_molecules(debug=debug)
        ncomplex = self.get_ncomplex(debug=debug) 
        if ncomplex == 0:
            print(f"CELL.SET_SPIN_METALS: there are no Transition Metal Complexes in this cell")
            return None
        ## Checks
        if isinstance(spins, int): spins = [spins] * ncomplex  ## If spin is an integer, then it assumes it aplies to all metals
        ## Verbose
        if debug > 0: 
            print(f"CELL.SET_SPIN_METALS: Preparing Spin Configuration for Specie {self.formula}")
            print(f"CELL.SET_SPIN_METALS: Received {spins=}")
        ## Main
        pointer = 0
        for mol in self.molecules:
            if mol.iscomplex: 
                if len(mol.metals) > 1: 
                    print(f"CELL.SET_SPIN_METALS: Molecule {mol.formula} has more than one metal atom.") 
                    print(f"Please set the spin for each molecule separately using mol.set_spin()")
                else: 
                    met = mol.metals[0]
                    if debug > 0: print(f"CELL.SET_SPIN_METALS: Setting spin={spins[pointer]} to metal atom {met.label}")
                    met.set_spin(spins[pointer]); pointer += 1
        return self.spin

    ###########
    ## Other ##
    ###########
    def save(self, filepath):
        """Serialize the cell to disk."""
        save_binary(self, filepath)

    def associate_cif(self, cif: object) -> None:
        """Attach a CIF object to the cell."""
        self.cif       = cif
        return self.cif

    def set_path(self, path: str) -> None:
        """Store a filesystem path for the cell."""
        self.path = path

    def get_ncomplex(self, debug: int=0):
        # Counts transition-metal complexes in `molecules`.
        if debug > 0: print(f"CELL.GET_NCOMPLEX checking fragmentation")
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"Found Fragmented molecules in the geometry of cell: {self.name}"
         
        if debug > 0: print(f"CELL.GET_NCOMPLEX getting molecules")
        if not hasattr(self,"molecules"): self.get_molecules(debug=debug)
        self.ncomplex = 0
        for mol in self.molecules:
            if mol.iscomplex: self.ncomplex += 1
        if debug > 0: print(f"CELL.GET_NCOMPLEX {self.ncomplex} complexes found in cell: {self.name}")
        return self.ncomplex

    def get_z(self, debug: int=0):
        """Compute and cache the unit-cell stoichiometric multiplicity (`Z`).

        Parameters:
            debug (int):                 Verbosity level.

        Returns:
            int: Stoichiometric multiplicity of the cell.
        """
        from scope.operations.vecs_and_mats import gcd_list

        if debug > 0: print(f"self.GET_Z checking fragmentation")
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"self.GET_Z found Fragmented molecules in the geometry of self: {self.name}"

        if not hasattr(self,"molecules"): 
            if debug > 0: print(f"self.GET_Z getting molecules")
            self.get_molecules(debug=debug)

        unique = [] 
        occurrences = []
        for mol in self.molecules:
            found = False
            for uni in unique:
                if mol == uni: found = True
            if not found: 
                unique.append(mol)
                occurrences.append(self.get_occurrence(mol, debug=debug))
        self._z = int(gcd_list(occurrences))
        return self._z 

    #############
    ## Parents ##
    #############
    ## This function mimics the specie-class function with the same name
    def check_parent(self, subtype):
        """Return whether this cell matches the requested parent subtype."""
        if subtype == "cell": return True
        else:                 return False

    ## This function mimics the specie-class function with the same name
    def get_parent(self, subtype):
        """Return this cell when `subtype == "cell"`."""
        if subtype == "cell": return self
        else:                 return None

    ## This function mimics the specie-class function with the same name
    def get_parent_indices(self, subtype: str):
        """Return atom indices for the cell parent subtype."""
        if subtype == "cell": return list(range(0,self.natoms))
        else:                 return None

    #########################
    ## Atoms and Molecules ##
    #########################
    def get_atoms(self, debug: int=0):
        """Collect atom objects from the stored molecules."""
        if not hasattr(self,"molecules"): 
            if debug > 0: print("CELL.GET_ATOMS: retrieving molecules")
            self.get_molecules()
        self.atoms = []
        for mol in self.molecules:
            for at in mol.atoms:
                self.atoms.append(at)
        return self.atoms

    ######
    def set_molecules(self, molecules: list) -> None:
        """Store the molecule list used by the cell."""
        self.molecules = molecules
    
    ######
    def get_molecules(self, overwrite: bool=False, cov_factor: float=1.3, metal_factor: float=1.0, debug: int=0):
        """Build or return the molecules detected in the cell.

        Parameters:
            overwrite (bool):            Recompute even if molecules already exist.
            cov_factor (float):          Covalent radii scaling factor.
            metal_factor (float):        Metal-specific scaling factor.
            debug (int):                 Verbosity level.

        Returns:
            list | None: Molecules found in the cell, or `None`.
        """
        ## Overwrite and Warning
        if not overwrite and hasattr(self,"molecules"): 
            if debug > 0: print(f"CELL.GET_MOLECULES. Molecules already exist and default is overwrite=False")
            return self.molecules

        ## Security
        if not hasattr(self,"labels") or not hasattr(self,"coord"): 
            if debug > 0: print(f"CELL.GET_MOLECULES. Labels or coordinates not found. Returning None")
            return None
        if len(self.labels) == 0 or len(self.coord) == 0:           
            if debug > 0: print(f"CELL.GET_MOLECULES. Empty labels or coordinates. Returning None")
            return None
        if debug > 0: print(f"CELL.GET_MOLECULES passed initial checks")

        ## If fragmented, then it reconstructs the unit cell
        if hasattr(self,"fragmented"): 
            if self.fragmented:
                print("CELL.GET_MOLECULES. Reconstructing unit cell")
                self.reconstruct(debug=debug)
                print("CELL.GET_MOLECULES. Arranging new coordinates from molecules")
                new_coord = self.fix_cell_coord()
                if new_coord is None:
                    print("CELL.GET_MOLECULES. NONE received from FIX_CELL_COORD. Stopping")
                    return None 

        blocklist = split_species(self.labels, self.coord, cov_factor=cov_factor, debug=debug)
        if blocklist is None: 
            return None
        else :
            if debug > 0: print(f"CELL.GET_MOLECULES: found {len(blocklist)} blocks")
            if debug > 0: print(f"CELL.GET_MOLECULES: {blocklist=}")
        
        self.molecules = []
        for b in blocklist:
            if debug > 0: print(f"CELL.GET_MOLECULES: doing block={b}")
            mol_labels      = extract_from_list(b, self.labels, dimension=1)
            mol_coord       = extract_from_list(b, self.coord, dimension=1)
            mol_frac_coord  = extract_from_list(b, self.frac_coord, dimension=1)
            # Creates Molecule Object
            newmolec        = Molecule(mol_labels, mol_coord, mol_frac_coord)
            # For debugging
            newmolec.origin = "cell.get_molecules"
            # Adds cell as parent of the molecule, with indices b
            newmolec.add_parent(self, indices=b)            
            # Store Adjacency Parameters
            newmolec.set_factors(cov_factor, metal_factor)
            # Creates The atom objects with adjacencies
            newmolec.set_atoms(create_adjacencies=True, debug=debug)
            # The split_complex must be below the frac_coord, so they are carried on to the ligands
            if newmolec.iscomplex: 
                if debug > 0: print(f"CELL.GET_MOLECULES: splitting complex")
                newmolec.split_complex(debug=debug)
            # Not needed here, as the reconstruction will take care of it
            self.molecules.append(newmolec)
        return self.molecules

    ######
    def get_occurrence(self, substructure: object, debug: int=0) -> int:
        """
        Count how many times a substructure appears in the cell.

        Parameters:
            substructure (object):       Substructure to search for.
            debug (int):                 Verbosity level.

        Returns:
            int: Number of occurrences found.
        """
        ## Finds how many times a substructure appears in self
        occurrence = 0

        if debug > 0: print(f"CELL.GET_OCCURRENCE checking fragmentation")
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"CELL.GET_OCCURRENCE found fragmented molecules in the geometry of cell: {self.name}"
         
        if debug > 0: print(f"CELL.GET_OCCURRENCE getting molecules")
        if not hasattr(self,"molecules"): self.get_molecules(debug=debug)

        ## Case of Species inside self
        if hasattr(substructure,"type"):
            if substructure.object_type == 'specie':
                for mol in self.molecules:
                    if mol.__eq__(substructure, with_graph=True): occurrence += 1
        return occurrence

    ##################
    ## Connectivity ##
    ##################
    def get_adjmatrix(self, adjust_factor: bool=False, debug: int=0):
        """Compute and cache the cell adjacency matrix."""
        isgood, adjmat, adjnum = get_adjmatrix(self.labels, self.coord, adjust_factor=adjust_factor, debug=debug)
        if isgood:
            self.adjmat = adjmat
            self.adjnum = adjnum
        else:
            self.adjmat = None
            self.adjnum = None
        return self.adjmat, self.adjnum

    ######
    def check_fragmentation(self, reconstruct: bool = False, debug: int=0):
        """Check whether the current molecular decomposition is fragmented.

        Parameters:
            reconstruct (bool):          Attempt reconstruction before returning.
            debug (int):                 Verbosity level.

        Returns:
            bool: `True` if fragmented, `False` otherwise.
        """
        if not hasattr(self,"molecules"): 
            if debug > 0: print("CELL.CHECK_FRAGMENTATION: retrieving molecules")
            self.get_molecules()

        if debug == 1: 
            print("CELL.CHECK_FRAGMENTATION: molecules available. Checking Fragmentation")
        elif debug == 2: 
            print("CELL.CHECK_FRAGMENTATION: molecules available (below). Checking Fragmentation")
            for mol in self.molecules:
                print(mol.formula)

        self.fragmented = False
        # First comparison with ref_molecules
        for mol in self.molecules:
            found = False
            for rmol in self.ref_molecules:
                issame = mol.__eq__(rmol, with_graph=False) # It is very unlikely that the graph is needed in this case. 
                if issame: found = True
                if debug > 0: print(f"CELL.CHECK_FRAGMENTATION: compared {mol.formula} and {rmol.formula} without graph analysis: {issame=}")
            if not found: 
                if debug > 0: print(f"CELL.CHECK_FRAGMENTATION: {mol.formula} not found in ref_molecules")
                self.fragmented = True

        if self.fragmented:
            if debug > 0: print(f"CELL.CHECK_FRAGMENTATION: cell is fragmented")
        else:
            if debug > 0: print(f"CELL.CHECK_FRAGMENTATION: cell is not fragmented")

        # If there are fragments and user wants reconstruction, it tries to reconstruct and checks the new molecules
        if self.fragmented and reconstruct:
            if debug > 0: print("CELL.CHECK_FRAGMENTATION: reconstructing")
            new_molecules = self.reconstruct(debug=debug)
            self.fragmented = False
            for mol in new_molecules:
                found = False
                for rmol in self.ref_molecules:
                    if rmol == mol: found = True
                if not found: self.fragmented = True
            if not self.fragmented: 
                self.molecules = new_molecules
               #self.set_geometry_from_molecules()
        return self.fragmented

    ######
    def fix_cell_coord(self, debug: int=0) -> None:
        """Rebuild cell labels and coordinates from the stored molecules."""
        ## In cell2mol, the cell object does not have the coordinates of the reconstructed cell.
        ## However, the molecule and atom objects are updated (i.e. reconstructed). We use this info to update the cell
        if not hasattr(self,"molecules"): 
            if debug > 0: print(f"CELL.FIX_CELL_COORD: creating molecules")
            self.get_molecules(debug=debug)

        self.labels = []
        self.coord  = []
        indices     = []
        for mol in self.molecules:
            for idx, a in enumerate(mol.atoms):
                self.labels.append(a.label)
                self.coord.append(a.coord)
                if hasattr(a,"parent_index"): indices.append(a.parent_index)
                elif hasattr(a,"index"):      indices.append(a.index)
                else:                         indices.append(idx)

        ## Below is to order the atoms as in the original cell, using the indices stored in the molecule object
        self.labels = [x for _, x in sorted(zip(indices, self.labels), key=lambda pair: pair[0])]
        self.coord  = [x for _, x in sorted(zip(indices, self.coord), key=lambda pair: pair[0])]
        assert len(self.labels) == len(self.coord)
        return self.coord

    ######
    def reconstruct(self, cov_factor: float=None, metal_factor: float=None, debug: int=0):
        """Reconstruct fragmented molecules in the periodic cell.

        Parameters:
            cov_factor (float | None):   Covalent radii scaling factor.
            metal_factor (float | None): Metal-specific scaling factor.
            debug (int):                 Verbosity level.

        Returns:
            list: Reconstructed molecule list.
        """
        from scope.reconstruct import classify_fragments, fragments_reconstruct

        if not hasattr(self,"fragmented"): self.check_fragmentation()
        if not self.fragmented:
            print("CELL.RECONSTRUCT. Cell is not fragmented")
            return self.molecules

        if not hasattr(self,"ref_molecules"): print("CELL.RECONSTRUCT. CELL missing list of reference molecules"); return
        if cov_factor is None:   cov_factor   = self.ref_molecules[0].cov_factor
        if metal_factor is None: metal_factor = self.ref_molecules[0].metal_factor

        ## Get the fragments, which is the molecules of a fragmented cell
        fragments = self.get_molecules(cov_factor=cov_factor, metal_factor=metal_factor, debug=debug)
        if fragments is None: self.error_get_fragments = True; return  
        else:                 self.error_get_fragments = False

        ## Classifies fragments
        molecules, fragments, hydrogens = classify_fragments(fragments, self.ref_molecules, debug=debug)
        if debug > 0: print(f"CELL.RECONSTRUCT: {len(molecules)} {molecules=}")
        if debug > 0: print(f"CELL.RECONSTRUCT: {len(fragments)} {fragments=}")
        if debug > 0: print(f"CELL.RECONSTRUCT: {len(hydrogens)} {hydrogens=}")

        ## Determines if Reconstruction is necessary
        if len(fragments) > 0 or len(hydrogens) > 0: self.is_fragmented = True
        else:                                        self.is_fragmented = False
        
        self.molecules = []
        if not self.is_fragmented: 
            for mol in molecules:
                self.molecules.append(mol)
            return self.molecules     
        else:
            reconstructed_molecules, Warning = fragments_reconstruct(molecules, fragments, hydrogens, self.ref_molecules, self.cell_vector, cov_factor, metal_factor)
            
            if Warning:
                self.is_fragmented = True
                self.error_reconstruction = True 

            else :
                self.is_fragmented = False
                self.error_reconstruction = False 

            ## For consistency, we create the molecules once again, even if mol is already a molecule-class object.
            ## One must follow the same structure as in self.get_molecules()
            if debug > 0: print(f"CELL.RECONSTRUCT: Creating and preparing molecules")
            for idx, mol in enumerate(reconstructed_molecules):
                if debug > 0: print(f"CELL.RECONSTRUCT: Doing molecule {idx} with {mol.formula=}")
                newmolec = Molecule(mol.labels, mol.coord)
                newmolec.origin = "cell.reconstruct"
                newmolec.set_factors(cov_factor, metal_factor)
                newmolec.set_atoms(create_adjacencies=True, debug=debug)
                newmolec.add_parent(self, mol.cell_indices, debug=debug) 
                if debug > 0: print(f"CELL.RECONSTRUCT: Setting fractional coordinates")
                newmolec.set_fractional_coord(mol.frac_coord)
                if newmolec.iscomplex: newmolec.split_complex()
                self.molecules.append(newmolec)         
            return self.molecules

    #########################################
    ### Functions to Interact with States ###
    #########################################
    def set_initial_state(self, name: str='initial', debug: int=0):
        """Create the initial state for this cell."""
        ini_state = self.add_state(name)
        ini_state.set_geometry(self.labels, self.coord)
        ini_state.set_cell(self.cell_vector, self.cell_param)
        ini_state.get_molecules(debug=debug)
        return ini_state
        
    ######
    def add_state(self, name: object, debug: int=0):
        """Create or return a state with the requested name."""
        from scope.classes_state import State
        if not hasattr(self,"states"): setattr(self,"states",list([]))
        exists, new_state = self.find_state(name)
        if exists:  
            if debug > 0: print(f"CELL.ADD_STATE. State with same {name=} found, returning it")
            return new_state
        else:
            if debug > 0: print("CELL.ADD_STATE. Creating new state, returning it")
            new_state = State(self, name, debug=debug)
            self.states.append(new_state)
        return new_state
    
    ######
    def remove_state(self, search_name: str, debug: int=0):
        from scope.classes_state import State
        if not hasattr(self,"states"): return False, None
        if debug > 0: print(f"CELL.REMOVE_STATE: Searching {search_name} state in CELL with {list(st.name for st in self.states)} states")
        found = False
        for idx, st in enumerate(self.states):
            if st.name == search_name and not found: 
                found = True; found_idx = idx
                if debug > 0: print(f"CELL.REMOVE_STATE: state {search_name} found. Removing it")
        if found: 
            del self.states[found_idx]
        else:
            if debug > 0: print(f"CELL.REMOVE_STATE: state {search_name} not found")

    ######
    def find_state(self, search_name: str, debug: int=0):
        """Find a state by name.

        Parameters:
            search_name (str):           Name of the state to search.
            debug (int):                 Verbosity level.

        Returns:
            tuple: `(exists, state_obj)`.
        """
        from scope.classes_state import State
        if not hasattr(self,"states"): return False, None
        if debug > 0: print(f"CELL.FIND_STATE: Searching {search_name} in Cell object with {len(self.states)} states")
        for sta in self.states:
            if debug > 0: print(f"CELL.FIND_STATE: Comparing {search_name} with {sta.name}")
            if sta.name == search_name: 
                if debug > 0: print(f"CELL.FIND_STATE: state {search_name} found")
                return True, sta
        if debug > 0: print(f"CELL.FIND_STATE: state {search_name} not found")
        return False, None
    
    ################################
    ## Visualization and Printing ##
    ################################
    def view(self, size: str='default'):
        """
        Visualize the cell with Plotly.

        Parameters:
            size (str):                  Figure size preset.

        Returns:
            None
        """

        import plotly.graph_objects as go
        from scope.read_write import set_scene
        from scope.elementdata import ElementData  
        elemdatabase = ElementData()

        size_map = {'default': (600, 600, 8, 9), 'small': (400, 400, 6, 7), 
                    'large': (800, 800, 10, 12), 'ultra': (1000, 1000, 11, 13)}
        width, height, marker_size, text_size = size_map.get(size.lower(), size_map['default'])

        if not hasattr(self, "adjmat"): self.get_adjmatrix(adjust_factor=True)
        fig = go.Figure()
        positions, symbols = np.array(self.coord), self.labels
        unique_bonds = {tuple(i) for i in np.argwhere(self.adjmat > 0)}

        fig.add_trace(go.Scatter3d(
            x=positions[:, 0], y=positions[:, 1], z=positions[:, 2],
            mode='markers', marker=dict(size=marker_size, 
            color=[elemdatabase.cpk_colors[l] for l in symbols], 
            line=dict(color='black', width=1)), hoverinfo='text', text=symbols, showlegend=False))

        fig.add_trace(go.Scatter3d(
            x=positions[:, 0], y=positions[:, 1], z=positions[:, 2],
            mode='text', text=symbols, textfont=dict(color='black', size=text_size), 
            hoverinfo='none', showlegend=False))

        for i, j in unique_bonds:
            fig.add_trace(go.Scatter3d(
                x=[positions[i, 0], positions[j, 0]], y=[positions[i, 1], positions[j, 1]],
                z=[positions[i, 2], positions[j, 2]], mode='lines',
                line=dict(color='gray', width=5), hoverinfo='none', showlegend=False))

        set_scene(fig, positions, width=width, height=height)
        fig.show()

    def __repr__(self, indirect: bool=False):
        """Return a formatted summary of the cell."""
        to_print = ''   
        if not indirect: to_print += '-------------------------------\n'
        if not indirect: to_print += '   >>> SCOPE CELL Object >>>   \n'
        if not indirect: to_print += '-------------------------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.object_type}\n'
        to_print += f' SubType               = {self.object_subtype}\n'
        to_print += f' Name                  = {self.name}\n'
        to_print += f' Num Atoms             = {self.natoms}\n'
        to_print += f' Cell Parameters a:c   = {self.cell_param[0:3]}\n'
        to_print += f' Cell Parameters al:ga = {self.cell_param[3:6]}\n'
        if hasattr(self,"volume"): to_print += f' Volume (Angs^3)       = {np.round(self.volume,4)}\n'
        to_print += f' Number of Units (Z)   = {self.z}\n'
        #if hasattr(self,"z"):      to_print += f' Z                     = {self.z}\n'
        if hasattr(self,"molecules"):  
            to_print += f' Num Molecules:        = {len(self.molecules)}\n'
            to_print += f' With Formulae:                               \n'
            for idx, m in enumerate(self.molecules):
                to_print += f'    {idx}: {m.formula} \n'
        to_print += '-------------------------------\n'
        if hasattr(self,"ref_molecules"):
            to_print += f' Num of Ref Molecules: = {len(self.ref_molecules)}\n'
            to_print += f' With Formulae:                                  \n'
            for idx, ref in enumerate(self.ref_molecules):
                to_print += f'    {idx}: {ref.formula} \n'
        if self.charge != 0:
            to_print += f' WARNING: Charge is {self.charge}. It should be 0\n'
        if not indirect: to_print += '\n'
        return to_print

######################
####    IMPORT    ####
######################
def import_cell(old_cell: object, debug: int=0) -> object:
    from scope.classes_specie  import import_molecule
    assert hasattr(old_cell,"labels") 
    assert hasattr(old_cell,"coord")   or hasattr(old_cell,"pos")
    assert hasattr(old_cell,"refcode") or hasattr(old_cell,"name")

    if hasattr(old_cell,"warning_list"):        ## In cell2mol cells, this variable contains any warning found during the interpretation
        if any(old_cell.warning_list):          ## For a cell to be valid, no error is allowed
            return None

    labels     = old_cell.labels
    if   hasattr(old_cell,"coord"):       coord       = old_cell.coord
    elif hasattr(old_cell,"pos"):         coord       = old_cell.pos

    if   hasattr(old_cell,"name"):        name        = old_cell.name
    elif hasattr(old_cell,"refcode"):     name        = old_cell.refcode

    if   hasattr(old_cell,"cellvec"):     cell_vector = old_cell.cellvec
    elif hasattr(old_cell,"cell_vector"): cell_vector = old_cell.cell_vector
    else:                                 cell_vector = None

    if   hasattr(old_cell,"cell_param"):  cell_param  = old_cell.cell_param
    elif hasattr(old_cell,"cellparam"):   cell_param  = old_cell.cellparam
    else:                                 cell_param  = None

    if   hasattr(old_cell,"moleclist"):   old_molecules = old_cell.moleclist
    elif hasattr(old_cell,"molecules"):   old_molecules = old_cell.molecules

    if   hasattr(old_cell,"refmoleclist"):    old_ref_molecules = old_cell.refmoleclist
    elif hasattr(old_cell,"ref_molecules"):   old_ref_molecules = old_cell.ref_molecules
    else:                                     old_ref_molecules = None

    # If cell_param and cell_vector are None, the creation of the cell will fail
    new_cell            = Cell(name, labels, coord, cell_vector, cell_param)
    new_cell.object_subtype = "cell"
    new_cell.origin     = "import_cell"
    if debug > 0: print(f"IMPORT CELL: importing cell {new_cell.name}")

    ## Molecules
    if debug > 0: print(f"IMPORT CELL: importing molecules")
    new_molecules = []
    for mol in old_molecules: 
        new_mol = import_molecule(mol, parent=new_cell, debug=debug)
        new_mol.set_bonds()
        new_mol.fix_ligands_rdkit_obj(debug=debug)
        new_molecules.append(new_mol)

    if debug > 0: 
        print(f"IMPORT CELL: prepared molecules: (formula, charge, spin)")
        for mol in new_molecules:
            print(mol.formula, mol.charge, mol.spin)

    new_cell.set_molecules(new_molecules)
    new_cell.fix_cell_coord(debug=debug)  ## In cell2mol, the cell object does not have the coordinates of the reconstructed cell. 
                                          ## However, the molecule and atom objects are updated (i.e. reconstructed). We use this info to update the cell

    ## Reference molecules
    new_cell.ref_molecules = []
    if old_ref_molecules is not None:
        if debug > 0: print(f"IMPORT CELL: importing ref_molecules")
        for rmol in old_ref_molecules:
            new_cell.ref_molecules.append(import_molecule(rmol, parent=new_cell))
    else:
        if debug > 0: print(f"IMPORT CELL: creating ref_molecules")
        for mol in new_cell.molecules:
            found = False
            for rmol in new_cell.ref_molecules:
                if mol == rmol: found = True 
            if not found: new_cell.ref_molecules.append(mol)

    ## Temporary things that I'd like to remove from the import once sorted 
    if hasattr(old_cell,"warning_list"): new_cell.warning_list = old_cell.warning_list 

    return new_cell
