import numpy as np
from scope.connectivity      import *
from scope.classes_data      import Collection, Data
from scope.classes_specie    import *
from scope.operations.dicts_and_lists import extract_from_list
from scope.elementdata       import ElementData
elemdatabase = ElementData()

##############
### STATES ###
##############
class State(object):
    """
    State class for representing a physical or chemical state associated with a source (cell or specie).
    This class provides methods for managing geometries, molecules, computations, vibrational normal modes (VNMs) and thermodynamic data

    Initiate

    _source : object                    The parent object (cell, molecule, or specie) to which this state belongs.
    name : str                          Name of the state.

    Attributes

    type : str                          Type of the object ("state").
    name : str                          Name of the state.
    results : dict                      Dictionary of results associated with the state.
    computations : list                 List of computations associated with the state.
    labels : list                       Atomic labels for the geometry.
    coord : list                        Atomic positions (Cartesian coordinates).
    natoms : int                        Number of atoms in the state.
    formula : str                       Chemical formula, derived from labels.
    radii : list                        List of Atomic radii.
    cell_vector : array-like            Cell vectors for periodic systems.
    cell_param : array-like             Cell parameters.
    forces : array-like                 Forces on atoms.
    atoms : list                        List of atom objects.
    moleclist : list                    List of molecule objects in the state.
    ncomplex : int                      Number of Transition Metal Complexes in the state.
    z : int                             Number of stoichiometric Units in the state (Z).
    fragmented : bool                   Whether the state geometry is fragmented.
    VNMs : list                         Vibrational normal modes.

    Methods

    set_geometry()                      Set atomic labels and positions for the state geometry.
    set_geometry_from_moleclist()       Set geometry from the list of molecule objects.
    set_cell()                          Set cell vectors and parameters.
    set_forces()                        Set atomic forces.
    set_atoms()                         Set atom objects from molecule list.
    get_atoms()                         Retrieve list of atom objects from molecules.
    get_moleclist()                     Generate molecule list from geometry.
    get_ncomplex()                      Count number of complex molecules.
    get_z()                             Count number of stoichiometric units in the state (Z)
    reconstruct()                       Attempt to reconstruct fragmented state geometry.
    check_fragmentation()               Check if the state geometry is fragmented.
    check_minimum()                     Check if the state is a minimum or almost a minimum.
    set_VNMs()                          Set vibrational normal modes and evaluate minimum status.
    find_computation()                  Find a computation in the state.
    add_computation()                   Add a computation to the state.
    sample_geometries()                 Sample geometries around the current state using VNMs and furthest point sampling.
    add_result()                        Add a result to the state.
    set_energy()                        Set electronic energy result.
    set_Helec()                         Set electronic energy per stoichiometric unit (Z)
    get_thermal_data()                  Compute and store thermodynamic data (Helec, Selec, Hvib, Svib, Gtot).
    """
    def __init__(self, _source: object, name: str, debug: int=0):
        self.type         = "state"
        self._source      = _source
        self.name         = name
        self.results      = dict()
        self.computations = []

############################################
#### Basic Functions to add information ####
############################################
    def set_geometry(self, labels, coord, debug: int=0):
        assert len(labels) == len(coord)
        self.labels      = labels
        self.coord       = coord
        self.natoms      = len(labels)
        self.formula     = labels2formula(self.labels)
        self.radii       = get_radii(labels)
        self.get_moleclist(debug=debug) ## Molecules require updated

    def set_geometry_from_moleclist(self, overwrite: bool=False, debug: int=0):
        if not hasattr(self,"moleclist"): self.get_moleclist(overwrite=overwrite, debug=debug)
        self.labels     = []
        self.coord      = []
        indices         = []
        for mol in self.moleclist:
            if not hasattr(mol,"atoms"): mol.set_atoms()
            for idx, at in enumerate(mol.atoms):
                self.labels.append(at.label)
                self.coord.append(at.coord)
                indices.append(mol.indices[idx])
        ## Below is to order the atoms as in the original cell, using the indices stored in the molecule object
        self.labels = [x for _, x in sorted(zip(indices, self.labels), key=lambda pair: pair[0])]
        self.coord = [x for _, x in sorted(zip(indices, self.coord), key=lambda pair: pair[0])]
        self.natoms = len(self.labels)
        self.formula = labels2formula(self.labels)
        assert len(self.labels) == len(self.coord)
         
    def set_cell(self, cell_vector: list=None, cell_param: list=None):
        if   cell_vector is None and cell_param is None:
            raise ValueError("STATE.SET_CELL: Either cell_vector or cell_param must be provided to set the cell")
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
        self.volume               = get_unit_cell_volume(*self.cell_param) 

    def set_forces(self, forces):
        self.forces      = forces

#########################
#### Charge and Spin ####
#########################
    ## The Spin and Charge of a State is always taken from the Source (Specie or Cell)
    @property
    def charge(self):
        return self._source.charge
    @property
    def atomic_charges(self):
        return self._source.atomic_charges
    @property
    def spin(self):
        return self._source.spin
    @property
    def atomic_spins(self):
        return self._source.atomic_spins
    @property
    def ismagnetic(self):
        return self._source.ismagnetic
    @property
    def spin_multiplicity(self):
        return self._source.spin_multiplicity

##########################
#### Other Properties ####
##########################
    @property
    def Z(self):
        return self.z

    @property
    def z(self):
        if not hasattr(self, "_z"):
            return self.get_z()
        return self._z

###################################
#### Operations with Molecules ####
###################################
    def get_moleclist(self, overwrite: bool=False, debug: int=0):
        from scope.classes_specie import Molecule

        # Overwrite
        if not overwrite and hasattr(self,"moleclist"): 
            if debug > 0: print(f"STATE.GET_MOLECLIST. Moleclist already exists and default is overwrite=False")
            return self.moleclist
        # Security
        if not hasattr(self,"labels") or not hasattr(self,"coord"): return None
        if len(self.labels) == 0 or len(self.coord) == 0: return None
        # Covalent Factor
        if hasattr(self._source,"factor"):       cov_factor = self._source.factor
        elif hasattr(self._source,"cov_factor"): cov_factor = self._source.cov_factor
        else: cov_factor = 1.3
        # Metal Factor
        if hasattr(self._source,"metal_factor"):   metal_factor = self._source.metal_factor
        else: metal_factor = 1.0

        if debug > 0: print(f"STATE.GET_MOLECLIST. Sending split_species with {cov_factor=} and {metal_factor=}") 
        blocklist = split_species(self.labels, self.coord, cov_factor=cov_factor, metal_factor=metal_factor, debug=debug)
        self.moleclist = [] 
        for b in blocklist:
            if debug > 0: print(f"STATE.GET_MOLECLIST: doing block={b}")
            mol_labels      = extract_from_list(b, self.labels, dimension=1)
            mol_coord       = extract_from_list(b, self.coord, dimension=1)
            if hasattr(self,"frac_coord"): mol_frac_coord  = extract_from_list(b, self.frac_coord, dimension=1)
            else:                          mol_frac_coord  = None
            # Creates Molecule Object
            newmolec    = Molecule(mol_labels, mol_coord, mol_frac_coord)
            # For debugging
            newmolec.origin = "state.get_moleclist"
            # Adds State as parent of the molecule, with indices b
            newmolec.add_parent(self._source, indices=b, debug=debug)            
            # Store Adjacency Parameters
            newmolec.set_factors(cov_factor, metal_factor)
            # Creates The atom objects with adjacencies
            newmolec.set_atoms(create_adjacencies=True, debug=debug)
            # Sets Fractional Coordinates
            if hasattr(self,"frac_coord"): 
                newmolec.set_fractional_coord(mol_frac_coord)
            elif self._source.type == "cell": 
                if debug > 0: print(f"STATE.GET_MOLECLIST: getting frac_coord from _source of type=cell")
                newmolec.get_fractional_coord(cell_vector = self._source.cell_vector)
            # The split_complex must be below the frac_coord, so they are carried on to the ligands    
            if newmolec.iscomplex: 
                if debug > 0: print(f"STATE.GET_MOLECLIST: splitting complex")
                newmolec.split_complex(debug=debug)
            self.moleclist.append(newmolec)
        return self.moleclist

    ######
    def get_ncomplex(self, debug: int=0):
        ## Returns the number of Transition Metal Complexes (TMC) in the unit cell. 
        ## Gradually replacing it by Z (computed with get_z).
        if debug > 0: print(f"STATE.GET_NCOMPLEX checking fragmentation")
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"Found Fragmented molecules in the geometry of state: {self.name}"
         
        if debug > 0: print(f"STATE.GET_NCOMPLEX getting moleclist")
        if not hasattr(self,"moleclist"): self.get_moleclist(debug=debug)
        self.ncomplex = 0
        for mol in self.moleclist:
            if mol.iscomplex: self.ncomplex += 1
        if debug > 0: print(f"State.get_ncomplex {self.ncomplex} complexes found in state: {self.name}")
        return self.ncomplex

    ######
    def get_z(self, debug: int=0):
        ## Returns the number of stoichiometric units in the unit cell. 
        ## Basically, how many times the same stoichiometry unit is repeated in the cell
        from scope.operations.vecs_and_mats import gcd_list
        if debug > 0: print(f"STATE.GET_Z checking fragmentation")
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"STATE.GET_Z found Fragmented molecules in the geometry of self: {self.name}"

        if debug > 0: print(f"STATE.GET_Z getting moleclist")
        if not hasattr(self,"moleclist"): self.get_moleclist(debug=debug)

        unique = [] 
        occurrences = []
        for mol in self.moleclist:
            found = False
            for uni in unique:
                if mol == uni: found = True
            if not found: 
                unique.append(mol)
                occurrences.append(self.get_occurrence(mol))
        self._z = int(gcd_list(occurrences))
        return self._z 
    
    ######
    def get_atoms(self, debug: int=0):
        ## Retrieves a list of atoms, extracted from the molecules in moleclist.
        ## The challenge is that the atoms do not necessarily appear in the same order as in labels/coord, so we need to reorder them
        if not hasattr(self,"moleclist"): 
            if debug > 0: print(f"STATE.GET_ATOMS: generating moleclist")
            self.get_moleclist(debug=debug)

        self.atoms = []
        tmp_indices = []
        for mol in self.moleclist:
            for at in mol.atoms:
                tmp_indices.append(at.get_parent_index(self._source.subtype))  ## We get the index of the atom in the source of this state
                self.atoms.append(at)
        self.atoms = [x for _, x in sorted(zip(tmp_indices, self.atoms), key=lambda pair: pair[0])]
        return self.atoms

    ######
    def get_occurrence(self, substructure: object, debug: int=0) -> int:
        """
        Counts the number of times a given substructure appears inside the State.
        Args:
            substructure (object): The substructure to search for within the State.
            debug (int, optional): Debug level for comparison functions. Defaults to 0.
        Returns:
            int: The number of times the substructure appears within the State.
        """
        ## Finds how many times a substructure appears in self
        occurrence = 0

        if debug > 0: print(f"STATE.GET_OCCURRENCE checking fragmentation")
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"STATE.GET_OCCURRENCE found fragmented molecules in the geometry of cell: {self.name}"
         
        if debug > 0: print(f"STATE.GET_OCCURRENCE getting moleclist")
        if not hasattr(self,"moleclist"): self.get_moleclist(debug=debug)

        ## Case of Species inside self
        if hasattr(substructure,"type"):
            if substructure.type == 'specie':
                for mol in self.moleclist:
                    if mol.__eq__(substructure, with_graph=True): occurrence += 1
        return occurrence

########################
#### Reconstruction ####
########################
    def reconstruct(self, debug: int=0):
        from scope.reconstruct import classify_fragments, fragments_reconstruct 
        if not self._source.type == "cell":          raise ValueError(f"STATE_RECONSTRUCT: state's source should by a CELL object") 
        if not hasattr(self,"cell_vector"):          raise ValueError(f"STATE_RECONSTRUCT: state should have a cell vector") 
        if not hasattr(self._source,"refmoleclist"): raise ValueError(f"STATE.RECONSTRUCT: state's source does not have a list of reference molecules"); return None
        from scope.read_write import HiddenPrints
        if debug > 0: print("STATE.RECONSTRUCT: reconstructing cell of state", self.name)
        with HiddenPrints():
            finished = False
            if not hasattr(self,"moleclist"): self.get_moleclist(debug=debug) 
            import itertools
            blocklist    = self.moleclist.copy()
            refmoleclist = self._source.refmoleclist.copy()
            cov_factor   = refmoleclist[0].cov_factor
            metal_factor = refmoleclist[0].metal_factor
            moleclist, fraglist, Hlist = classify_fragments(blocklist, refmoleclist, debug=debug) 
            if len(fraglist) > 0 or len(Hlist) > 0: 
                moleclist, finalmols, Warning = fragments_reconstruct(moleclist,fraglist,Hlist,refmoleclist,self.cell_vector,cov_factor,metal_factor, debug=debug)
                moleclist.extend(finalmols)
                self.moleclist = moleclist
                self.set_geometry_from_moleclist()
                finished = True
        if debug > 0 and finished: print("STATE.RECONSTRUCT: state reconstructed succesfully")
        return self.moleclist

    def check_fragmentation(self, reconstruct: bool = False, debug: int=0):
        ## If the source is a specie, then in principle there should only be one. So with moleclist is enough to find
        if self._source.type == "specie":
            if not hasattr(self,"moleclist"): self.get_moleclist(debug=debug)
            if len(self.moleclist) > 1: self.fragmented = True
            else:                       self.fragmented = False
            if debug > 0: print(f"STATE.CHECK_FRAGMENTATION: source type=specie. {self.fragmented=}")
            return self.fragmented
        ## If it is a unit cell, then we need a list of molecules that should in principle be there. This is refmoleclist
        ## If the cell is created by cell2mol, then this list is already stored in the .cell object
        elif self._source.type == "cell": 
            assert hasattr(self,"cell_vector")
            assert hasattr(self._source,"refmoleclist")
            if not hasattr(self,"moleclist"): self.get_moleclist(debug=debug)
            self.fragmented = False
            # First comparison with current moleclist
            for mol in self.moleclist:
                found = False
                for rmol in self._source.refmoleclist:
                    if mol.__eq__(rmol,with_graph=False): found = True # Graph cannot be used for rmol, since it doesn't have rdkit object
                if not found: self.fragmented = True
            # If there are fragments and user wants reconstruction, it tries to reconstruct and checks the new moleclist
            if self.fragmented and reconstruct:
                new_moleclist = self.reconstruct(debug=debug)
                self.fragmented = False
                for mol in new_moleclist:
                    found = False
                    for rmol in self._source.refmoleclist:
                        if mol.__eq__(rmol,with_graph=False): found = True # Graph cannot be used for rmol, since it doesn't have rdkit object
                    if not found: self.fragmented = True
                if not self.fragmented: 
                    self.moleclist = new_moleclist
                    self.set_geometry_from_moleclist(debug=debug)
        else: print(f"STATE.CHECK_FRAGMENTATION: Unknown source Type {self._source.type}")
        return self.fragmented

##############################
#### Connection with VNMs ####
##############################
    def check_minimum(self, debug: int=0):
        if hasattr(self,"isminimum"): 
            if self.isminimum or self.almost_minimum: return True
            else:                                     return False
        else:
            if not hasattr(self,"VNMs"): 
                if debug > 0: 
                    if hasattr(self._source,"name"): print(f"STATE.check_minimum: state {self.name} of {self._source.name} does not have VNMs")
                    else:                            print(f"STATE.check_minimum: state {self.name} of {self._source.formula} does not have VNMs")
                return False
            else:
                self.set_VNMs(self.VNMs)
                if self.isminimum or self.almost_minimum: return True
                else:                                     return False

    def set_VNMs(self, VNMs):
        self.VNMs       = VNMs
        self.freqs_cm   = [vnm.freq_cm for vnm in VNMs]
        if all(vnm.freq_cm >= 0.0 for vnm in self.VNMs): self.isminimum = True
        else:                                            self.isminimum = False
        ## If it is not a minimum, evaluates if, at least, is close
        if not self.isminimum:
            self.num_neg_freqs = 0
            for vnm in self.VNMs: 
                if vnm.freq_cm < 0.0: self.num_neg_freqs += 1
            if self.num_neg_freqs <= 3 and VNMs[0].freq_cm > -50: self.almost_minimum = True
            else:                                                 self.almost_minimum = False

########################################
#### Connection with Excited States ####
########################################
    def set_exc_states(self, exc_states, debug: int=0):
        self.set_exc_states = exc_states
        return self.set_exc_states

##################################
#### Connection with Workflow ####
##################################
    def find_computation(self, job_keyword: str='', step: int=1, run_number: int=1, debug: int=0):
        for idx, comp in enumerate(self.computations):
            if comp._job.keyword == job_keyword and comp.step == step and comp.run_number == run_number: this_comp = comp; return True, this_comp
        return False, None

    def add_computation(self, computation: object, debug: int=0):
        found, comp = self.find_computation(computation._job.keyword, computation.step, computation.run_number)
        if not found: 
            if debug > 0: print("STATE.ADD_COMPUTATION: same computation wasn't found. So adding it to state")
            self.computations.append(computation)
            computation.add_state(self)
        else:
            if debug > 0: print("STATE.ADD_COMPUTATION: same computation was already found in state. Ignoring")

#############################
#### Sampling Geometries ####
#############################
    def sample_geometries(self, ngeoms: int, n_aux_geoms=100, n_fps_rounds=0, temp: float=300, sigma_damp_factor: float=1, freq_bottom_limit: float=50, debug: int=0):
        """
        - The sampling parameters (temperature, number of rounds, samples per round) can be adjusted.
        Samples geometries around the current state geometry. It uses vibrational normal modes (VNMs), and explores the geometries in Q space. 
        If n_fps_rounds = 0, the distribution of output geometries should naturally follow a Wigner-like energetic criterion and be representative of the molecular motion
        If n_fps_rounds > 0, the distribution of output geometries will prioritize distant geometries, selected using a furthest point sampling (FPS) on the Q-displacement arrays.
        Parameters
        ----------
        ngeoms: int
            Number of geometries to select per round (default is 10). This is also the final number of geometries you will receive
        temp : float, optional
            Temperature parameter for sampling (default is 300).
        n_aux_geoms: int, optional
            Number of geometries to sample per geometry selected in FPS round (default is 100). If n_fps_rounds == 0, it won't be used
        n_fps_rounds : int, optional
            Number of FPS sampling rounds (default is 0).
        debug : int, optional
            Debug level for verbose output (default is 0).

        Returns
        -------
        current_q_disp : list of np.ndarray
            List of displacement vectors (in normal mode coordinates) for the selected geometries.
        current_geoms : list of np.ndarray
            List of sampled geometries (Cartesian coordinates).

        Raises
        ------
        ValueError
            If the state is not a minimum or if VNMs do not have eigenvectors.

        Notes
        -----
        - The method performs several rounds of sampling, each time generating new geometries by perturbing the current ones along their normal modes.
        - Furthest point sampling (FPS) is used to select a diverse subset of geometries in each round.
        - The sampling parameters (temperature, number of rounds, samples per round) depend on the `typ` argument.
        - Requires that the state is a minimum and that VNMs have eigenvectors parsed.
        """
        from scope.vnm_tools import geom_sampling_from_vnm, euclidean_q_distance, custom_q_distance, beta_distance
        from scope.other import furthest_point_sampling

        if not self.check_minimum():
            raise ValueError("State is not a minimum")
        if not hasattr(self.VNMs[0],"has_mode"):
            raise ValueError("VNMs do not have Eigenvectors. Please parse them")

        if debug > 0:
            if n_fps_rounds == 0: 
                print(f"-------------------------------------------------------------------------------------------------")
                print(f"STATE.SAMPLE_GEOMETRIES: {ngeoms} geometries will be generated directly from the State's geometry")
                print(f"-------------------------------------------------------------------------------------------------")
            else:                 
                print(f"-------------------------------------------------------------------------------------------------")
                print(f"STATE.SAMPLE_GEOMETRIES: {ngeoms} geometries will be generated in two steps:")
                print(f"STATE.SAMPLE_GEOMETRIES: 1) An initial sampling starting from the State's geometry, generating {ngeoms} geometries")
                print(f"STATE.SAMPLE_GEOMETRIES: 2) From each of the resulting {ngeoms} geometries, another sampling will be performed in which {n_aux_geoms} will be generated")
                print(f"STATE.SAMPLE_GEOMETRIES: 3) Step 2 will be repeated for {n_fps_rounds} rounds. At each round, {n_aux_geoms**2} will be created.")
                print(f"STATE.SAMPLE_GEOMETRIES: 4) After the last round, {ngeoms} will be selected")
                print(f"-------------------------------------------------------------------------------------------------")

        q_min = np.zeros((len(self.VNMs)))
        current_geoms       = [] 
        current_q_disp      = [] 
        current_energies    = [] 
        current_geoms.append(self.coord)
        current_q_disp.append(q_min)
        current_energies.append(float(0.0))

        ## Main Loop
        for nr in range(n_fps_rounds+1):
            q_fps, g_fps, e_fps = [], [], [] 
            count = 0
            for c, q in zip(current_geoms, current_q_disp):
                count += 1

                if nr == 0: geoms_this_round = ngeoms
                else:       geoms_this_round = n_aux_geoms

                if debug > 0: print(f"STATE.SAMPLE_GEOMETRIES: Running sampling of initial geometry, with: {geoms_this_round=} and {debug=}") 
                geoms, q_disp, energies = geom_sampling_from_vnm(self.labels, c, self.VNMs, qini=q, T=temp, n_samples=geoms_this_round, sigma_damp_factor=sigma_damp_factor, freq_bottom_limit=freq_bottom_limit, check_adjacencies=True, debug=debug)
                if debug > 0: 
                    if nr == 0: 
                        print(f"STATE.SAMPLE_GEOMETRIES: Initial structure {count}/{len(current_geoms)} sampled {len(q_disp)} geometries")
                    else:
                        print(f"STATE.SAMPLE_GEOMETRIES: Initial structure {count}/{len(current_geoms)} of FPS round {nr}/{n_fps_rounds+1} sampled {len(q_disp)} geometries")

                # Data for FPS
                if nr < n_fps_rounds:  ## Minimum (i.e. initial structure, with Q=0) is added in every round except the last one
                    q_fps.append(q_min)
                    g_fps.append(self.coord)
                    e_fps.append(float(0.0))
                q_fps.extend(q_disp)
                g_fps.extend(geoms)
                e_fps.extend(energies)

            if len(q_fps) > ngeoms:
                #Run FPS
                if debug > 0: print(f"STATE.SAMPLE_GEOMETRIES: Entering FPS selection with {len(q_fps)} geometries. Selecting {ngeoms}")
                idxs = furthest_point_sampling(q_fps, ngeoms, euclidean_q_distance)
                if debug > 0: print(f"STATE.SAMPLE_GEOMETRIES: FPS of round {nr+1}/{n_fps_rounds+1} kept {len(idxs)} geometries, with indices:{idxs}")
            else:
                if debug > 0 and n_fps_rounds > 0: print(f"STATE.SAMPLE_GEOMETRIES: Sampling of round {nr+1}/{n_fps_rounds+1} failed to generate enough samples for FPS. Taking all available to the next round")
                idxs = list(range(len(q_disp)))

            # Prepares Next Round 
            current_geoms    = [] 
            current_q_disp   = [] 
            current_energies = [] 
            for idx in idxs:
                current_geoms.append(g_fps[idx])
                current_q_disp.append(q_fps[idx])
                current_energies.append(e_fps[idx])

        return current_q_disp, current_geoms, current_energies 

#######################################
#### Results associated with State ####
#######################################
    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  
            self.results[result.key] = result

    def set_energy(self, energy, units, overwrite: bool=True):
        self.add_result(Data("energy",energy,units,"state.set_energy()"), overwrite=overwrite)

    def set_Helec(self, overwrite: bool=True, debug: int=0):
        assert "energy" in self.results
        if not hasattr(self,"z"): self.get_z(debug=debug)
        self.add_result(Data("Helec",self.results["energy"].value/self.z,self.results["energy"].units,"state.set_Helec()"), overwrite=overwrite)
        
################################
#### Get Thermodynamic Data ####
################################
    def get_thermal_data(self, Trange: range=range(10,501,1), Helec=None, Selec=None, Hvib=None, Svib=None, Gtot=None, overwrite: bool=False, debug: int=0):
        from scope.thermal_corrections import get_Selec, get_Hvib, get_Svib, get_Gibbs

        if not hasattr(self,"z"): self.get_z(debug=debug)
        if debug > 0:           print(f"STATE.GET_THERMAL_DATA: found {self.z} stoichiometric units")

        if Hvib is None and Svib is None:
            assert hasattr(self,"isminimum"), f"I can't compute thermal data on this state. Missing VNMs"
        if not self.results["energy"]:
            raise ValueError(f"STATE.GET_THERMAL_DATA: missing State energy value")
    
        ############## Helec ##############
        if Helec is None:   ### One can provide specific values for Helec, Selec, Hvib, Svib and Gtot 
            if overwrite or not "Helec" in self.results.keys():
                self.add_result(Data("Helec",self.results["energy"].value/self.z,self.results["energy"].units,"state.get_thermal_data()"), overwrite=overwrite)
        else: 
            if isinstance(Helec, Data):
                if overwrite or not "Helec" in self.results.keys():
                    self.add_result(Data("Helec",Helec.value,Helec.units,"enforced in state.get_thermal_data()"), overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Helec. It must be a DATA-class object")
        if debug > 0: print(f"Helec is {self.results['Helec']}")

        ############## Selec ##############
        if Selec is None:
            if overwrite or not "Selec" in self.results.keys():
                self.add_result(get_Selec(self.spin_multiplicity, outunits='au', nmol=self.z), overwrite=overwrite)
        else: 
            if isinstance(Selec, Data):
                if overwrite or not "Selec" in self.results.keys():
                    self.add_result(Data("Selec",Selec.value,Selec.units,"enforced in state.get_thermal_data()"), overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Selec. It must be a DATA-class object")
        if debug > 0: print(f"Selec is {self.results['Selec']}")

        ############## Hvib ##############
        if Hvib is None:
            if overwrite or not "Hvib" in self.results.keys():
                Hvib = Collection("Hvib", "Temperature")
                for temp in Trange:
                    Hvib.add_data(get_Hvib(np.abs(self.freqs_cm), temp, freq_units='cm', outunits='au', nmol=self.z))
                self.add_result(Hvib, overwrite=overwrite)
        else: 
            if isinstance(Hvib, Collection):
                if overwrite or not "Hvib" in self.results.keys():
                    self.add_result(Hvib, overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Hvib. It must be a COLLECTION-class object")
        if debug > 0: print(f"Hvib is {self.results['Hvib']}")

        ############## Svib ##############
        if Svib is None:
            if overwrite or not "Svib" in self.results.keys():
                Svib = Collection("Svib", "Temperature")
                for temp in Trange:
                    Svib.add_data(get_Svib(np.abs(self.freqs_cm), temp, freq_units='cm', outunits='au', nmol=self.z))
                self.add_result(Svib, overwrite=overwrite)
        else: 
            if isinstance(Svib, Collection):
                if overwrite or not "Svib" in self.results.keys():
                    self.add_result(Svib, overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Svib. It must be a COLLECTION-class object")
        if debug > 0: print(f"Svib is {self.results['Svib']}")

        ############## Gtot ##############
        if Gtot is None:
            if overwrite or not "Gtot" in self.results.keys():
                Gtot = Collection("Gtot", "Temperature")
                for temp in Trange:
                    # Retrieve data (not value)
                    Helec = self.results["Helec"]
                    Selec = self.results["Selec"]
                    Hvib_i = Hvib.find_value_with_property("temperature", temp)
                    Svib_i = Svib.find_value_with_property("temperature", temp)
                    assert Helec.units == Selec.units == Hvib_i.units == Svib_i.units, f"{Helec.units=}, {Selec.units=}, {Hvib_i.units=}, {Svib_i.units=}"
                    key = "Gtot"
                    value = get_Gibbs(Helec.value, Hvib_i.value, Selec.value, Svib_i.value, temp)
                    units = Helec.units
                    function = "state.get_thermal_data()"
                    new_data = Data(key, value, units, function)
                    new_data.add_property("temperature", temp, overwrite=overwrite)
                    Gtot.add_data(new_data)
                self.add_result(Gtot, overwrite=overwrite)
        else: 
            if isinstance(Gtot, Collection):
                if overwrite or not "Gtot" in self.results.keys():
                    self.add_result(Gtot, overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Gtot. It must be a COLLECTION")
        if debug > 0: print(f"Gtot is {self.results['Gtot']}")

    ######
    def compute_PV_term(self, pressure: float = 101.325, overwrite: bool=False, debug: int=0):
        from scope import constants
        # This function computes the PV term in kJ/mol, given a pressure in kilo-pascal
        # It is only valid for states whose source is a cell, as it uses its volume. 

        # Volume in angs^3
        # Pressure in kilo-pascal (10e3 Pa). The default is 1 atm = 101.325 kPa

        # It gets the stoichiometry number (Z)
        if not hasattr(self,"z"): self.get_z(debug=debug)

        # If the source is not a cell, it doesn't make sense to compute the PV term, so we return 0.0 kJ as a default value.
        if self._source.type != 'cell': 
            data = Data("PV",float(0.0),'kj',"state.compute_PV_term()")
            return data                     
        
        vm3 = self.volume * 1e-30 * constants.bohr2angs**3              ## Convert volume to m^3
        ppa = float(pressure) * 1e+6                                    ## Convert pressure to Pa 
        pv  = (ppa * vm3)                                               ## [Pa·m3] = [Joule] 
        pv *= constants.avogadro / 1000 / self.z                        ## kJ/molecule
        if overwrite or not "PV" in self.results.keys():
            data = Data("PV",pv,'kj',"state.compute_PV_term()")
            self.add_result(data, overwrite=overwrite)
        return data 

#######################
#### Visualization ####
#######################
    def __repr__(self) -> None:
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   STATE                                           \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Name                  = {self.name}\n'
        if hasattr(self._source,"name"):   to_print += f' Source Name           = {self._source.name}\n'
        if hasattr(self._source,"type"):   to_print += f' Source Type           = {self._source.type}\n'
        if hasattr(self,"labels"):         to_print += f' Labels                = {self.labels[0]}...\n'
        if hasattr(self,"coord"):          to_print += f' Coord                 = {self.coord[0]}...\n'
        if hasattr(self,"z"):              to_print += f' Number of Units (Z)   = {self.z}\n' 
        if hasattr(self,"isminimum"):      to_print += f' Is Minimum            = {self.isminimum}\n'
        if hasattr(self,"almost_minimum"): to_print += f' Almost a Minimum      = {self.almost_minimum}\n'
        if hasattr(self,"freq_cm"):        to_print += f' First Frequency (cm-1)= {self.freq_cm[0]}...\n'
        if hasattr(self,"moleclist"):  
            to_print += f' # Molecules:          = {len(self.moleclist)}\n'
            to_print += f' With Formulae:                               \n'
            for idx, m in enumerate(self.moleclist):
                to_print += f'    {idx}: {m.formula} \n'
        return to_print

##############################################################################
## Generic Find_State Function. There are specific CELL and SPECIE class functions ### 
##############################################################################
def find_state(source: object, search_name: str, debug: int=0):
    if debug >= 1: print("FIND_STATE: enters",search_name," with", len(source.states),"states in source")
    if not hasattr(source,"states"): return False, None
    else: 
        for idx, sta in enumerate(source.states):
            if sta.name == search_name: 
                if debug >= 1: print(f"FIND STATE: state {search_name} found")
                return True, sta
        if debug >= 1: print(f"FIND STATE: state {search_name} not found")
        return False, None
