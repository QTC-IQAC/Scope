import numpy as np
from Scope.Connectivity      import *
from Scope.Classes_Data      import collection, data
from Scope.Classes_Specie    import *
from Scope.Reconstruct       import *
from Scope.Elementdata       import ElementData
elemdatabase = ElementData()

##############
### STATES ###
##############
class state(object):
    """
    State class for representing a physical or chemical state associated with a source (cell, molecule, or specie).
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
    spin_config : object                Spin configuration.
    atoms : list                        List of atom objects.
    moleclist : list                    List of molecule objects in the state.
    ncomplex : int                      Number of Transition Metal Complexes in the state.
    fragmented : bool                   Whether the state geometry is fragmented.
    VNMs : list                         Vibrational normal modes.

    Methods

    set_geometry()                      Set atomic labels and positions for the state geometry.
    set_geometry_from_moleclist()       Set geometry from the list of molecule objects.
    set_cell()                          Set cell vectors and parameters.
    set_forces()                        Set atomic forces.
    set_spin_config()                   Set spin configuration.
    set_atoms()                         Set atom objects from molecule list.
    get_moleclist()                     Generate molecule list from geometry.
    get_ncomplex()                      Count number of complex molecules.
    get_SCO_geom()                      Print spin crossover geometry for complex molecules.
    reconstruct()                       Attempt to reconstruct fragmented state geometry.
    check_fragmentation()               Check if the state geometry is fragmented.
    check_minimum()                     Check if the state is a minimum or almost a minimum.
    set_VNMs()                          Set vibrational normal modes and evaluate minimum status.
    find_computation()                  Find a computation in the state.
    add_computation()                   Add a computation to the state.
    sample_geometries()                 Sample geometries around the current state using VNMs and furthest point sampling.
    add_result()                        Add a result to the state.
    set_energy()                        Set electronic energy result.
    set_Helec()                         Set electronic energy per complex.
    get_thermal_data()                  Compute and store thermodynamic data (Helec, Selec, Hvib, Svib, Gtot).

    - The class is tightly integrated with other classes in Scope (e.g., Classes_Scecie).
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
        self.labels      = []
        self.coord       = []
        indices     = []
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
         
    def set_cell(self, cell_vector, cell_param):
        self.cell_vector  = cell_vector
        self.cell_param   = cell_param

    def set_forces(self, forces):
        self.forces      = forces

#########################
#### Charge and Spin ####
#########################
    def set_spin_config(self, spins: list | int, typ: str='metals', debug: int=0):
        ## Verbose
        if debug > 0: 
            print(f"STATE.SET_SPIN_CONFIG: Preparing Spin Configuration for State {self.name}")
            print(f"STATE.SET_SPIN_CONFIG: Received spins", spins)
            print(f"STATE.SET_SPIN_CONFIG: Received typ", typ)
        ## Checks
        if typ == 'metals':
            if isinstance(spins, int): spins = [spins] * self.get_ncomplex()
            assert len(spins) == self.get_ncomplex(), f"STATE.SET_SPIN_CONFIG: number of spins provided ({len(spins)}) does not match number of complexes in state ({self.get_ncomplex()})"
        elif typ != 'metals':
            assert len(spins) == len(self.moleclist), f"STATE.SET_SPIN_CONFIG: number of spins provided ({len(spins)}) does not match number of molecules in state ({len(self.moleclist)})"

        ## Allocates the list of spins to the molecules in moleclist. Each item in spins corresponds to a molecule in the cell 
        pointer = 0
        if not hasattr(self,"moleclist"): self.get_moleclist()
        for idx, mol in enumerate(self.moleclist):
            if typ == 'metals' and mol.iscomplex: 
                print(f"STATE.SET_SPIN_CONFIG: Setting spin={spins[pointer]} to the metal of molecule {mol.formula} in index {idx}")
                print(f"STATE.SET_SPIN_CONFIG: Formal SPIN is added to the first metal of the molecule: {mol.metals[0].label}")
                mol.metals[0].set_spin(spins[pointer]); pointer += 1  
            elif typ != 'metals':                 
                print(f"STATE.SET_SPIN_CONFIG: Setting spin={spins[pointer]} to molecule {mol.formula} in index {idx}")
                print(f"STATE.SET_SPIN_CONFIG: Formal SPIN is added to the first atom of the molecule: {mol.atoms[0].label}")
                mol.atoms[0].set_spin(spins[pointer]); pointer += 1

    @property
    def charge(self):
        if not hasattr(self,"atoms"): self.get_atoms()
        if not all(hasattr(at,"charge") for at in self.atoms): return int(0)
        return np.sum(at.charge for at in self.atoms)

    @property
    def atomic_charges(self):
        if not hasattr(self,"atoms"): self.get_atoms()
        if not all(hasattr(at,"charge") for at in self.atoms): return [0]*self.natoms
        return [at.charge for at in self.atoms]

    @property
    def spin(self):
        if not hasattr(self,"atoms"): self.get_atoms()
        if not all(hasattr(at,"spin") for at in self.atoms): return int(0)
        return np.sum(at.spin for at in self.atoms)

    @property
    def atomic_spins(self):
        if not hasattr(self,"atoms"): self.get_atoms()
        if not all(hasattr(at,"spin") for at in self.atoms): return [0]*self.natoms
        return [at.spin for at in self.atoms]

    @property
    def ismagnetic(self):
        self.ismagnetic = False
        for at_s in self.atomic_spins:
            if at_s != 0: return True
        return False
   
    @property
    def spin_multiplicity(self):
        self.spin_multiplicity = int(self.spin + 1)
        return self.spin_multiplicity 

###################################
#### Operations with Molecules ####
###################################
    def get_moleclist(self, overwrite: bool=False, debug: int=0):
        from .Classes_Specie import molecule

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
            newmolec    = molecule(mol_labels, mol_coord, mol_frac_coord)
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
    def get_atoms(self):
        if not hasattr(self,"moleclist"): self.get_moleclist()
        self.atoms = []
        for mol in self.moleclist:
            for at in mol.atoms:
                self.atoms.append(at)
        return self.atoms

####################################
#### Specific to Spin Crossover ####
####################################
    def get_SCO_geom(self, debug: int=0):
        from .Spin_Crossover.SCO_Structure import geom_sco_from_xyz
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"Found Fragmented molecules in the geometry of state: {self.name}"
        if not hasattr(self,"moleclist"): self.get_moleclist(debug=debug)
        for mol in self.moleclist:
            if mol.iscomplex: print(geom_sco_from_xyz(self.labels, self.coord, debug=debug)) 

########################
#### Reconstruction ####
########################
    def reconstruct(self, debug: int=0):
        assert hasattr(self,"cell_vector")
        if not hasattr(self._source,"refmoleclist"): print("CLASS STATE.RECONSTRUCT: _source does not have refmoleclist"); return None
        from .Read_Write import HiddenPrints
        if debug > 0: print("CLASS_STATE.RECONSTRUCT: reconstructing cell of state", self)
        with HiddenPrints():
            finished = False
            if hasattr(self._source,"type") and hasattr(self,"cell_vector"):
                if not hasattr(self,"moleclist"): self.get_moleclist(debug=debug) 
                if self._source.type.lower() == "cell":
                    import itertools
                    blocklist = self.moleclist.copy()
                    refmoleclist = self._source.refmoleclist.copy()
                    cov_factor = refmoleclist[0].cov_factor
                    metal_factor = refmoleclist[0].metal_factor
                    moleclist, fraglist, Hlist = classify_fragments(blocklist, refmoleclist, debug=debug) 
                    if len(fraglist) > 0 or len(Hlist) > 0: 
                        moleclist, finalmols, Warning = fragments_reconstruct(moleclist,fraglist,Hlist,refmoleclist,self.cell_vector,cov_factor,metal_factor, debug=debug)
                        moleclist.extend(finalmols)
                        self.moleclist = moleclist
                        self.set_geometry_from_moleclist()
                        finished = True
     
                else: print("WARNING: reconstruct state, _source is not a cell. I will not reconstruct"); return None
            else: print("WARNING: reconstruct state, _source does not have 'type' or 'cell_vector' variables"); return None
        if debug > 0 and finished: print("CLASS STATE.RECONSTRUCT: state reconstructed to", self)
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
                    issame = compare_species(mol, rmol)  
                    if issame: found = True
                if not found: self.fragmented = True
            # If there are fragments and user wants reconstruction, it tries to reconstruct and checks the new moleclist
            if self.fragmented and reconstruct:
                new_moleclist = self.reconstruct(debug=debug)
                self.fragmented = False
                for mol in new_moleclist:
                    found = False
                    for rmol in self._source.refmoleclist:
                        issame = compare_species(mol, rmol)  
                        if issame: found = True
                    if not found: self.fragmented = True
                if not self.fragmented: 
                    self.moleclist = new_moleclist
                    self.set_geometry_from_moleclist(debug=debug)
        elif self._source.type == "perxyz": 
            if not hasattr(self,"moleclist"): self.get_moleclist(debug=debug)
            if len(self.moleclist) > 1:       self.fragmented = True
            else:                             self.fragmented = False
            if debug > 0: print(f"STATE.CHECK_FRAGMENTATION: source type={self._source.type}. {self.fragmented=}")
            return self.fragmented
        else:
            print(f"STATE.CHECK_FRAGMENTATION: Unknown source Type {self._source.type}")
        return self.fragmented

##############################
#### Connection with VNMs ####
##############################
    def check_minimum(self):
        if hasattr(self,"isminimum"): 
            if self.isminimum or self.almost_minimum: return True
            else:                                     return False
        else:
            if not hasattr(self,"VNMs"): 
                print(f"STATE.check_minimum: state does not have VNMs")
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
    def sample_geometries(self, n_selected: int=10, temp: float=300, n_samples_round=100, n_rounds=2, debug: int=0):
        """
        - The sampling parameters (temperature, number of rounds, samples per round) can be adjusted.
        Samples geometries around the current state geometry using vibrational normal modes (VNMs) and furthest point sampling (FPS).

        Parameters
        ----------
        n_selected : int, optional
            Number of geometries to select per round (default is 10).
        temp : float, optional
            Temperature parameter for sampling (default is 300).
        n_samples_round : int, optional
            Number of geometries to sample per round (default is 100).
        n_rounds : int, optional
            Number of sampling rounds (default is 2).
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
        from Scope.VNM_tools import geom_sampling_from_vnm, euclidean_q_distance, custom_q_distance, beta_distance
        from Scope.Other import furthest_point_sampling
        #if   typ.lower() == 'light':    temp=100; n_rounds=2; n_samples_round=100
        #elif typ.lower() == 'default':  temp=200; n_rounds=5; n_samples_round=300
        #elif typ.lower() == 'heavy':    temp=300; n_rounds=8; n_samples_round=500
        #else: print(f"STATE.SAMPLE_GEOMETRIES: could not understand {typ=}. Options are light/default/heavy")

        if not self.check_minimum():
            raise ValueError("State is not a minimum")
        if not hasattr(self.VNMs[0],"has_mode"):
            raise ValueError("VNMs do not have Eigenvectors. Please parse them")

        #if debug > 0: print(f"STATE.SAMPLE_GEOMETRIES: you selected {typ=}. Parameters are:")
        #if debug > 0: print(f"\t Temperature (Relevant for Sampling):  {temp}")
        #if debug > 0: print(f"\t Number of Rounds:                     {n_rounds}")
        #if debug > 0: print(f"\t Number of Selected Samples per Round: {n_samples_round}")

        q_min = np.zeros((len(self.VNMs)))
        current_geoms = [] 
        current_q_disp = [] 
        current_geoms.append(self.coord)
        current_q_disp.append(q_min)

        ## Main Loop
        for nr in range(n_rounds):
            q_fps, g_fps = [], [] 
            count = 0
            for c, q in zip(current_geoms, current_q_disp):
                count += 1
                geoms, q_disp, energies = geom_sampling_from_vnm(self.labels, c, self.VNMs, qini=q, T=temp, n_samples=n_samples_round, check_adjacencies=True, debug=debug)
                if debug > 0: print(f"STATE.SAMPLE_GEOMETRIES: Initial structure {count}/{len(current_geoms)} of round {nr+1}/{n_rounds} sampled {len(q_disp)} starting geometries")

                # Data for FPS
                if nr+1 < n_rounds:  ## Minimum (i.e. initial structure, with Q=0) is not added in the last round
                    q_fps.append(q_min)
                    g_fps.append(self.coord)
                q_fps.extend(q_disp)
                g_fps.extend(geoms)

            if len(q_fps) > n_selected:
                #Run FPS
                if debug > 0: print(f"STATE.SAMPLE_GEOMETRIES: Entering FPS selection with {len(q_fps)} geometries. Keeping {n_selected}")
                #idxs = furthest_point_sampling(q_fps, self.freqs_cm, n_selected, beta_distance)
                idxs = furthest_point_sampling(q_fps, n_selected, euclidean_q_distance)
                if debug > 0: print(f"STATE.SAMPLE_GEOMETRIES: FPS of round {nr+1}/{n_rounds} kept {len(idxs)} geometries, with indices:{idxs}")
            else:
                if debug > 0: print(f"STATE.SAMPLE_GEOMETRIES: Sampling of round {nr+1}/{n_rounds} failed to generate enough samples for FPS. Taking all available to the next round")
                idxs = list(range(len(q_disp)))

            # Prepares First Sequential Round 
            current_geoms = [] 
            current_q_disp = [] 
            for idx in idxs:
                current_geoms.append(g_fps[idx])
                current_q_disp.append(q_fps[idx])

        return current_q_disp, current_geoms

#######################################
#### Results associated with State ####
#######################################
    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  
            self.results[result.key] = result

    def set_energy(self, energy, units, overwrite: bool=True):
        self.add_result(data("energy",energy,units,"state.set_energy()"), overwrite=overwrite)

    def set_Helec(self, overwrite: bool=True, debug: int=0):
        assert "energy" in self.results
        if not hasattr(self,"ncomplex"): self.get_ncomplex(debug=debug)
        self.add_result(data("Helec",self.results["energy"].value/self.ncomplex,self.results["energy"].units,"state.set_Helec()"), overwrite=overwrite)

################################
#### Get Thermodynamic Data ####
################################
    def get_thermal_data(self, Trange: range=range(10,501,1), Helec=None, Selec=None, Hvib=None, Svib=None, Gtot=None, overwrite: bool=False, debug: int=0):
        from Scope.Thermal_Corrections import get_Selec, get_Hvib, get_Svib, get_Gibbs

        if not hasattr(self,"ncomplex"): self.get_ncomplex(debug=debug)
        if debug > 0: print(f"STATE.GET_THERMAL_DATA: found {self.ncomplex} complex molecules")

        if Hvib is None and Svib is None:
            assert hasattr(self,"isminimum"), f"I can't compute thermal data on this state. Missing VNMs"
    
        ############## Helec ##############
        if Helec is None:   ### One can provide specific values for Helec, Selec, Hvib, Svib and Gtot 
            if overwrite or not "Helec" in self.results.keys():
                self.add_result(data("Helec",self.results["energy"].value/self.ncomplex,self.results["energy"].units,"state.get_thermal_data()"), overwrite=overwrite)
        else: 
            if isinstance(Helec, data):
                if overwrite or not "Helec" in self.results.keys():
                    self.add_result(data("Helec",Helec.value,Helec.units,"enforced in state.get_thermal_data()"), overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Helec. It must be a DATA-class object")
        if debug > 0: print(f"Helec is {self.results['Helec']}")

        ############## Selec ##############
        if Selec is None:
            if overwrite or not "Selec" in self.results.keys():
                self.add_result(get_Selec(self._source.spin, outunits='au', nmol=self.ncomplex), overwrite=overwrite)
        else: 
            if isinstance(Selec, data):
                if overwrite or not "Selec" in self.results.keys():
                    self.add_result(data("Selec",Selec.value,Selec.units,"enforced in state.get_thermal_data()"), overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Selec. It must be a DATA-class object")
        if debug > 0: print(f"Selec is {self.results['Selec']}")

        ############## Hvib ##############
        if Hvib is None:
            if overwrite or not "Hvib" in self.results.keys():
                Hvib = collection("Hvib", "Temperature")
                for temp in Trange:
                    Hvib.add_data(get_Hvib(np.abs(self.freqs_cm), temp, freq_units='cm', outunits='au', nmol=self.ncomplex))
                self.add_result(Hvib, overwrite=overwrite)
        else: 
            if isinstance(Hvib, collection):
                if overwrite or not "Hvib" in self.results.keys():
                    self.add_result(Hvib, overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Hvib. It must be a COLLECTION-class object")
        if debug > 0: print(f"Hvib is {self.results['Hvib']}")

        ############## Svib ##############
        if Svib is None:
            if overwrite or not "Svib" in self.results.keys():
                Svib = collection("Svib", "Temperature")
                for temp in Trange:
                    Svib.add_data(get_Svib(np.abs(self.freqs_cm), temp, freq_units='cm', outunits='au', nmol=self.ncomplex))
                self.add_result(Svib, overwrite=overwrite)
        else: 
            if isinstance(Svib, collection):
                if overwrite or not "Svib" in self.results.keys():
                    self.add_result(Svib, overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Svib. It must be a COLLECTION-class object")
        if debug > 0: print(f"Svib is {self.results['Svib']}")

        ############## Gtot ##############
        if Gtot is None:
            if overwrite or not "Gtot" in self.results.keys():
                Gtot = collection("Gtot", "Temperature")
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
                    new_data = data(key, value, units, function)
                    new_data.add_property("temperature", temp, overwrite=overwrite)
                    Gtot.add_data(new_data)
                self.add_result(Gtot, overwrite=overwrite)
        else: 
            if isinstance(Gtot, collection):
                if overwrite or not "Gtot" in self.results.keys():
                    self.add_result(Gtot, overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Gtot. It must be a COLLECTION")
        if debug > 0: print(f"Gtot is {self.results['Gtot']}")

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
        if hasattr(self,"ncomplex"):       to_print += f' Number of Complexes   = {self.ncomplex}\n' 
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
