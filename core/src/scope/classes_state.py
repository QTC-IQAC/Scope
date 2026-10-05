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
    Represent a geometry- and result-bearing state of a source object.

    Attributes:
        object_type (str):              Object category (`"state"`).
        name (str):                     State name.
        _source (object):               Parent molecule, specie, or cell.
        results (dict):                 Registered results for the state.
        computations (list):            Computations linked to the state.
        labels (list):                  Atomic symbols for the geometry.
        coord (list):                   Cartesian coordinates.

    Methods:
        set_geometry():                 Store labels and coordinates.
        set_cell():                     Attach unit-cell metadata.
        get_molecules():                Build molecular fragments from the geometry.
        set_VNMs():                     Register vibrational normal modes.
        get_thermal_data():             Compute thermodynamic quantities.
    """
    def __init__(self, _source: object, name: str, debug: int=0):
        self.object_type    = "state"
        self.object_subtype = "state"
        self._source        = _source
        self.name           = name
        self.results        = dict()
        self.computations   = []

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
        if hasattr(self,"cell_vector"):
            self.frac_coord = cart2frac(self.coord, self.cell_vector)
        self.get_molecules(overwrite=True, debug=debug) ## Molecules require updated

    def set_geometry_from_molecules(self, overwrite: bool=False, debug: int=0):
        if not hasattr(self,"molecules"): self.get_molecules(overwrite=overwrite, debug=debug)
        self.labels     = []
        self.coord      = []
        indices         = []
        for mol in self.molecules:
            if not hasattr(mol,"atoms"): mol.set_atoms()
            if mol.check_parent("state", search_by="type"):
                mol_indices = mol.get_parent_indices("state", search_by="type")
            elif hasattr(mol,"cell_indices"):
                mol_indices = mol.cell_indices
            else:
                mol_indices = mol.indices
            for idx, at in enumerate(mol.atoms):
                self.labels.append(at.label)
                self.coord.append(at.coord)
                indices.append(mol_indices[idx])
        ## Below is to order the atoms as in the original cell, using the indices stored in the molecule object
        self.labels  = [x for _, x in sorted(zip(indices, self.labels), key=lambda pair: pair[0])]
        self.coord   = [x for _, x in sorted(zip(indices, self.coord), key=lambda pair: pair[0])]
        self.natoms  = len(self.labels)
        self.formula = labels2formula(self.labels)
        assert len(self.labels) == len(self.coord)
        if hasattr(self,"cell_vector"):
            self.frac_coord = cart2frac(self.coord, self.cell_vector)
         
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
        self.get_molecules(overwrite=True)

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
    def get_molecules(self, overwrite: bool=False, cov_factor: float=1.3, metal_factor: float=1.0, smart: bool=False, bond_margin: float=0.1, debug: int=0):
        from scope.classes_specie import Molecule

        # Overwrite
        if not overwrite and hasattr(self,"molecules"): 
            if debug > 0: print(f"STATE.GET_MOLECULES. Molecules already exist and default is overwrite=False")
            return self.molecules
        # Security
        if not hasattr(self,"labels") or not hasattr(self,"coord"): 
            if debug > 0: print(f"STATE.GET_MOLECULES. State labels and coordinates not found. Returning None")
            return None
        if len(self.labels) == 0 or len(self.coord) == 0: 
            if debug > 0: print(f"STATE.GET_MOLECULES. State labels and coordinates are empty. Returning None")
            return None

        # State connectivity must describe the current geometry, which may differ from the source after optimization.
        if debug > 0: print(f"STATE.GET_MOLECULES: Constructing connectivity with {cov_factor=} and {metal_factor=}")
        blocklist = split_species(self.labels, self.coord, cov_factor=cov_factor, metal_factor=metal_factor, smart=smart, bond_margin=bond_margin, debug=debug)
        self.molecules = [] 

        # Creates a molecule for each disconnected fragment
        for b in blocklist:
            if debug > 0: print(f"STATE.GET_MOLECULES: doing block={b}")
            mol_labels      = extract_from_list(b, self.labels, dimension=1)
            mol_coord       = extract_from_list(b, self.coord, dimension=1)
            if hasattr(self,"frac_coord"):      mol_frac_coord = extract_from_list(b, self.frac_coord, dimension=1)
            elif hasattr(self,"cell_vector"):   mol_frac_coord = cart2frac(mol_coord, self.cell_vector)
            else:                               mol_frac_coord = None
            # Creates Molecule Object
            newmolec    = Molecule(mol_labels, mol_coord, mol_frac_coord)
            # For debugging
            newmolec.origin = "state.get_molecules"
            # Adds State as parent of the molecule, with indices b
            newmolec.add_parent(self, indices=b, debug=debug)
            # Construct both regular and metal-only adjacency matrices from the State geometry.
            newmolec.get_adjmatrix(smart=smart, cov_factor=cov_factor, metal_factor=metal_factor, bond_margin=bond_margin, debug=debug)
            # Creates The atom objects with adjacencies
            newmolec.set_atoms(create_adjacencies=True, debug=debug)
            # The split_complex must be below the frac_coord, so they are carried on to the ligands    
            if newmolec.iscomplex: 
                if debug > 0: print(f"STATE.GET_MOLECULES: splitting complex")
                newmolec.split_complex(debug=debug)
            self.molecules.append(newmolec)
        return self.molecules

    ######
    def get_ncomplex(self, debug: int=0):
        ## Returns the number of Transition Metal Complexes (TMC) in the unit cell. 
        ## Gradually replacing it by Z (computed with get_z).
        if debug > 0: print(f"STATE.GET_NCOMPLEX checking fragmentation")
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"Found Fragmented molecules in the geometry of state: {self.name}"
         
        if debug > 0: print(f"STATE.GET_NCOMPLEX getting molecules")
        if not hasattr(self,"molecules"): self.get_molecules(debug=debug)
        self.ncomplex = 0
        for mol in self.molecules:
            if mol.iscomplex: self.ncomplex += 1
        if debug > 0: print(f"State.get_ncomplex {self.ncomplex} complexes found in state: {self.name}")
        return self.ncomplex

    ######
    def get_z(self, debug: int=0):
        ## Returns the number of stoichiometric units in the unit cell. 
        ## Basically, how many times the same stoichiometry unit is repeated in the cell
        from scope.operations.vecs_and_mats import gcd_list
        if debug > 0: print(f"STATE.GET_Z: checking fragmentation")
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"STATE.GET_Z: found Fragmented molecules in the geometry of self: {self.name}"

        if debug > 0: print(f"STATE.GET_Z: getting molecules")
        if not hasattr(self,"molecules"): self.get_molecules(debug=debug)
        if debug > 0: print(f"STATE.GET_Z: received {len(self.molecules)} molecules")

        unique = [] 
        occurrences = []
        for mol in self.molecules:
            found = False
            for uni in unique:
                if mol == uni: found = True
            if not found: 
                unique.append(mol)
                occurrences.append(self.get_occurrence(mol, debug=debug))
        if debug > 0: print(f"STATE.GET_Z: {occurrences=}")
        self._z = int(gcd_list(occurrences))
        return self._z 
    
    ######
    def get_occurrence(self, substructure: object, debug: int=0) -> int:
        """
        Count how many times a substructure appears in the state.

        Parameters:
            substructure (object):       Substructure to search for.
            debug (int):                 Verbosity level.

        Returns:
            int: Number of occurrences found.
        """
        ## Finds how many times a substructure appears in self
        occurrence = 0

        if debug > 0: print(f"STATE.GET_OCCURRENCE checking fragmentation")
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"STATE.GET_OCCURRENCE found fragmented molecules in the geometry of cell: {self.name}"
         
        if debug > 0: print(f"STATE.GET_OCCURRENCE getting molecules")
        if not hasattr(self,"molecules"): self.get_molecules(debug=debug)

        ## Case of Species inside self
        if hasattr(substructure,"object_type"):
            if substructure.object_type == 'specie':
                for mol in self.molecules:
                    if mol.__eq__(substructure, with_graph=True): occurrence += 1
        return occurrence

########################
#### Reconstruction ####
########################
    def reconstruct(self, cov_factor: float=1.3, metal_factor: float=1.0, debug: int=0):
        from scope.reconstruct import classify_fragments, fragments_reconstruct 
        if not self._source.object_type == "cell":    raise ValueError(f"STATE_RECONSTRUCT: state's source should by a CELL object") 
        if not hasattr(self,"cell_vector"):           raise ValueError(f"STATE_RECONSTRUCT: state should have a cell vector") 
        if not hasattr(self._source,"ref_molecules"): raise ValueError(f"STATE.RECONSTRUCT: state's source does not have a list of reference molecules"); return None
        from scope.read_write import HiddenPrints
        if debug > 0: print("STATE.RECONSTRUCT: reconstructing cell of state", self.name)
        with HiddenPrints():
            finished = False
            if not hasattr(self,"molecules"): self.get_molecules(debug=debug) 
            import itertools
            blocklist    = self.molecules.copy()
            ref_molecules = self._source.ref_molecules.copy()
            molecules, fraglist, Hlist = classify_fragments(blocklist, ref_molecules, debug=debug) 
            if len(fraglist) > 0 or len(Hlist) > 0: 
                molecules, finalmols, Warning = fragments_reconstruct(molecules,fraglist,Hlist,ref_molecules,self.cell_vector,cov_factor,metal_factor, debug=debug)
                molecules.extend(finalmols)
                self.molecules = molecules
                for mol in self.molecules:
                    if hasattr(mol,"cell_indices"):
                        mol_indices = mol.cell_indices
                    elif mol.check_parent("state", search_by="type"):
                        mol_indices = mol.get_parent_indices("state", search_by="type")
                    else:
                        mol_indices = mol.indices
                    mol.add_parent(self, mol_indices, debug=debug)
                    mol.set_atoms(atomlist=mol.atoms, debug=debug)
                    if mol.iscomplex: mol.split_complex(debug=debug)
                self.set_geometry_from_molecules()
                finished = True
        if debug > 0 and finished: print("STATE.RECONSTRUCT: state reconstructed succesfully")
        return self.molecules

    def check_fragmentation(self, reconstruct: bool = False, debug: int=0):
        ## If the source is a specie, then in principle there should only be one. So with molecules is enough to find
        if self._source.object_type == "specie":
            if not hasattr(self,"molecules"): self.get_molecules(debug=debug)
            if len(self.molecules) > 1: self.fragmented = True
            else:                       self.fragmented = False
            if debug > 0: print(f"STATE.CHECK_FRAGMENTATION: source type=specie. {self.fragmented=}")
            return self.fragmented
        ## If it is a unit cell, then we need a list of molecules that should in principle be there. This is ref_molecules
        ## If the cell is created by cell2mol, then this list is already stored in the .cell object
        elif self._source.object_type == "cell": 
            assert hasattr(self,"cell_vector")
            assert hasattr(self._source,"ref_molecules")
            if not hasattr(self,"molecules"): self.get_molecules(debug=debug)
            self.fragmented = False
            # First comparison with current molecules
            for mol in self.molecules:
                found = False
                for rmol in self._source.ref_molecules:
                    if mol.__eq__(rmol,with_graph=False): found = True # Graph cannot be used for rmol, since it doesn't have rdkit object
                if not found: self.fragmented = True
            # If there are fragments and user wants reconstruction, it tries to reconstruct and checks the new molecules
            if self.fragmented and reconstruct:
                new_molecules = self.reconstruct(debug=debug)
                self.fragmented = False
                for mol in new_molecules:
                    found = False
                    for rmol in self._source.ref_molecules:
                        if mol.__eq__(rmol,with_graph=False): found = True # Graph cannot be used for rmol, since it doesn't have rdkit object
                    if not found: self.fragmented = True
                if not self.fragmented: 
                    self.molecules = new_molecules
                    self.set_geometry_from_molecules(debug=debug)
        else: print(f"STATE.CHECK_FRAGMENTATION: Unknown source Type {self._source.object_type}")
        return self.fragmented

##############################
#### Connection with VNMs ####
##############################
    def check_minimum(self, debug: int=0):
        ## I think this function could be removed
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

        ## Checks if it is a minimum energy structure: all positive frequencies 
        if all(vnm.freq_cm >= 0.0 for vnm in self.VNMs): self.isminimum = True
        else:                                            self.isminimum = False

        ## Checks if it is a TS: only one very negative frequency 
        if sum(vnm.freq_cm < 0.0 for vnm in self.VNMs) == 1 and min(self.freqs_cm) < -50: self.is_ts = True
        else:                                                                             self.is_ts = False

        ## If it is not a minimum, evaluates if, at least, is close
        if not self.isminimum:
            self.num_neg_freqs = 0
            for vnm in self.VNMs: 
                if vnm.freq_cm < 0.0: self.num_neg_freqs += 1
            if self.num_neg_freqs <= 3 and VNMs[0].freq_cm > -50: self.almost_minimum = True
            else:                                                 self.almost_minimum = False

    def get_ir_spectrum(self, vmin=None, vmax=None, function: str='gaussian', sigma: float=10, debug: int=0):
        from scope.operations.vecs_and_mats import build_spectrum
        """
        Build a simulated IR spectrum from the stored VNMs.

        Parameters:
            vmin (float | None):         Lower bound of the frequency range.
            vmax (float | None):         Upper bound of the frequency range.
            function (str):              Broadening kernel.
            sigma (float):               Broadening width.
            debug (int):                 Verbosity level.

        Returns:
            tuple: Frequency grid and broadened intensity.
        """

        # Extract frequencies and intensities
        freqs       = np.array([v.freq_cm for v in self.VNMs])
        intensities = np.array([v.IR_int for  v in self.VNMs])

        # Define range
        if vmin is None:  vmin = -100
        if vmax is None:  vmax = freqs.max() + 100
        xrange = np.linspace(vmin, vmax, 2000)
        self.ir_spec_x, self.ir_spec_y = build_spectrum(xrange, freqs, intensities, function=function, sigma=sigma, debug=debug)
        return self.ir_spec_x, self.ir_spec_y

    def plot_ir_spectrum(self, vmin=None, vmax=None, function: str='gaussian', sigma: float=10, debug: int=0):
        import matplotlib.pyplot as plt
        if not hasattr(self,"ir_spec_x"): self.get_ir_spectrum(vmin=vmin, vmax=vmax, function=function, sigma=sigma, debug=debug)

        x, y = self.ir_spec_x, self.ir_spec_y
        spectrum = np.zeros_like(x)

        freqs       = np.array([v.freq_cm for v in self.VNMs])
        intensities = np.array([v.IR_int for  v in self.VNMs])
        
        # Plot
        plt.figure(figsize=(8,4))
        plt.plot(x, y, 'k-', lw=1.5)
        plt.fill_between(x, 0, spectrum, color="grey", alpha=0.4)
        # Add sticks
        for f, I in zip(freqs, intensities):
            plt.vlines(f, 0, I, color="r", lw=1, linestyle="--")
        plt.xlabel("Wavenumber (cm$^{-1}$)")
        plt.ylabel("Intensity (a.u.)")
        plt.title("Simulated IR Spectrum")
        plt.tight_layout()
        plt.show()

########################################
#### Connection with Excited States ####
########################################
    def set_exc_states(self, exc_states, debug: int=0):
        self.exc_states = exc_states
        return self.exc_states
    
    def shift_exc_states_wl(self, shift: float, debug: int=0):
        # Shifts the wavelength of the Excited states, applied in nm
        for es in self.exc_states:
            es.shift_wavelength(shift, debug=debug)

    def shift_exc_states_energy(self, shift: float, debug: int=0):
        # Shifts the wavelength of the Excited states, applied in eV
        for es in self.exc_states:
            es.shift_energy(shift, debug=debug)

    def restore_exc_states(self): 
        # Reverts any changes in wavelength or energy
        for es in self.exc_states:
            es.restore()

    def get_abs_spectrum(self, lmin: float=200, lmax: float=1000, function: str='gaussian', sigma: float=0.2, as_cross_section: bool=False, debug: int=0):
        """
        Build an absorption spectrum from the stored excited states.

        Parameters:
            lmin (float):                Lower wavelength bound in nm.
            lmax (float):                Upper wavelength bound in nm.
            function (str):              Broadening kernel.
            sigma (float):               Broadening width.
            as_cross_section (bool):     Normalize as a cross-section-like spectrum.
            debug (int):                 Verbosity level.

        Returns:
            tuple: Wavelength grid and intensity values.
        """
        from scope.operations.vecs_and_mats import build_spectrum
        from scope import constants

        # Check if TDDFT data exists.
        if not hasattr(self, 'exc_states'): raise ValueError('AZO.GET_ABS_SPECTRUM: [WARNING] No TDDFT data found in this state')

        # Collects Values
        energies = [es.energy for es in self.exc_states] # in eV
        fosc     = [es.fosc for es in self.exc_states]   # oscillator strength   
        if debug > 0: print(f'STATE_AZO.GET_ABS_SPECTRUM: energies {energies}')
        if debug > 0: print(f'STATE_AZO.GET_ABS_SPECTRUM: osc. strengths {fosc}')

        ## Convert desired range in nm (lrange) to energies (erange)
        lrange = np.linspace(lmin, lmax, lmax-lmin)
        erange = constants.hc/lrange[::-1]
        if debug > 0: print(f'STATE_AZO.GET_ABS_SPECTRUM: erange {np.min(erange):6.4f}-{np.max(erange):6.4f}')

        # Builds the spectrum from discrete values, using Gaussian broadening
        # If normalize = True, the spectrum is normalized and the resulting 'y' has units of energy-1.
        # If normalize = False, the resulting 'y' is just a deconvoluted spectrum in fosc units
        normalize = True if as_cross_section else False
        x, y = build_spectrum(erange, energies, fosc, function=function, sigma=sigma, normalize=normalize, debug=debug)

        self.abs_spec_x = constants.hc/x[::-1]  # Converts the result to a range of nm values   
        self.abs_spec_y = y[::-1]
        return self.abs_spec_x, self.abs_spec_y

    def get_cross_section(self, lmin: float=200, lmax: float=1000, function: str='gaussian', sigma: float=0.2, debug: int=0):
        """
        Build an absorption cross section from the stored excited states.

        Parameters:
            lmin (float):                Lower wavelength bound in nm.
            lmax (float):                Upper wavelength bound in nm.
            function (str):              Broadening kernel.
            sigma (float):               Broadening width.
            debug (int):                 Verbosity level.

        Returns:
            tuple: Wavelength grid and cross-section values.
        """
        from scope.operations.vecs_and_mats import build_spectrum
        from scope import constants

        # Check if TDDFT data exists.
        if not hasattr(self, 'exc_states'): raise ValueError('AZO.GET_CROSS_SECTION: [WARNING] No TDDFT data found in this state')

        # Collects Values
        energies = [es.energy for es in self.exc_states] # in eV
        fosc     = [es.fosc for es in self.exc_states]   # oscillator strength   
        if debug > 0: print(f'STATE_AZO.GET_CROSS_SECTION: energies {energies}')
        if debug > 0: print(f'STATE_AZO.GET_CROSS_SECTION: osc. strengths {fosc}')

        ## Convert desired range in nm (lrange) to energies (erange)
        lrange = np.linspace(lmin, lmax, lmax-lmin)
        erange = constants.hc/lrange[::-1]
        if debug > 0: print(f'STATE_AZO.GET_CROSS_SECTION: erange {np.min(erange):6.4f}-{np.max(erange):6.4f}')

        # Builds the spectrum from discrete values, using Gaussian broadening
        x, y = build_spectrum(erange, energies, fosc, function=function, sigma=sigma, normalize=True, debug=debug)
                
        # Here, we convert from energy-1 in eV, to absorption cross section in m2
        K = (constants.planck_Js * constants.elem_charge) / (4 * constants.epsilon_0 * constants.speed_light * constants.electron_mass)
        y = y * K

        self.cross_sec_x = constants.hc/x[::-1]  # Converts the result to a range of nm values   
        self.cross_sec_y = y[::-1]
        return self.cross_sec_x, self.cross_sec_y

    def plot_abs_spectrum(self, lmin: float=200, lmax: float=1000, function: str='gaussian', sigma: float=0.2, debug: int=0):
        import matplotlib.pyplot as plt
        x, y = self.get_abs_spectrum(lmin=lmin, lmax=lmax, function=function, sigma=sigma, debug=debug)
        fig, ax = plt.subplots(figsize=(3, 2), dpi=200)
        ax.plot(x, y, color='black')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel(r'osc. strength (a.u.)')  
        ax.set_xlim(lmin, lmax)
        plt.show()
    
##################################
#### Connection with Workflow ####
##################################
    def find_computation(self, job_name: str='', step: int=1, run_number: int=1, debug: int=0):
        for idx, comp in enumerate(self.computations):
            if comp._job.name == job_name and comp.step == step and comp.run_number == run_number: this_comp = comp; return True, this_comp
        return False, None

    def add_computation(self, computation: object, debug: int=0):
        found, comp = self.find_computation(computation._job.name, computation.step, computation.run_number)
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
        Sample geometries around the current state using vibrational modes.

        Parameters:
            ngeoms (int):                Number of geometries to keep.
            n_aux_geoms (int):           Number of trial geometries per FPS round.
            n_fps_rounds (int):          Number of furthest-point-sampling rounds.
            temp (float):                Sampling temperature.
            sigma_damp_factor (float):   Damping factor for mode amplitudes.
            freq_bottom_limit (float):   Lower frequency cutoff.
            debug (int):                 Verbosity level.

        Returns:
            tuple: Selected Q displacements and Cartesian geometries.
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

    def remove_result(self, key: str):
        return self.results.pop(key, None)

    def set_energy(self, energy, units, overwrite: bool=True):
        self.add_result(Data("energy",energy,units,"state.set_energy()"), overwrite=overwrite)

    def set_Helec(self, overwrite: bool=True, debug: int=0):
        assert "energy" in self.results
        if not hasattr(self,"z"): self.get_z(debug=debug)
        self.add_result(Data("Helec",self.results["energy"].value/self.z,self.results["energy"].units,"state.set_Helec()"), overwrite=overwrite)
        
################################
#### Get Thermodynamic Data ####
################################
    def get_thermal_data(self, temp: float=298.15, Helec=None, Selec=None, Hvib=None, Svib=None, Gtot=None, overwrite: bool=False, vib_options: dict=None, debug: int=0):
        """Computes and Handles Thermochemistry results

        Parameters:
            temp:                       Temperature (in K) as a float or an iterable of temperatures.
            Helec, Selec:               Optional enforced Data objects.
            Hvib, Svib, Gtot:           Optional enforced Collections with Temperature as variable.
            overwrite (bool):           Force replacement, including unchanged settings.
            vib_options (dict):         model ('HO' or 'QRRHO'), FR_cutoff (cm-1), FR_alpha, and imaginary treatment.
            debug (int):                Verbosity; child functions receive debug - 1.

        The option "imaginary" in vib_options sets the policy towards imaginary frequencies. 
        It applies when computing both Hvib and Svib:
        - 'ignore' (default): excludes negative frequencies from vibrational terms;
          the usual treatment for a transition-state reaction coordinate.
        - 'absolute': uses their absolute values as real, positive frequencies;
          a diagnostic comparison, not the usual transition-state treatment.
        - 'raise': stops with ValueError upon encountering a negative frequency;
          useful for checking intended minima, but expected to reject transition
          states with imaginary modes. It does not simply skip the offending mode.

        Settings changes replace affected collections automatically. Gibbs energies
        are refreshed when contributing results are replaced. Use overwrite=True
        after changing source energies or frequencies without replacing results.
        """
        from scope.thermodynamics import get_Selec, get_Hvib, get_Svib, get_Gibbs, normalize_vib_options
        Svib_options  = normalize_vib_options(vib_options)
        Hvib_options  = {'imaginary': Svib_options['imaginary']}
        Svib_settings = {'model': Svib_options['model']}
        if Svib_options['model'] == 'QRRHO':
            Svib_settings.update(fr_cutoff=Svib_options['FR_cutoff'], fr_alpha=Svib_options['FR_alpha'])
        Svib_settings['imaginary'] = Svib_options['imaginary']
        child_debug = max(debug - 1, 0)

        # 0) Checks inputs before changing stored results
        if isinstance(temp, (int, float)): temperatures = [temp]
        else:
            try: temperatures = list(temp)
            except TypeError as exc: raise TypeError("STATE.GET_THERMAL_DATA: temp must be numeric or an iterable of temperatures") from exc
        if not temperatures: raise ValueError("STATE.GET_THERMAL_DATA: temp cannot be empty")
        if not all(isinstance(t, (int, float)) for t in temperatures): raise TypeError("STATE.GET_THERMAL_DATA: all temperatures must be numeric")
        if (Hvib is None or Svib is None) and not hasattr(self, "VNMs"): raise ValueError("STATE.GET_THERMAL_DATA: Missing VNMs")
        if 'energy' not in self.results or self.results['energy'] is None: raise ValueError("STATE.GET_THERMAL_DATA: Missing State energy")
        for key, supplied in [('Helec', Helec), ('Selec', Selec)]:
            if supplied is not None and not isinstance(supplied, Data): raise TypeError(f"STATE.GET_THERMAL_DATA: Provided {key} must be a Data object")
        for key, supplied in [('Hvib', Hvib), ('Svib', Svib), ('Gtot', Gtot)]:
            if supplied is None: continue
            if not isinstance(supplied, Collection): raise TypeError(f"STATE.GET_THERMAL_DATA: Provided {key} must be a Collection")
            if supplied.variable.lower() != 'temperature': raise ValueError(f"STATE.GET_THERMAL_DATA: Provided {key} must scan temperature")
            for temperature in temperatures:
                if supplied.find_value_with_property('temperature', temperature) is None: raise ValueError(f"STATE.GET_THERMAL_DATA: Provided {key} lacks temperature {temperature}")
        if not hasattr(self, "z"): self.get_z(debug=child_debug)

        # 1) Stores electronic contributions
        if overwrite or 'Helec' not in self.results:
            if Helec is None: Helec = Data('Helec', self.results['energy'].value/self.z, self.results['energy'].units, 'state.get_thermal_data()')
            else: Helec = Data('Helec', Helec.value, Helec.units, 'enforced in state.get_thermal_data()')
            self.add_result(Helec, overwrite=True)
        if overwrite or 'Selec' not in self.results:
            if Selec is None: Selec = get_Selec(self.spin_multiplicity, outunits='au', nmol=self.z)
            else: Selec = Data('Selec', Selec.value, Selec.units, 'enforced in state.get_thermal_data()')
            self.add_result(Selec, overwrite=True)

        # 2) Reuses compatible vibrational results and fills missing temperatures
        for key, supplied, settings, options, calculate in [('Hvib', Hvib, Hvib_options, Hvib_options, get_Hvib), ('Svib', Svib, Svib_settings, Svib_options, get_Svib)]:
            stored = self.results.get(key)
            if supplied is not None:
                if overwrite or stored is None: self.add_result(supplied, overwrite=True)
                continue
            replace = overwrite or not isinstance(stored, Collection) or not stored.check_settings(settings)
            if replace:
                stored = Collection(key, 'temperature')
                if debug > 0: print(f"STATE.GET_THERMAL_DATA: Computing {key} with {settings}")
            for temperature in temperatures:
                if stored.find_value_with_property('temperature', temperature) is None:
                    stored.add_data(calculate(self.freqs_cm, temperature, freq_units='cm', outunits='au', nmol=self.z, **options, debug=child_debug))
            self.add_result(stored, overwrite=True)

        # 3) Refreshes Gibbs energies when settings or contributing results change
        if Gtot is not None:
            if overwrite or 'Gtot' not in self.results: self.add_result(Gtot, overwrite=True)
        else:
            Helec = self.results['Helec']
            Selec = self.results['Selec']
            Hvib  = self.results['Hvib']
            Svib  = self.results['Svib']
            Gtot_settings = {setting: getattr(Svib, setting) for setting in getattr(Svib, 'settings', [])}
            stored_Gtot = self.results.get('Gtot')
            replace = overwrite or not isinstance(stored_Gtot, Collection) or not stored_Gtot.check_settings(Gtot_settings)
            if not replace:
                for data in stored_Gtot.datas:
                    inputs = (Helec, Selec, Hvib.find_value_with_property('temperature', data.temperature), Svib.find_value_with_property('temperature', data.temperature))
                    previous = getattr(data, '_thermal_inputs', ())
                    if len(previous) != len(inputs) or any(old is not current for old, current in zip(previous, inputs)):
                        replace = True
                        break
            if replace: stored_Gtot = Collection('Gtot', 'temperature')
            for temperature in temperatures:
                if stored_Gtot.find_value_with_property('temperature', temperature) is not None: continue
                Hvib_i = Hvib.find_value_with_property('temperature', temperature)
                Svib_i = Svib.find_value_with_property('temperature', temperature)
                assert Helec.units == Selec.units == Hvib_i.units == Svib_i.units, f"{Helec.units=}, {Selec.units=}, {Hvib_i.units=}, {Svib_i.units=}"
                value = get_Gibbs(Helec.value, Hvib_i.value, Selec.value, Svib_i.value, temperature)
                data = Data('Gtot', value, Helec.units, 'state.get_thermal_data()')
                data.add_property('temperature', temperature)
                for setting, setting_value in Gtot_settings.items(): data.add_setting(setting, setting_value)
                data.vib_options = Svib_i.vib_options.copy() if hasattr(Svib_i, 'vib_options') else None
                data._thermal_inputs = (Helec, Selec, Hvib_i, Svib_i)
                stored_Gtot.add_data(data)
            self.add_result(stored_Gtot, overwrite=True)
        if debug > 0: print(f"STATE.GET_THERMAL_DATA: Thermal results available at {temperatures} K")

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
        if self._source.object_type != 'cell': 
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
    def __repr__(self, indirect: bool=False) -> None:
        to_print = ''
        if not indirect: to_print += f'--------------------------------\n'
        if not indirect: to_print += f'------ SCOPE STATE Object ------\n'  
        if not indirect: to_print += f'--------------------------------\n'
        to_print += f' Name                  = {self.name}\n'
        if hasattr(self._source,"name"):        to_print += f' Source Name           = {self._source.name}\n'
        if hasattr(self._source,"object_type"): to_print += f' Source Type           = {self._source.object_type}\n'
        if hasattr(self,"labels"):         to_print += f' Labels                = {self.labels[0]}...\n'
        if hasattr(self,"coord"):          to_print += f' Coord                 = {self.coord[0]}...\n'
        if hasattr(self,"z"):              to_print += f' Number of Units (Z)   = {self.z}\n' 
        to_print += f'\n' 
        if hasattr(self,"VNMs"):           to_print += f' Has VNMs              = YES\n'
        if hasattr(self,"isminimum"):      to_print += f' Is Minimum            = {self.isminimum}\n'
        if hasattr(self,"almost_minimum"): to_print += f' Is Almost a Minimum   = {self.almost_minimum}\n'
        if hasattr(self,"is_ts"):          to_print += f' Is a Transition State = {self.is_ts}\n'
        if hasattr(self,"freqs_cm"):       to_print += f' First Freq (cm-1)     = {self.freqs_cm[0]}\n'
        to_print += f'\n' 

        if hasattr(self,"exc_states"):     to_print += f' Has Excited States    = YES\n'

        if hasattr(self,"molecules"):  
            to_print += f' Num of Molecules:     = {len(self.molecules)}\n'
            to_print += f' With Formulae:                               \n'
            for idx, m in enumerate(self.molecules):
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
