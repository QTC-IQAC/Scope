import numpy as np
from Scope.Adapted_from_cell2mol import * 
from Scope.Other import get_metal_idxs, get_dist
from Scope.Unit_cell_tools import * 
#from Scope.Reconstruct import * #error when loading here 

################################
####  BASIS FOR CELL2MOL 2  ####
################################
class specie(object):
    def __init__(self, labels: list, coord: list, frac_coord: list=None, radii: list=None) -> None:

       # Sanity Checks
        assert len(labels) == len(coord)
        if frac_coord is not None:        
            assert len(coord) == len(frac_coord)
            self.frac_coord = frac_coord
            
        # Optional Information
        if radii   is not None: self.radii   = radii
        else:                   self.radii   = get_radii(labels)

        self.type      = "specie"
        self.version   = "2.0"
        self.labels    = labels
        self.coord     = coord
        self.formula   = labels2formula(labels)
        self.eleccount = labels2electrons(labels)
        self.natoms    = len(labels)
        self.iscomplex = any((elemdatabase.elementblock[l] == "d") or (elemdatabase.elementblock[l] == "f") for l in self.labels)

        self.parents           = []
        self.parents_indices   = []
        self.indices = [*range(0,self.natoms,1)]

        ## Defaults
        self.cov_factor   = 1.3
        self.metal_factor = 1.0

    ############
    def add_parent(self, parent: object, indices: list, overwrite: bool=True, debug: int=0):
        ## associates a parent specie to self. The atom indices of self in parent are given in "indices"
        ## if parent of the same subtype already in self.parent then it is overwritten
        ## this is to avoid having a substructure (e.g. a ligand) in more than one superstructure (e.g. a molecule) 

        # 1st-evaluates parent
        append = True
        for idx, p in enumerate(self.parents):
            if p.subtype == parent.subtype:
                if overwrite: 
                    self.parents[idx]         = parent
                    self.parents_indices[idx] = indices
                append = False
        if append: 
            self.parents.append(parent)
            self.parents_indices.append(indices)
            if debug > 0: print(f"SPECIE.ADD_PARENT: added parent with subtype={parent.subtype}. Indices are {indices}")

        # 2nd-evaluates parents of parent
        if hasattr(parent,"parents"):
            for jdx, p2 in enumerate(parent.parents):
                append = True
                new_indices = extract_from_list(indices, parent.get_parent_indices(p2.subtype), dimension=1)
                for idx, p in enumerate(self.parents):
                    if p.subtype == p2.subtype:
                        if overwrite: 
                            self.parents[idx]         = p2
                            #self.parents_indices[idx] = parent.get_parent_indices(p2.subtype)
                            self.parents_indices[idx] = new_indices
                        append = False
                if append: 
                    self.parents.append(p2)
                    #self.parents_indices.append(parent.get_parent_indices(p2.subtype))
                    self.parents_indices.append(new_indices)
                    if debug > 0: print(f"SPECIE.ADD_PARENT: added parent of parent with subtype {p2.subtype}. Indices are {new_indices}")
                    #if debug > 0: print(f"SPECIE.ADD_PARENT: added parent of parent with {p2.subtype=}. Indices are {parent.get_parent_indices(p2.subtype)}")

    ############
    def check_parent(self, subtype: str):
        ## checks if parent of a given subtype exists
        for p in self.parents:
            if p.subtype == subtype: return True
        return False

    ############
    def get_parent(self, subtype: str):
        ## retrieves parent of a given subtype 
        for p in self.parents:
            if hasattr(p, "subtype"):  
                if p.subtype.lower() == subtype.lower(): return p
            else:  
                print(f"Warning. Parent with type: {p.type} does not have subtype")
        return None

    ############
    def get_parent_indices(self, subtype: str):
        ## retrieves parent indices of a given subtype 
        for idx, p in enumerate(self.parents):
            if hasattr(p, "subtype"):  
                if p.subtype.lower() == subtype.lower(): return self.parents_indices[idx]
        return None

    ############
    def get_centroid(self):
        self.centroid = compute_centroid(self.coord)
        if hasattr(self,"frac_coord"): self.frac_centroid = compute_centroid(self.frac_coord)
        return self.centroid

    ############
    def set_fractional_coord(self, frac_coord: list, debug: int=0) -> None:
        if debug > 0: print(f"SPECIE.SET_FRACTIONAL_COORD: length comparison: {len(frac_coord)} vs {len(self.coord)}")
        if debug > 0: print(f"SPECIE.SET_FRACTIONAL_COORD: length comparison: {np.shape(frac_coord)} vs {np.shape(self.coord)}")
        if debug > 0: print(f"SPECIE.SET_FRACTIONAL_COORD: 1st items: {frac_coord[0]} vs {self.coord[0]}")
        #assert len(frac_coord) == len(self.coord)  ## Removed 'cos it was complaining despite both vectors having the same length
        self.frac_coord = frac_coord 

    ############
    def get_fractional_coord(self, cell_vector=None, debug: int=0) -> None:
        from Scope.Reconstruct import cart2frac

        # If cell_vector is provided, then its easy
        if cell_vector is not None:
            self.frac_coord = cart2frac(self.coord, cell_vector)
            assert len(self.frac_coord) == len(self.coord)
        # Otherwise, looks for a cell parent, and gets the cell_vector from there
        elif self.check_parent("cell"):
            par = self.get_parent("cell")
            if hasattr(par,"cell_vector"): 
                self.frac_coord = cart2frac(self.coord, par.cell_vector)
                assert len(self.frac_coord) == len(self.coord)
            else:  print("SPECIE.GET_FRACTIONAL_COORD. Parent cell is missing the cell vector"); return None
        else:  print("SPECIE.GET_FRACTIONAL_COORD. Please provide the cell_vector"); return None
        return self.frac_coord

    ############
    def get_atomic_numbers(self):
        if not hasattr(self,"atoms"): self.set_atoms()
        self.atnums = []
        for at in self.atoms:
            self.atnums.append(at.atnum)
        return self.atnums

    ############
    def set_element_count(self, heavy_only: bool=False):
        self.element_count = get_element_count(self.labels, heavy_only=heavy_only)
        return self.element_count

    ############
    def set_adj_types(self):
        if not hasattr(self,"adjmat"): self.get_adjmatrix()
        self.adj_types = get_adjacency_types(self.labels, self.adjmat)
        return self.adj_types

    ############
    # Replaces by set_factors
    #def set_adjacency_parameters(self, cov_factor: float, metal_factor: float) -> None:
    #    # Stores the covalentradii factor and metal factor that were used to generate the molecule
    #    self.cov_factor   = cov_factor
    #    self.metal_factor = metal_factor

    ############
    def reset_charge(self):
        if hasattr(self,"totcharge"):      delattr(self,"totcharge")
        if hasattr(self,"atomic_charges"): delattr(self,"atomic")
        if hasattr(self,"smiles"):         delattr(self,"smiles")
        if hasattr(self,"rdkit_obj"):      delattr(self,"rdkit_obj")
        if hasattr(self,"poscharges"):     delattr(self,"poscharges") 
        for a in self.atoms:               
            a.reset_charge() 

    ############
    def set_charges(self, totcharge: int=None, atomic_charges: list=None) -> None:
        ## Sets total charge  
        if totcharge is not None:                              self.totcharge = int(totcharge)
        elif totcharge is None and atomic_charges is not None: self.totcharge = int(np.sum(atomic_charges))
        elif totcharge is None and atomic_charges is None:     self.totcharge = "Unknown" 
        ## Sets atomic charges
        if atomic_charges is not None:
            self.atomic_charges = atomic_charges
            if not hasattr(self,"atoms"): self.set_atoms()
            for idx, a in enumerate(self.atoms):
                a.set_charge(self.atomic_charges[idx])

        ############
    def set_atoms(self, atomlist=None, create_adjacencies: bool=False, debug: int=0):
        ## If the atom objects already exist, and you want to set them in self from a different specie
        if atomlist is not None: 
            if debug > 0: print(f"SPECIE.SET_ATOMS: received {atomlist=}")
            self.atoms = atomlist.copy()
            for idx, at in enumerate(self.atoms):
                at.add_parent(self, index=idx)
                if debug > 0: print(f"SPECIE.SET_ATOMS: set parent {self.subtype} to atom, with index={idx}")

        ## If not, that is, if the atom objects must be created from scratch....
        else: 
            self.atoms = []
            for idx, l in enumerate(self.labels):
                if debug > 0: print(f"SPECIE.SET_ATOMS: creating atom for label {l}")
                ## For each l in labels, creates an atom class object.
                ismetal = elemdatabase.elementblock[l] == "d" or elemdatabase.elementblock[l] == "f"
                # non transition metals
                if len(get_non_transition_metal_idxs([l])) > 0: ismetal = True
                if ismetal and debug > 0: print(f"SPECIE.SET_ATOMS: {l}")
                if debug > 0:             print(f"SPECIE.SET_ATOMS: {ismetal=}")

                # Prepares fractional coordiantes
                if not hasattr(self,"frac_coord"): self.get_fractional_coord(debug=debug)
                if hasattr(self,"frac_coord"):
                    if ismetal: newatom = metal(l, self.coord[idx], self.frac_coord[idx], radii=self.radii[idx])
                    else:       newatom =  atom(l, self.coord[idx], self.frac_coord[idx],radii=self.radii[idx])
                else :
                    if ismetal: newatom = metal(l, self.coord[idx], radii=self.radii[idx])
                    else:       newatom =  atom(l, self.coord[idx], radii=self.radii[idx])
                if debug > 0: print(f"SPECIE.SET_ATOMS: added atom to specie: {self.formula}")
                newatom.add_parent(self, index=idx)
                self.atoms.append(newatom)
        
        if create_adjacencies:
            if not hasattr(self,"adjmat"):  self.get_adjmatrix()
            if not hasattr(self,"madjmat"): self.get_metal_adjmatrix()
            if self.adjmat is not None and self.madjmat is not None: 
                for idx, at in enumerate(self.atoms): 
                    at.set_adjacencies(self.adjmat[idx],self.madjmat[idx],self.adjnum[idx],self.madjnum[idx])

    #######################################################
    def inherit_adjmatrix(self, parent_subtype: str, debug: int=0):
        exists  = self.check_parent(parent_subtype)
        if not exists: 
            print(f"SPECIE.INHERIT. {parent_subtype=} does not exist")
            return None
        parent  = self.get_parent(parent_subtype)
        indices = self.get_parent_indices(parent_subtype)
        if not hasattr(parent,"madjnum"): 
            print(f"SPECIE.INHERIT. {parent_subtype=} does not have madjnum")
            return None 
        self.madjmat = np.stack(extract_from_list(indices, parent.madjmat, dimension=2), axis=0)
        self.madjnum = np.stack(extract_from_list(indices, parent.madjnum, dimension=1), axis=0)
        self.adjmat  = np.stack(extract_from_list(indices, parent.adjmat, dimension=2), axis=0)
        self.adjnum  = np.stack(extract_from_list(indices, parent.adjnum, dimension=1), axis=0)

    #######################################################
    def set_factors(self, cov_factor: float=1.3, metal_factor: float=1.0) -> None:
        self.cov_factor   = cov_factor
        self.metal_factor = metal_factor
        if hasattr(self,"atoms"):
            for at in self.atoms:
                at.set_factors(cov_factor=self.cov_factor, metal_factor=self.metal_factor)

    #######################################################
    def get_adjmatrix(self):
        isgood, adjmat, adjnum = get_adjmatrix(self.labels, self.coord, self.cov_factor, self.metal_factor, radii=self.radii)
        if isgood:
            self.adjmat = adjmat
            self.adjnum = adjnum
        else:
            self.adjmat = None
            self.adjnum = None
        return self.adjmat, self.adjnum

    def get_metal_adjmatrix(self):
        isgood, madjmat, madjnum = get_adjmatrix(self.labels, self.coord, self.cov_factor, self.metal_factor, radii=self.radii, metal_only=True)
        if isgood:
            self.madjmat = madjmat
            self.madjnum = madjnum
        else:
            self.madjmat = None
            self.madjnum = None
        return self.madjmat, self.madjnum

    def get_occurrence(self, substructure: object) -> int:
        occurrence = 0
        ## Ligands in Complexes or Groups in Ligands
        done = False
        if hasattr(substructure,"subtype") and hasattr(self,"subtype"):
            if substructure.subtype == 'ligand' and self.subtype == 'molecule':
                if not hasattr(self,"ligands"): self.split_complex()
                if self.ligands is not None:
                    for l in self.ligands:
                        issame = compare_species(substructure, l, debug=1)
                        if issame: occurrence += 1
                    done = True 
            elif substructure.subtype == 'group' and self.subtype == 'ligand':
                if not hasattr(self,"ligands"): self.split_complex()
                if self.ligands is not None:
                    for l in self.ligands:
                        if not hasattr(l,"groups"): self.split_ligand()
                        for g in l.groups:
                            issame = compare_species(substructure, g, debug=1)
                            if issame: occurrence += 1
                done = True 
        ## Atoms in Species
        if not done:
            if substructure.type == 'atom' and self.type == 'specie':
                if not hasattr(self,"atoms"): self.set_atoms()
                for at in self.atoms:
                    issame = compare_atoms(substructure, at)
                    if issame: occurrence += 1
        return occurrence

    def magnetism(self, spin: int) -> None:
        self.spin = spin

    def check_fragmentation(self, cov_factor: float=1.3, metal_factor: float=None, debug: int=0):
        blocklist = split_species(self.labels, self.coord, cov_factor=cov_factor)
        if len(blocklist) > 1: self.isfragmented = True
        else:                  self.isfragmented = False
        return self.isfragmented

    def print_xyz(self):
        print(self.natoms)
        print("")
        for idx, l in enumerate(self.labels):
            print("%s  %.6f  %.6f  %.6f" % (l, self.coord[idx][0], self.coord[idx][1], self.coord[idx][2]))

    def view(self, show_indices: bool=False, size: str='default'):
        import plotly.graph_objects as go
        from Scope.Read_Write import set_scene
        from Scope.Elementdata import ElementData  
        elemdatabase = ElementData()

        ### Adjusts size
        if   size.lower() == 'default': width=600; height=600;   marker_size=8;  text_size=9
        elif size.lower() == 'small':   width=400; height=400;   marker_size=6;  text_size=7 
        elif size.lower() == 'large':   width=800; height=800;   marker_size=10; text_size=12
        elif size.lower() == 'ultra':   width=1000; height=1000; marker_size=11; text_size=13

        if not hasattr(self,"adjmat"): self.get_adjmatrix()
        fig             = go.Figure()

        # Gather Data
        positions       = np.array(self.coord)
        symbols         = self.labels
        adjacencies     = self.adjmat
        atom_indices    = list(range(len(self.labels)))

        # Gets bonds from adjacency matrix
        indices = np.argwhere(self.adjmat >0)
        unique_bonds    = set()

        for i in indices:
            unique_bonds.add(tuple(i))

        # Plot atoms as markers
        fig.add_trace(go.Scatter3d(
            x           = positions[:, 0],
            y           = positions[:, 1],
            z           = positions[:, 2],
            mode        ='markers',
            marker      = dict(
                size        = marker_size,
                color       = [elemdatabase.cpk_colors[l] for l in symbols],
                line        = dict(color='black', width=1),
            ),
            hoverinfo   = 'text',
            text        = symbols,
            showlegend  = False
        ))

        if show_indices:
            fig.add_trace(go.Scatter3d(
                x           = positions[:, 0],
                y           = positions[:, 1],
                z           = positions[:, 2],
                mode        = 'text',
                text        = [str(i) for i in range(len(positions))],
                textfont    = dict(color='black', size=text_size),
                hoverinfo   = 'none',
                showlegend  = False
            ))
        else:            
            fig.add_trace(go.Scatter3d(
                x           = positions[:, 0],
                y           = positions[:, 1],
                z           = positions[:, 2],
                mode        = 'text',
                text        = self.labels,
                textfont    = dict(color='black', size=text_size),
                hoverinfo   = 'none',
                showlegend  = False
            ))

        ## Plot bonds as lines and calculate midpoints
        #midpoints   = []
        bond_pairs  = []

        for i, j in unique_bonds:
            # Add bond trace
            fig.add_trace(go.Scatter3d(
                x           = [positions[i, 0], positions[j, 0]],
                y           = [positions[i, 1], positions[j, 1]],
                z           = [positions[i, 2], positions[j, 2]],
                mode        = 'lines',
                line        = dict(color='gray', width=5),
                hoverinfo   = 'none',
                showlegend  = False
            ))

        #    # Calculate midpoints
        #    midpoint = (positions[i] + positions[j]) / 2
        #    midpoints.append(midpoint)
            bond_pairs.append((i, j))

        #midpoints = np.array(midpoints)

        set_scene(fig, np.array(self.coord), width=width, height=height)
        fig.show()


    ## To be implemented
    def __add__(self, other):
        if not isinstance(other, type(self)): return self
        return self

    ############
    def __repr__(self, indirect: bool=False):
        to_print = ""
        if not indirect: to_print  += f'------------- SCOPE SPECIE Object --------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        if hasattr(self,'subtype'): to_print += f' Sub-Type              = {self.subtype}\n'
        to_print += f' Number of Atoms       = {self.natoms}\n'
        to_print += f' Formula               = {self.formula}\n'
        if hasattr(self,"adjmat"):     to_print += f' Has Adjacency Matrix  = YES\n'
        else:                          to_print += f' Has Adjacency Matrix  = NO \n'
        if hasattr(self,"parents"):    to_print += f' Number of Parents     = {len(self.parents)}\n'
        if hasattr(self,"occurrence"): to_print += f' Occurrence in Parent  = {self.occurrence}\n'
        if hasattr(self,"totcharge"):  to_print += f' Total Charge          = {self.totcharge}\n'
        if hasattr(self,"spin"):       to_print += f' Spin                  = {self.spin}\n'
        if hasattr(self,"smiles"):     to_print += f' SMILES                = {self.smiles}\n'
        if not indirect: to_print += '------------------------------------------------\n'
        return to_print

###############
### MOLECULE ##
###############
class molecule(specie):
    def __init__(self, labels: list, coord: list, frac_coord: list=None, radii: list=None) -> None:
        self.subtype = "molecule"
        specie.__init__(self, labels, coord, frac_coord, radii)

    ############
    def reset_charge(self):
        specie.reset_charge(self)      ## First uses the generic specie class function for itself and its atoms
        if hasattr(self,"ligands"):    ## Second removes for the child classes 
            for lig in self.ligands:
                lig.reset_charge()
        if hasattr(self,"metals"):    
            for met in self.metals:
                met.reset_charge()

    ############
    def split_complex(self, debug: int=0):
        if not hasattr(self,"atoms"): self.set_atoms()
        if not self.iscomplex:        self.ligands = None; self.metals = None
        else: 
            self.ligands = []
            self.metals  = []
            # Identify Metals and the rest
            metal_idx = list([self.indices[idx] for idx in get_metal_idxs(self.labels, debug=debug)])
            non_transition_metals_idx = list([self.indices[idx] for idx in get_non_transition_metal_idxs(self.labels, debug=debug)])
            if len(non_transition_metals_idx) > 0:
                print(f"MOLECULE.SPLIT COMPLEX: Found non-transition metals in the molecule {self.formula}")
                print(f"MOLECULE.SPLIT COMPLEX: Non-transition metals found: {[self.labels[idx] for idx in non_transition_metals_idx]}")
                metal_idx.extend(non_transition_metals_idx)

            rest_idx  = list(idx for idx in self.indices if idx not in metal_idx) 
            if debug > 0 :  print(f"MOLECULE.SPLIT COMPLEX: labels={self.labels}")
            if debug > 0 :  print(f"MOLECULE.SPLIT COMPLEX: metal_idx={metal_idx}")
            if debug > 0 :  print(f"MOLECULE.SPLIT COMPLEX: rest_idx={rest_idx}")

            # Split the "rest" to obtain the ligands
            rest_labels  = extract_from_list(rest_idx, self.labels, dimension=1)
            rest_coord   = extract_from_list(rest_idx, self.coord, dimension=1)
            if hasattr(self,"frac_coord"): rest_frac    = extract_from_list(rest_idx, self.frac_coord, dimension=1)
            rest_indices = extract_from_list(rest_idx, self.indices, dimension=1)
            rest_radii   = extract_from_list(rest_idx, self.radii, dimension=1)
            rest_atoms   = extract_from_list(rest_idx, self.atoms, dimension=1)
            if debug >= 2: 
                print(f"SPLIT COMPLEX: rest labels: {rest_labels}")
                print(f"SPLIT COMPLEX: rest indices: {rest_indices}")
                print(f"SPLIT COMPLEX: rest radii: {rest_radii}")

            if debug > 0: print(f"SPLIT COMPLEX: splitting species with {len(rest_labels)} atoms in block")
            if hasattr(self,"cov_factor"): blocklist = split_species(rest_labels, rest_coord, radii=rest_radii, cov_factor=self.cov_factor, debug=debug)
            else:                          blocklist = split_species(rest_labels, rest_coord, radii=rest_radii, cov_factor=self.cov_factor, debug=debug)      
            if debug > 0: print(f"SPLIT COMPLEX: received {len(blocklist)} blocks")
            
            ## Arranges Ligands
            for b in blocklist:
                if debug > 0: print(f"PREPARING BLOCK: {b}")
                lig_indices = extract_from_list(b, rest_indices, dimension=1)
                lig_labels  = extract_from_list(b, rest_labels, dimension=1) 
                lig_coord   = extract_from_list(b, rest_coord, dimension=1) 
                if hasattr(self,"frac_coord"): lig_frac_coord = extract_from_list(b, rest_frac, dimension=1)
                lig_radii   = extract_from_list(b, rest_radii, dimension=1) 
                lig_atoms   = extract_from_list(b, rest_atoms, dimension=1) 
                
                if debug > 0: print(f"CREATING LIGAND: {labels2formula(lig_labels)}")
                # Create Ligand Object
                if hasattr(self,"frac_coord"): newligand   = ligand(lig_labels, lig_coord, lig_frac_coord, radii=lig_radii)
                else:                          newligand   = ligand(lig_labels, lig_coord, radii=lig_radii)
                # For debugging
                newligand.origin = "split_complex"
                # Define the molecule as parent of the ligand. Bottom-Up hierarchy
                newligand.add_parent(self, indices=lig_indices)
                
                if self.check_parent("unitcell"):
                    cell_indices = [a.get_parent_index("unitcell") for a in lig_atoms]
                    newligand.add_parent(self.get_parent("unitcell"), indices=cell_indices)

                if self.check_parent("reference"):
                    ref_indices = [a.get_parent_index("reference") for a in lig_atoms]
                    newligand.add_parent(self.get_parent("reference"), indices=ref_indices)

                # Update the ligand with the covalent and metal factors 
                newligand.set_factors(self.cov_factor, self.metal_factor)
                # Pass the molecule atoms to the ligand
                newligand.set_atoms(atomlist=lig_atoms)
                # Inherit the adjacencies from molecule
                newligand.inherit_adjmatrix("molecule")
                # Add ligand to the list. Top-Down hierarchy
                # newligand.evaluate_as_nitrosyl()
                self.ligands.append(newligand)

            ## Arranges Metals
            for m in metal_idx:
                ## We were creating the metal again, but it is already in the list of molecule.atoms
                # newmetal    = metal(self.labels[m], self.coord[m], self.frac_coord[m], self.radii[m])
                # newmetal.add_parent(self, index=self.indices[m])
                # self.metals.append(newmetal)                            
                self.metals.append(self.atoms[m])                            
        return self.ligands, self.metals
        
    #######################################################
    def get_hapticity(self, debug: int=0):
        if not hasattr(self,"ligands"): self.split_complex(debug=debug)
        self.is_haptic = False 
        self.haptic_type = []
        if self.iscomplex: 
            for lig in self.ligands:
                if not hasattr(lig,"is_haptic"):  lig.get_hapticity(debug=debug)
                if lig.is_haptic: self.is_haptic = True
                for entry in lig.haptic_type:
                    if entry not in self.haptic_type: self.haptic_type.append(entry)
        return self.haptic_type
   
    def save(self, path):
        import pickle
        print(f"SAVING Molecule object to {path}")
        with open(path, "wb") as fil:
            pickle.dump(self,fil)

    def __repr__(self):
        to_print = ""
        to_print += f'------------- SCOPE MOLECULE Object --------------\n'
        to_print += specie.__repr__(self, indirect=True)
        if hasattr(self,"ligands"):  
            if self.ligands is not None: to_print += f' # Ligands             = {len(self.ligands)}\n'
        if hasattr(self,"metals"):   
            if self.metals is not None:  to_print += f' # Metals              = {len(self.metals)}\n'
        to_print += '------------------------------------------------\n'
        return to_print

############

###############
### LIGAND ####
###############
class ligand(specie):
    def __init__(self, labels: list, coord: list, frac_coord: list=None, radii: list=None) -> None:
        self.subtype  = "ligand"
        if frac_coord is not None:        
            assert len(frac_coord) == len(coord)
            self.frac_coord = frac_coord
        specie.__init__(self, labels, coord, frac_coord, radii)

    #######################################################
    def set_hapticity(self, hapttype):
        self.hapticity = True 
        self.hapttype  = hapttype 

    #######################################################
    def get_connected_metals(self, debug: int=0):
        # metal.groups will be used for the calculation of the relative metal radius 
        # and define the coordination geometry of the metal /hapicitiy/ hapttype    
        self.metals = []
        mol = self.get_parent("molecule")
        for met in mol.metals:
            tmplabels = self.labels.copy()
            tmpcoord  = self.coord.copy()
            tmplabels.append(met.label)
            tmpcoord.append(met.coord)
            isgood, tmpadjmat, tmpadjnum = get_adjmatrix(tmplabels, tmpcoord, metal_only=True)
            if isgood and any(tmpadjnum) > 0: self.metals.append(met)
        return self.metals
    
    #######################################################
    def get_connected_idx(self, debug: int=0):
        ## Remember madjmat should not be computed at the ligand level. Since the metal is not there.
        ## Now we operate at the molecular level. We get the parent molecule, and the indices of the ligand atoms in the molecule
        self.connected_idx = [] 
        if not hasattr(self,"madjnum"): self.inherit_adjmatrix("molecule")
        for idx, con in enumerate(self.madjnum):
            if con > 0: self.connected_idx.append(idx)
        return self.connected_idx 

    #######################################################
    def get_connected_atoms(self, debug: int=0):
        if not hasattr(self,"atoms"):         self.set_atoms()
        if not hasattr(self,"connected_idx"): self.get_connected_idx()
        self.connected_atoms = []
        for idx, at in enumerate(self.atoms):
            if idx in self.connected_idx and at.madjnum > 0: 
                self.connected_atoms.append(at) 
            elif idx in self.connected_idx and at.madjnum == 0:
                print("WARNING: Atom appears in connected_idx, but has madjnum=0")
        return self.connected_atoms

    #######################################################
    def split_ligand(self, debug: int=0):
        # Split the "ligand to obtain the groups
        self.groups = []
        # Identify Connected and Unconnected atoms (to the metal)
        if not hasattr(self,"connected_idx"): self.get_connected_idx()

        ## Creates the list of variables
        connected_idx     = self.connected_idx
        if debug > 0: print(f"\nLIGAND.SPLIT_LIGAND: splitting {self.formula} into groups")
        if debug >= 2:
            print(f"\tLIGAND.SPLIT_LIGAND: {self.indices=}") 
            print(f"\tLIGAND.SPLIT_LIGAND: {connected_idx=}")
        conn_labels     = extract_from_list(connected_idx, self.labels, dimension=1)
        conn_coord      = extract_from_list(connected_idx, self.coord, dimension=1)
        if hasattr(self,"frac_coord"): conn_frac_coord = extract_from_list(connected_idx, self.frac_coord, dimension=1)
        conn_radii      = extract_from_list(connected_idx, self.radii, dimension=1)
        conn_atoms      = extract_from_list(connected_idx, self.atoms, dimension=1)
        if debug >= 2: print(f"\tLIGAND.SPLIT_LIGAND: {conn_labels=}")

        if hasattr(self,"cov_factor"): blocklist = split_species(conn_labels, conn_coord, radii=conn_radii, cov_factor=self.cov_factor, debug=debug)
        else:                          blocklist = split_species(conn_labels, conn_coord, radii=conn_radii, debug=debug)      
        if debug >= 2: print(f"\tLIGAND.SPLIT_LIGAND: {blocklist=}")
        ## Arranges Groups 
        for b in blocklist:
            if debug >= 2 : print(f"\tLIGAND.SPLIT_LIGAND: block={b}")
            gr_indices = extract_from_list(b, connected_idx, dimension=1, debug=debug)
            if debug > 1: print(f"\tLIGAND.SPLIT_LIGAND: {gr_indices=}")
            gr_labels       = extract_from_list(b, conn_labels, dimension=1, debug=debug)
            gr_coord        = extract_from_list(b, conn_coord, dimension=1)
            if hasattr(self,"frac_coord"): gr_frac_coord   = extract_from_list(b, conn_frac_coord, dimension=1)
            gr_radii        = extract_from_list(b, conn_radii, dimension=1)
            gr_atoms        = extract_from_list(b, conn_atoms, dimension=1)
            # Create Group Object
            if hasattr(self,"frac_coord"): newgroup = group(gr_labels, gr_coord, gr_frac_coord, radii=gr_radii)
            else:                          newgroup = group(gr_labels, gr_coord, radii=gr_radii)
            # For debugging
            newgroup.origin = "split_ligand"
            # Define the ligand as parent of the group. Bottom-Up hierarchy
            newgroup.add_parent(self, indices=gr_indices)
            # Pass the ligand atoms to the groud
            newgroup.set_atoms(atomlist=gr_atoms)
            # Inherit the adjacencies from molecule
            newgroup.inherit_adjmatrix("ligand")
            # Associate the Groups with the Metals
            newgroup.get_connected_metals(debug=debug)
            newgroup.get_closest_metal(debug=debug)
            newgroup.get_hapticity(debug=debug)
            self.groups.append(newgroup)
        if debug > 0 : print(f"\tLIGAND.SPLIT_LIGAND: found groups {[ group.formula for group in self.groups]}")
        if debug > 2 : print(f"{self.groups}")
        return self.groups

    #######################################################
    def get_hapticity(self, debug: int=0):
        if not hasattr(self,"groups"): self.split_ligand(debug=debug)
        self.is_haptic = False 
        self.haptic_type = []
        for gr in self.groups:
            if not hasattr(gr,"is_haptic"): gr.get_hapticity(debug=debug)
            if gr.is_haptic: self.is_haptic = True; self.haptic_type = gr.haptic_type
            for entry in gr.haptic_type:
                if entry not in self.haptic_type: self.haptic_type.append(entry)
        return self.haptic_type

    #######################################################
    def get_denticity(self, debug: int=0):
        if not hasattr(self,"groups"):      self.split_ligand(debug=debug)
        if debug > 1: print(f"LIGAND.Get_denticity: checking connectivity of ligand {self.formula}")
        if debug > 1: print(f"LIGAND.Get_denticity: initial connectivity is {len(self.connected_idx)}")
        self.denticity = 0
        for g in self.groups:
            #if debug > 0: print(f"LIGAND.Get_denticity: checking denticity of group \n{g}\n{g.madjnum=}\n{g.madjmat=}")
            self.denticity += g.get_denticity(debug=debug)      ## A check is also performed at the group level
        if debug > 0: print(f"LIGAND.Get_denticity: final connectivity of ligand {self.formula} is {self.denticity}")
        return self.denticity 

    #######################################################
    def __repr__(self):
        to_print = ""
        to_print += f'------------- SCOPE LIGAND Object --------------\n'
        to_print += specie.__repr__(self, indirect=True)
        if hasattr(self,"rdkit_obj"):  to_print += f' # HAS RDKIT OBJECT    = YES\n'
        else:                          to_print += f' # HAS RDKIT OBJECT    = NO\n'
        to_print += '------------------------------------------------\n'
        return to_print

###############
#### GROUP ####
###############
class group(specie):
    def __init__(self, labels: list, coord: list, frac_coord: list=None, radii: list=None) -> None:
        self.subtype = "group"
        if frac_coord is not None:        
            assert len(frac_coord) == len(coord)
            self.frac_coord = frac_coord
        specie.__init__(self, labels, coord, frac_coord, radii)

    #######################################################
    def remove_atom(self, index: int, debug: int=0):
        if debug > 0: print(f"GROUP.REMOVE_ATOM: deleting atom {index=} from group with {self.natoms} atoms")
        if index > self.natoms: return None
        if not hasattr(self,"atoms"): self.set_atoms()
        self.atoms.pop(index)
        self.labels.pop(index)
        self.coord.pop(index)
        self.radii.pop(index)
        self.formula   = labels2formula(self.labels)
        self.eleccount = labels2electrons(self.labels)   ### Assuming neutral specie (so basically this is the sum of atomic numbers)
        self.natoms    = len(self.labels)
        self.iscomplex = any((elemdatabase.elementblock[l] == "d") or (elemdatabase.elementblock[l] == "f") for l in self.labels)
        if debug > 0: print("GROUP.REMOVE_ATOM. Group after removing atom:")
        if debug > 0: print(self)
        if self.natoms > 0:
            if hasattr(self,"closest_metal"): self.get_closest_metal()
            if hasattr(self,"is_haptic"):     self.get_hapticity()
            if hasattr(self,"centroid"):      self.get_centroid()
            if hasattr(self,"frac_coord"):    self.frac_coord.pop(index)
            if hasattr(self,"adjmat"):        self.get_adjmatrix()
            if hasattr(self,"madjmat"):       self.get_metal_adjmatrix()

    #######################################################
    def get_closest_metal(self, debug: int=0):
        apos = compute_centroid(np.array(self.coord))
        dist = []
        mol  = self.get_parent("molecule")
        for met in mol.metals:
            bpos = np.array(met.coord)
            dist.append(np.linalg.norm(apos - bpos))
        # finds the closest Metal Atom (tgt)
        self.closest_metal = mol.metals[np.argmin(dist)]
        return self.closest_metal

    #######################################################
    def get_connected_metals(self, debug: int=0):
        # metal.groups will be used for the calculation of the relative metal radius 
        # and define the coordination geometry of the metal /hapicitiy/ hapttype    
        self.metals = []
        lig = self.get_parent("ligand")
        if not hasattr(lig,"metals"): lig.get_connected_metals()
        for met in lig.metals:
            tmplabels = self.labels.copy()
            tmpcoord  = self.coord.copy()
            tmplabels.append(met.label)
            tmpcoord.append(met.coord)
            isgood, tmpadjmat, tmpadjnum = get_adjmatrix(tmplabels, tmpcoord, metal_only=True)
            if isgood and any(tmpadjnum) > 0: self.metals.append(met)
        return self.metals

    #######################################################
    def get_hapticity(self, debug: int=0):
        if not hasattr(self,"atoms"): self.set_atoms()
        self.is_haptic   = False ## old self.hapticity
        self.haptic_type = []    ## old self.hapttype
    
        numC  = self.labels.count("C")  # Carbon is the most common connected atom in ligands with hapticity
        numAs = self.labels.count("As") # I've seen one case of a Cp but with As instead of C (VENNEH, Fe dataset)
        numP  = self.labels.count("P")  
        numO  = self.labels.count("O")  # For h4-Enone
        numN  = self.labels.count("N")
    
        ## Carbon-based Haptic Ligands
        if   numC == 2:                   self.haptic_type = ["h2-Benzene", "h2-Butadiene", "h2-ethylene"]; self.is_haptic = True
        elif numC == 3 and numO == 0:     self.haptic_type = ["h3-Allyl", "h3-Cp"];                         self.is_haptic = True
        elif numC == 3 and numO == 1:     self.haptic_type = ["h4-Enone"];                                  self.is_haptic = True
        elif numC == 4:                   self.haptic_type = ["h4-Butadiene", "h4-Benzene"];                self.is_haptic = True
        elif numC == 5:                   self.haptic_type = ["h5-Cp"];                                     self.is_haptic = True
        elif numC == 6:                   self.haptic_type = ["h6-Benzene"];                                self.is_haptic = True
        elif numC == 7:                   self.haptic_type = ["h7-Cycloheptatrienyl"];                      self.is_haptic = True
        elif numC == 8:                   self.haptic_type = ["h8-Cyclooctatetraenyl"];                     self.is_haptic = True
        # Other less common types of haptic ligands
        elif numC == 0 and numAs == 5:    self.haptic_type = ["h5-AsCp"];                                   self.is_haptic = True
        elif numC == 0 and numP == 5:     self.haptic_type = ["h5-Pentaphosphole"];                         self.is_haptic = True
        elif numC == 1 and numP == 1:     self.haptic_type = ["h2-P=C"];                                    self.is_haptic = True
        return self.haptic_type 

    #######################################################
    def get_denticity(self, debug: int=0):
        self.denticity = 0
        for a in self.atoms: 
            if debug > 0: print(f"GROUP.GET_DENTICITY. Evaluating Atom with {a.madjnum=} and so far {self.denticity=}")
            self.denticity += a.madjnum      
        return self.denticity

    #######################################################
    def __repr__(self):
        to_print = ""
        to_print += f'------------- SCOPE GROUP Object --------------\n'
        to_print += specie.__repr__(self, indirect=True)
        to_print += '------------------------------------------------\n'
        return to_print

###############
### BOND ######
###############
class bond(object):
    def __init__(self, atom1: object, atom2: object, bond_order: int=1):
        self.type       = "bond"
        self.version    = "2.0"
        self.atom1      = atom1
        self.atom2      = atom2
        self.order      = bond_order
        self.distance   = np.linalg.norm(np.array(atom1.coord) - np.array(atom2.coord))

    def __repr__(self):
        to_print += f'------------- SCOPE BOND Object --------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        to_print += f' Atom 1                = {self.atom1.parent_index}\n'
        to_print += f' Atom 2                = {self.atom2.parent_index}\n'
        to_print += f' Bond Order            = {self.order}\n'
        to_print += f' Distance              = {round(self.distance,3)}\n'
        to_print += '-------------------------------------------------\n'
        return to_print

###############
### ATOM ######
###############
class atom(object):
    def __init__(self, label: str, coord: list, frac_coord: list=None, radii: float=None) -> None:
        self.type            = "atom"
        self.subtype         = "atom"
        self.version         = "2.0"
        self.label           = label
        self.coord           = coord
        self.atnum           = elemdatabase.elementnr[label]
        self.block           = elemdatabase.elementblock[label]
        self.parents         = []
        self.parents_index   = []
        self.formula         = label

        if frac_coord is not None:        self.frac_coord = frac_coord
        if radii is None:                 self.radii = get_radii(label)
        else:                             self.radii = radii

    #######################################################
    def check_connectivity(self, other: object, debug: int=0):
        ## Checks whether two atoms are connected (through the adjacency)
        if not isinstance(other, type(self)): return False
        labels = list([self.label,other.label])
        coords = list([self.coord,other.coord])
        isgood, adjmat, adjnum = get_adjmatrix(labels, coords)
        if isgood and adjnum[0] > 0: return True
        else:                        return False

    ############
    def add_parent(self, parent: object, index: int, overwrite: bool=True, debug: int=0):
        ## associates a parent specie to self. The atom indices of self in parent are given in "indices"
        ## if parent of the same subtype already in self.parent then it is overwritten
        ## this is to avoid having a substructure (e.g. a ligand) in more than one superstructure (e.g. a molecule) 
        append = True
        for idx, p in enumerate(self.parents):
            if p.subtype == parent.subtype:
                if overwrite: 
                    self.parents[idx]         = parent
                    self.parents_index[idx]   = index
                append = False
        if append: 
            self.parents.append(parent)
            self.parents_index.append(index)
            if debug > 0: print(f"ATOM.ADD_PARENT: added parent with subtype {parent.subtype}. Atom index is {index=}")

    ############
    def check_parent(self, subtype: str):
        ## checks if parent of a given subtype exists
        for p in self.parents:
            if p.subtype == subtype: return True
        return False

    ############
    def get_parent(self, subtype: str):
        ## retrieves parent of a given subtype 
        for p in self.parents:
            if p.subtype == subtype: return p
        return None
    
    ############
    def get_parent_index(self, subtype: str):
        ## retrieves parent of a given subtype 
        for idx, p in enumerate(self.parents):
            if p.subtype == subtype: return self.parents_index[idx]
        return None

    #######################################################
    def inherit_connectivity(self, parent_subtype: str, debug: int=0):
        exists  = self.check_parent(parent_subtype)
        if not exists:
            print(f"ATOM.INHERIT. {parent_subtype=} does not exist")
            return None
        parent  = self.get_parent(parent_subtype)
        index   = self.get_parent_index(parent_subtype)
        index   = list([index])
        assert hasattr(parent, "madjnum")
        assert hasattr(parent, "adjnum")
        if debug > 0: print(f"ATOM.INHERIT. Identified {parent_subtype=} and {index=}") 
        if not hasattr(parent,"madjnum"):
            print(f"ATOM.INHERIT. {parent_subtype=} does not have madjnum")
            return None
        self.madjnum = np.stack(extract_from_list(index, parent.madjnum, dimension=1), axis=0)[0]
        self.adjnum  = np.stack(extract_from_list(index, parent.adjnum, dimension=1), axis=0)[0]

    #######################################################
    def check_connectivity(self, other: object, debug: int=0):
        ## Checks whether two atoms are connected (through the adjacency)
        if not isinstance(other, type(self)): return False
        labels = list([self.label,other.label]) 
        coords = list([self.coord,other.coord]) 
        isgood, adjmat, adjnum = get_adjmatrix(labels, coords)
        if isgood and adjnum[0] > 0: return True
        else:                        return False

    #######################################################
    def add_bond(self, newbond: object, debug: int=0):
        if not hasattr(self,"bonds"): self.bonds = []
        at1 = newbond.atom1
        at2 = newbond.atom2
        found = False
        for b in self.bonds:
            if (b.atom1 == at1 and b.atom2 == at2) or (b.atom1 == at2 and b.atom2 == at1): 
                if debug > 0: print(f"ATOM.ADD_BOND found the same bond with atoms:") 
                if debug > 0: print(f"atom1: {b.atom1}") 
                if debug > 0: print(f"atom2: {b.atom2}") 
                found = True    ### It means that the same bond has already been defined
        if not found: self.bonds.append(newbond)

    #######################################################
    # Replaced by set_factors
    #def set_adjacency_parameters(self, cov_factor: float, metal_factor: float) -> None:
    #    self.cov_factor   = cov_factor
    #    self.metal_factor = metal_factor

    #######################################################
    def reset_charge(self) -> None:
        if hasattr(self,"charge"):    delattr(self,"charge")
        if hasattr(self,"poscharges"): delattr(self,"charge")

    #######################################################
    def set_charge(self, charge: int) -> None:
        self.charge = int(charge)

    #######################################################
    def set_adjacencies(self, adjmat, madjmat, connectivity: int, metal_connectivity: int=0):
        self.adjnum  = int(connectivity)
        self.madjnum = int(metal_connectivity)
        self.adjacency       = []
        self.metal_adjacency = []
        for idx, c in enumerate(adjmat):   ## The atom only receives one row of adjmat, so this is not a matrix anymore. Keep in mind that the idx are the indices of parent
            if c > 0: self.adjacency.append(idx)
        for idx, c in enumerate(madjmat):  ## The atom only receives one row of madjmat, so this is not a matrix anymore
            if c > 0: self.metal_adjacency.append(idx)

    #######################################################
    def get_connected_metals(self, metalist: list, debug: int=0):
        self.metals = []
        for met in metalist:
            tmplabels = self.label.copy()
            tmpcoord  = self.coord.copy()
            tmplabels.append(met.label)
            tmpcoord.append(met.coord)
            isgood, tmpadjmat, tmpadjnum = get_adjmatrix(tmplabels, tmpcoord, metal_only=True)
            if isgood and any(tmpadjnum) > 0: self.metals.append(met)
        return self.metals

    #######################################################
    def get_closest_metal(self, debug: int=0):
        ## Here, the list of metal atoms must be provided
        apos = self.coord
        dist = []
        mol = self.get_parent("molecule")
        for met in mol.metals:
            bpos = np.array(met.coord)
            dist.append(np.linalg.norm(apos - bpos))
        self.closest_metal = mol.metals[np.argmin(dist)]
        return self.closest_metal

    #######################################################
    def set_factors(self, cov_factor: float=1.3, metal_factor: float=1.0) -> None:
        self.cov_factor   = cov_factor
        self.metal_factor = metal_factor

    #######################################################
    def __repr__(self, indirect: bool=False):
        to_print = ""
        if not indirect: to_print += f'------------- SCOPE ATOM Object ----------------\n'
        to_print += f' Version                      = {self.version}\n'
        to_print += f' Type                         = {self.type}\n'
        if hasattr(self,'subtype'): to_print += f' Sub-Type                     = {self.subtype}\n'
        to_print += f' Label                        = {self.label}\n'
        to_print += f' Atomic Number                = {self.atnum}\n'
        if hasattr(self,"occurrence"): to_print += f' Occurrence in Parent         = {self.occurrence}\n'
        if hasattr(self,"charge"):     to_print += f' Atom Charge                  = {self.charge}\n'

        # Adjacency and Metal Adjacency
        if hasattr(self,"madjnum"):    to_print += f' Metal Adjacency (madjnum)    = {self.madjnum}\n'
        elif hasattr(self,"madjnum"):  to_print += f' Metal Adjacency (madjnum)    = {self.madjnum}\n'
        else:                          to_print += f' No Metal Adjacency Info\n'
        if hasattr(self,"connec"):     to_print += f' Regular Adjacencies (connec) = {self.adjnum}\n'
        elif hasattr(self,"adjnum"):   to_print += f' Regular Adjacencies (adjnum) = {self.adjnum}\n'
        else:                          to_print += f' No Adjacency Info\n'

        if not indirect: to_print += '-------------------------------------------------\n'
        return to_print

    #######################################################
    def reset_madjnum(self, met, diff: int=-1, debug: int=0):
        if debug > 0: print(f"ATOM.RESET_MCONN: resetting madjnum (and connec) for atom {self.label=}")
        if debug > 0: print(f"ATOM.RESET_MCONN: initial {self.adjnum=} {self.madjnum=}")
        if debug > 0 : print(f"ATOM.RESET_MCONN: initial = {self.adjacency=} {self.metal_adjacency=}")
        self.madjnum += diff
        self.adjnum  += diff

        if debug > 0: print(f"ATOM.RESET_MCONN: initial {met.adjnum=} {met.madjnum=}")
        if debug > 0 : print(f"ATOM.RESET_MCONN: initial = {met.adjacency=} {met.metal_adjacency=}")

        # Correct Metal Data
        met.madjnum += diff                             # Corrects data of metal object
        met.adjnum  += diff                             # Corrects data of metal object


        exists = self.check_parent("ligand")
        if exists:
            lig     = self.get_parent("ligand")
            lig_idx = self.get_parent_index("ligand")

            if debug > 0: print(f"ATOM.RESET_MCONN: resetting madjnum (and connec) for atom {self.label=} in ligadn {lig_idx=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: updating ligand atoms and madjnum")
            if debug > 0: print(f"ATOM.RESET_MCONN: {lig.natoms=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: {lig.labels=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.madjnum=} {len(lig.madjnum)}") 
            # Nothing in madjmat of the ligand object, all zeros
            # if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.madjmat=} {(lig.madjmat).shape}")
            if debug > 0: print(f"ATOM.RESET_MCONN: updating ligand atoms and adjnum")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.adjnum=} {len(lig.adjnum)}") 
            # if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.adjmat=} {(lig.adjmat).shape}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.atoms[lig_idx].adjnum=} {lig.atoms[lig_idx].madjnum=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.madjnum[lig_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.adjnum[lig_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {met.adjnum=} {met.madjnum=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.madjmat[lig_idx]=} {lig.adjmat[lig_idx]=}")
            # Correct Ligand Data
            lig.madjnum[lig_idx] += diff                    # Corrects data in metal_adjacency number of the ligand class
            #lig.madjmat[lig_idx,met_idx] += diff            # Corrects data in metal_adjacency matrix
            #lig.madjmat[met_idx,lig_idx] += diff            # Corrects data in metal_adjacency matrix
            lig.adjnum[lig_idx]  += diff                    # Corrects data in adjacency number of the ligand class

            lig.atoms[lig_idx].set_adjacencies(lig.adjmat[lig_idx], lig.madjmat[lig_idx], lig.adjnum[lig_idx], lig.madjnum[lig_idx])

            # lig.adjmat[lig_idx,met_idx]  += diff            # Corrects data in adjacency matrix
            # lig.adjmat[met_idx,lig_idx]  += diff            # Corrects data in adjacency matrix
            # we should delete the adjacencies, but not a priority 
            if debug > 0: print(f"ATOM.RESET_MCONN: final {lig.madjnum=}")
            # if debug > 0: print(f"ATOM.RESET_MCONN: final {lig.madjmat=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {lig.adjnum=}")  
            # if debug > 0: print(f"ATOM.RESET_MCONN: final {lig.adjmat=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {lig.atoms[lig_idx].adjnum=} {lig.atoms[lig_idx].madjnum=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {lig.atoms[lig_idx].adjacency=} {lig.atoms[lig_idx].metal_adjacency=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {lig.madjnum[lig_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {lig.adjnum[lig_idx]=}")
            lig.get_connected_idx(debug=debug)
            lig.get_connected_atoms(debug=debug)

        exists = self.check_parent("molecule")
        if exists:
            mol     = self.get_parent("molecule")
            mol_idx = self.get_parent_index("molecule")
            met_idx = met.get_parent_index("molecule")
            if debug > 0: print(f"ATOM.RESET_MCONN: resetting madjnum (and connec) for atom {self.label=} in molecule {mol_idx=} with metal {met_idx=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: updating molecule atoms and madjnum")
            if debug > 0: print(f"ATOM.RESET_MCONN: {mol.natoms=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: {mol.labels=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.madjnum=} {len(mol.madjnum)}") 
            #if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.madjmat=} {(mol.madjmat).shape}") # Nothing in madjmat of the ligand object, all zeros
            if debug > 0: print(f"ATOM.RESET_MCONN: updating molecule atoms and adjnum")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.adjnum=} {len(mol.adjnum)}") 
            #if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.adjmat=} {(mol.adjmat).shape}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.atoms[mol_idx].madjnum=} {mol.atoms[mol_idx].adjnum=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {met.madjnum=} {met.adjnum=}")

            if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.madjnum[mol_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.adjnum[mol_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.madjnum[met_idx]=} {mol.madjnum[mol_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.adjnum[met_idx]=} {mol.adjnum[mol_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.madjmat[met_idx,mol_idx]=} {mol.madjmat[mol_idx,met_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.adjmat[met_idx,mol_idx]=} {mol.adjmat[mol_idx,met_idx]=}")

            # Correct Molecule Data
            mol.madjnum[mol_idx] += diff                    # Corrects data in metal_adjacency number of the molecule class
            mol.madjnum[met_idx] += diff                    # Corrects data in metal_adjacency number of the molecule class

            mol.madjmat[mol_idx,met_idx] += diff            # Corrects data in metal_adjacency matrix
            mol.madjmat[met_idx,mol_idx] += diff            # Corrects data in metal_adjacency matrix
            
            mol.adjnum[mol_idx]  += diff                    # Corrects data in adjacency number of the molecule class
            mol.adjnum[met_idx]  += diff                    # Corrects data in adjacency number of the molecule class
        
            mol.adjmat[mol_idx,met_idx]  += diff            # Corrects data in adjacency matrix
            mol.adjmat[met_idx,mol_idx]  += diff            # Corrects data in adjacency matrix

            self.set_adjacencies(mol.adjmat[mol_idx], mol.madjmat[mol_idx], mol.adjnum[mol_idx], mol.madjnum[mol_idx])

            met.set_adjacencies(mol.adjmat[met_idx], mol.madjmat[met_idx], mol.adjnum[met_idx], mol.madjnum[met_idx])

            if debug > 0: print(f"ATOM.RESET_MCONN: final {mol.atoms[mol_idx].adjnum=} {mol.atoms[mol_idx].madjnum=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {mol.atoms[mol_idx].adjacency=} {mol.atoms[mol_idx].metal_adjacency=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {met.adjnum=} {met.madjnum=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {met.adjacency=} {met.metal_adjacency=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {mol.madjnum[mol_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {mol.adjnum[mol_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {mol.madjnum[met_idx]=} {mol.madjnum[mol_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {mol.adjnum[met_idx]=} {mol.adjnum[mol_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {mol.madjmat[met_idx,mol_idx]=} {mol.madjmat[mol_idx,met_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {mol.adjmat[met_idx,mol_idx]=} {mol.adjmat[mol_idx,met_idx]=}")

###############
#### METAL ####
###############
class metal(atom):
    def __init__(self, label: str, coord: list, frac_coord: list=None, radii: float=None) -> None:
        atom.__init__(self, label, coord, frac_coord=frac_coord, radii=radii)
        self.subtype = "metal"

    #######################################################
    def get_valence_elec (self, m_ox: int):
        """ Count valence electrons for a given transition metal and metal oxidation state """
        v_elec = elemdatabase.valenceelectrons[self.label] - m_ox
        if v_elec >= 0 :  self.valence_elec = v_elec
        else :            self.valence_elec = elemdatabase.elementgroup[self.label] - m_ox
        return self.valence_elec

    #######################################################
    def get_coord_sphere(self, debug: int=0):
        if not self.check_parent("molecule"): 
            print(f"METAL.Get_coord_sphere. Metal does not have parent molecule")
            return None
        mol = self.get_parent("molecule")
        pidx = self.get_parent_index("molecule")
        if not hasattr(mol,"adjmat"): mol.get_adjmatrix()
        adjmat = mol.adjmat.copy()
        
        ## Cordination sphere defined as a collection of atoms
        self.coord_sphere = []
        for idx, at in enumerate(adjmat[pidx]):
            if at >= 1: self.coord_sphere.append(mol.atoms[idx])
        return self.coord_sphere

    #######################################################
    def get_coord_sphere_formula(self, debug: int=0):
        if not hasattr(self,"coord_sphere"): self.get_coord_sphere(debug=debug)
        self.coord_sphere_formula = labels2formula(list([at.label for at in self.coord_sphere])) 
        if debug > 0: print(f"METAL.Get_coord_sphere_formula: {self.get_parent_index('molecule')} {self.label} {self.coord_sphere_formula}")
        return self.coord_sphere_formula 

    #######################################################
    def get_connected_groups(self, debug: int=2):
        from Scope.Adapted_from_cell2mol import split_group
        # metal.groups will be used for the calculation of the relative metal radius 
        # and define the coordination geometry of the metal /hapicitiy/ hapttype    
        if not self.check_parent("molecule"): return None
        mol = self.get_parent("molecule")
        self.groups = []
        for lig in mol.ligands:
            for group in lig.groups:
                if debug > 1: print(group.formula)
                ligand_indices = [ a.get_parent_index("ligand") for a in group.atoms ]
                tmplabels = []
                tmpcoord  = []
                tmplabels.append(self.label)
                tmpcoord.append(self.coord)
                tmplabels.extend(group.labels)
                tmpcoord.extend(group.coord)
                if debug > 1: print(tmplabels, tmpcoord)
                isgood, tmpadjmat, tmpadjnum = get_adjmatrix(tmplabels, tmpcoord, metal_only=True)
                # if isgood and any(tmpadjnum) > 0: self.groups.append(group)
                if isgood:
                    if debug > 1: print(group.formula, tmpadjmat, tmpadjnum)
                    if all(tmpadjnum[1:]): 
                        self.groups.append(group)
                    elif any(tmpadjnum[1:]): 
                        
                        if debug > 1: print(f"Metal {self.label} is connected to {group.formula} but not all atoms are connected")
                        conn_idx = [ idx for idx, num in enumerate(tmpadjnum[1:]) if num == 1 ]
                        conn_ligand_indices = [ ligand_indices[idx] for idx, num in enumerate(tmpadjnum[1:]) if num == 1 ]
                        if debug > 1: print(f"get_connected_groups {tmpadjnum[1:]=} {conn_idx=} {conn_ligand_indices=} {ligand_indices=}")
                        splitted_groups = split_group(group, conn_idx, conn_ligand_indices, debug=debug)
                        for g in splitted_groups:
                            self.groups.append(g)
                            if debug > 1: print(f"Metal {self.label} is connected to {g.formula}")
                    else:
                        if debug > 1: print(f"Metal {self.label} is not connected to {group.formula}")
        return self.groups

    #######################################################
    def get_relative_metal_radius(self, debug: int=0):
        if not hasattr(self,"groups"): self.get_connected_groups(debug=debug)
        diff_list = []
        for group in self.groups:
            if group.is_haptic == False :
                for atom in group.atoms:
                    diff = round(get_dist(self.coord, atom.coord) - elemdatabase.CovalentRadius3[atom.label], 3)
                    diff_list.append(diff)
            else :
                haptic_center_label = "C"
                haptic_center_coord = compute_centroid(np.array([atom.coord for atom in group.atoms]))
                diff = round(get_dist(self.coord, haptic_center_coord) - elemdatabase.CovalentRadius3[haptic_center_label], 3)
                diff_list.append(diff)     
        average = round(np.average(diff_list), 3)   
        if debug > 1: 
            print(f"METAL.Get_relative_metal_radius: {diff_list=}")
            print(f"METAL.Get_relative_metal_radius: {average=}") 
        self.rel_metal_radius = round(average/elemdatabase.CovalentRadius3[self.label], 3)
        return self.rel_metal_radius

#######################################################
    def get_connected_metals(self, debug: int=2):
        self.metals = []
        mol = self.get_parent("molecule")
        for met in mol.metals:
            if met == self : continue
            tmplabels = []
            tmpcoord  = []
            tmplabels.append(self.label)
            tmpcoord.append(self.coord)
            tmplabels.append(met.label)
            tmpcoord.append(met.coord)
            if debug > 1: print(tmplabels, tmpcoord)
            isgood, tmpadjmat, tmpadjnum = get_adjmatrix(tmplabels, tmpcoord, metal_only=True)
            if isgood:
                if debug > 1: print(met.label, tmpadjmat, tmpadjnum)
                if all(tmpadjnum[1:]): 
                    self.metals.append(met)
                else:
                    if debug > 1: print(f"Metal {self.label} is not connected to {met.label}")
        if debug >= 2 : print(f"METAL.Get_connected_metals: {self.label} connected to {len(self.metals)} metals {[m.label for m in self.metals]}")
        return self.metals
    
    ############
    def reset_charge(self):
        atom.reset_charge(self)     ## First uses the generic atom class function for itself
        if hasattr(self,"poscharges"):   delattr(self,"poscharge")

    def __repr__(self):
        to_print = ""
        to_print += f'------------- SCOPE METAL Object --------------\n'
        to_print += atom.__repr__(self, indirect=True)
        to_print += '-------------------------------------------------\n'
        return to_print


##############
#### CELL ####
##############
class cell(object):
    def __init__(self, name: str, labels: list, pos: list, frac_coord: list, cell_vector: list, cell_param: list) -> None:
        self.version    = "2.0"
        self.type       = "cell"
        self.subtype    = "cell"
        self.name       = name
        self.refcode    = name
        self.labels     = labels 
        self.coord      = pos
        self.frac_coord = frac_coord
        self.cell_vector = cell_vector
        self.cell_param  = cell_param
        self.natoms     = len(labels)

    #######################################################    
    def set_subtype(self, subtype):
        self.subtype    = subtype

    #######################################################
    ## This function mimics the specie-class function with the same name
    def check_parent(self, subtype):
        if subtype == "cell": return True
        else:                 return False

    #######################################################
    ## This function mimics the specie-class function with the same name
    def get_parent(self, subtype):
        if subtype == "cell": return self
        else:                 return None

    #######################################################
    ## This function mimics the specie-class function with the same name
    def get_parent_indices(self, subtype: str):
        if subtype == "cell": return list(range(0,self.natoms))
        else:                 return None

    #######################################################
    def get_unique_species(self, debug: int=0): 
        #if not hasattr(self,"is_fragmented"): self.reconstruct(debug=debug)  
        #if self.is_fragmented: return None # Stopping. self.is_fragmented must be false to determine the charges of the cell
        
        if debug >= 0: print(f"Getting unique species in {self.subtype}")
        self.unique_species = []
        self.unique_indices = []
        self.species_list = []
        typelist_mols = [] # temporary variable  
        typelist_ligs = [] # temporary variable
        typelist_mets = [] # temporary variable

        specs_found = -1
        if self.type == "cell" and not hasattr(self,"subtype"): self.set_subtype("cell")
        if self.subtype == "reference": moleclist = self.refmoleclist
        else:                           moleclist = self.moleclist
        for idx, mol in enumerate(moleclist):
            if debug >= 2: print(f"Molecule {idx} formula={mol.formula}")
            if not mol.iscomplex:
                found = False
                for ldx, typ in enumerate(typelist_mols):   # Molecules
                    issame = compare_species(mol, typ[0], debug=0)
                    if issame :
                        found = True ; kdx = typ[1]
                        if debug >= 2: print(f"Molecule {idx} is the same with {ldx} in typelist")
                if not found:
                    specs_found += 1 ; kdx = specs_found
                    typelist_mols.append(list([mol, kdx]))
                    self.unique_species.append(mol)
                    if debug >= 2: print(f"New molecule found with: formula={mol.formula} and added in position {kdx}")
                self.unique_indices.append(kdx)
                mol.unique_index = kdx
                self.species_list.append(mol)
            else:
                if not hasattr(mol,"ligands"): mol.split_complex(debug=debug)
                for jdx, lig in enumerate(mol.ligands):     # ligands
                    found = False
                    for ldx, typ in enumerate(typelist_ligs):
                        if not hasattr(lig, "is_nitrosyl"): lig.evaluate_as_nitrosyl()
                        if not hasattr(typ[0], "is_nitrosyl"): typ[0].evaluate_as_nitrosyl()
                        if lig.is_nitrosyl and typ[0].is_nitrosyl: 
                            if lig.NO_type == typ[0].NO_type: 
                                issame = True
                            else:
                                issame = False
                        else:                            
                            issame = compare_species(lig, typ[0], debug=0)
                        if issame :
                            found = True ; kdx = typ[1]
                            if debug >= 2: print(f"ligand {jdx} is the same with {ldx} in typelist")
                    if not found:
                        specs_found += 1 ; kdx = specs_found
                        typelist_ligs.append(list([lig, kdx]))
                        self.unique_species.append(lig)
                        if debug >= 2: print(f"New ligand found with: formula {lig.formula} added in position {kdx}")
                    self.unique_indices.append(kdx)
                    lig.unique_index = kdx
                    self.species_list.append(lig)
                for jdx, met in enumerate(mol.metals):      #  metals
                    found = False
                    for ldx, typ in enumerate(typelist_mets):
                        issame = compare_metals(met, typ[0], debug=0)
                        if issame :
                            found = True ; kdx = typ[1]
                            if debug >= 2: print(f"Metal {jdx} is the same with {ldx} in typelist")
                    if not found:
                        specs_found += 1 ; kdx = specs_found
                        typelist_mets.append(list([met, kdx]))
                        self.unique_species.append(met)
                        if debug >= 2: print(f"New Metal Center found with: labels {met.label} and added in position {kdx}")
                    self.unique_indices.append(kdx)
                    met.unique_index = kdx
                    self.species_list.append(met)
        return self.unique_species

    ######################################################
    def get_reference_molecules(self, ref_labels: list, ref_fracs: list, cov_factor: float=1.3, metal_factor: float=1.0, debug: int=0):
        if debug >= 0:
            print("#########################################")
            print("  GETREFS: Generate reference molecules  ")
            print("#########################################")

        # Convert fractional coordinates to cartesian
        ref_pos = frac2cart_fromparam(ref_fracs, self.cell_param)
        
        # Define reference cell
        refcell = cell(self.name, ref_labels, ref_pos, ref_fracs, self.cell_vector, self.cell_param)
        refcell.set_subtype("reference")
        # Get reference molecules
        blocklist = split_species(ref_labels, ref_pos, cov_factor=cov_factor)
        self.refmoleclist = []
        for b in blocklist:
            mol_labels       = extract_from_list(b, ref_labels, dimension=1)
            mol_coord        = extract_from_list(b, ref_pos, dimension=1)
            mol_frac_coord   = extract_from_list(b, ref_fracs, dimension=1)
            newmolec         = molecule(mol_labels, mol_coord, mol_frac_coord)
            newmolec.add_parent(self, indices=b)
            newmolec.add_parent(refcell, indices=b)
            newmolec.set_factors(cov_factor, metal_factor)
            newmolec.set_atoms(create_adjacencies=True, debug=debug)
            for atom, idx in zip(newmolec.atoms, b):
                atom.add_parent(refcell, index=idx)
            # This must be below the frac_coord, so they are carried on to the ligands
            if newmolec.iscomplex: 
                newmolec.split_complex()
            else:
                newmolec.add_parent(newmolec, indices=[*range(0,newmolec.natoms,1)])
            self.refmoleclist.append(newmolec)
        
        if debug >= 0: print(f"GETREFS: found {len(self.refmoleclist)} reference molecules")
        if debug >= 0: print(f"GETREFS:", [ref.formula for ref in self.refmoleclist])
        # Checks for isolated atoms, and retrieves warning if there is any. Except if it is H, halogen (group 17) or alkalyne (group 2)
        isgood = True 
        for ref in self.refmoleclist:
            if ref.natoms == 1:
                label = ref.atoms[0].label
                group = elemdatabase.elementgroup[label]
                if label == "H" or label == "D": 
                    isgood = False
                else: #(group == 1 or group == 2 or group == 17)
                    if debug >= 0: print(f"GETREFS: found ref molecule with only one atom {ref.labels}")

        # If all good, then works with the reference molecules
        if isgood:
            self.has_isolated_H = False
            for ref in self.refmoleclist:
                if debug >= 0: print(f"GETREFS: working with {ref.formula}")
                if ref.iscomplex: 
                    ref.get_hapticity(debug=debug)
                    for lig in ref.ligands:
                        lig.get_denticity(debug=debug)
                    for met in ref.metals:
                        met.get_connected_metals(debug=debug)                         
                        met.get_coord_sphere_formula()
        else:      
            self.has_isolated_H = True
            
        return self.refmoleclist

    ######################################################
    def get_adjmatrix(self):
        isgood, adjmat, adjnum = get_adjmatrix(self.labels, self.coord, self.cov_factor, self.metal_factor, get_radii(self.labels))
        if isgood:
            self.adjmat = adjmat
            self.adjnum = adjnum
        else:
            self.adjmat = None
            self.adjnum = None
        return self.adjmat, self.adjnum

    ######################################################
    def set_moleclist(self, moleclist: list) -> None:
        self.moleclist = moleclist

    ######################################################
    def check_fragmentation(self, reconstruct: bool = False, debug: int=0):

        if not hasattr(self,"moleclist"): self.get_moleclist()
        self.fragmented = False

        # First comparison with current moleclist
        for mol in self.moleclist:
            found = False
            for rmol in self.refmoleclist:
                issame = compare_species(mol, rmol)  
                if issame: found = True
            if not found: self.fragmented = True

        # If there are fragments and user wants reconstruction, it tries to reconstruct and checks the new moleclist
        if self.fragmented and reconstruct:
            new_moleclist = self.reconstruct(debug=debug)
            self.fragmented = False
            for mol in new_moleclist:
                found = False
                for rmol in self.refmoleclist:
                    issame = compare_species(mol, rmol)  
                    if issame: found = True
                if not found: self.fragmented = True
            if not self.fragmented: 
                self.moleclist = new_moleclist
               #self.set_geometry_from_moleclist()

        return self.fragmented

    ######################################################
    def get_moleclist(self, overwrite: bool=False, cov_factor: float=1.3, metal_factor: float=1.0, debug: int=0):

        ## Overwrite and Warning
        if not overwrite and hasattr(self,"moleclist"): 
            if debug > 0: print(f"CELL.MOLECLIST. Moleclist already exists and default is overwrite=False")
            return self.moleclist

        ## Security
        if not hasattr(self,"labels") or not hasattr(self,"coord"): 
            if debug > 0: print(f"CELL.MOLECLIST. Labels or coordinates not found. Returning None")
            return None
        if len(self.labels) == 0 or len(self.coord) == 0:           
            if debug > 0: print(f"CELL.MOLECLIST. Empty labels or coordinates. Returning None")
            return None
        if debug > 0: print(f"CELL.MOLECLIST passed initial checks")

        ## If fragmented, then it reconstructs the unit cell
        if not hasattr(self,"fragmented"): self.check_fragmentation()
        if self.fragmented:
            print("CELL.MOLECLIST. Reconstructing unit cell")
            self.reconstruct(debug=debug)
            print("CELL.MOLECLIST. Arranging new coordinates from moleclist")
            new_coord = self.arrange_cell_coord()
            if new_coord is None:
                print("CELL.MOLECLIST. WARNING received from ARRANGE_CELL_COORD. Stopping")
                return 
        
        if debug > 0: print("CELL.MOLECLIST. Getting Ions")
        ions_idx = []
        for ref in self.refmoleclist:
            if ref.natoms == 1:
                label = ref.atoms[0].label
                ions_idx.extend([idx for idx, l in enumerate(self.labels) if l == label])
                if debug > 1: print(f"CELL.MOLECLIST: {ions_idx=} with {label=}")        
        
        if len(ions_idx) > 0: 
            cell_indices = [*range(0,len(self.labels),1)]
            if debug > 1: print(f"CELL.MOLECLIST: found {len(ions_idx)} ions")
            rest_idx  = list(idx for idx in cell_indices if idx not in ions_idx)
            if debug > 1: print(f"CELL.MOLECLIST: {rest_idx=}")
            rest_labels  = extract_from_list(rest_idx, self.labels, dimension=1)
            rest_coord   = extract_from_list(rest_idx, self.coord, dimension=1)
            rest_indices = extract_from_list(rest_idx, cell_indices, dimension=1)
            if debug > 1: print(f"CELL.MOLECLIST: {rest_labels=}")
            if debug > 1: print(f"CELL.MOLECLIST: {rest_coord=}")
            if debug > 1: print(f"CELL.MOLECLIST: {rest_indices=}")
            blocklist = split_species(rest_labels, rest_coord, indices=rest_indices, cov_factor=cov_factor, debug=debug)
            for idx in ions_idx: blocklist.append([idx])    
        else :
            blocklist = split_species(self.labels, self.coord, cov_factor=cov_factor, debug=debug)
        
        if blocklist is None: 
            return None
        else :
            if debug > 0: print(f"CELL.MOLECLIST: found {len(blocklist)} blocks")
            if debug > 0: print(f"CELL.MOLECLIST: {blocklist=}")
        
        self.moleclist = []
        for b in blocklist:
            if debug > 0: print(f"CELL.MOLECLIST: doing block={b}")
            mol_labels  = extract_from_list(b, self.labels, dimension=1)
            mol_coord   = extract_from_list(b, self.coord, dimension=1)
            mol_frac_coord  = extract_from_list(b, self.frac_coord, dimension=1)
            # Creates Molecule Object
            newmolec    = molecule(mol_labels, mol_coord, mol_frac_coord)
            # For debugging
            newmolec.origin = "cell.get_moleclist"
            # Adds cell as parent of the molecule, with indices b
            newmolec.add_parent(self, indices=b)            
            # Store Adjacency Parameters
            newmolec.set_factors(cov_factor, metal_factor)
            # Creates The atom objects with adjacencies
            newmolec.set_atoms(create_adjacencies=True, debug=debug)
            # The split_complex must be below the frac_coord, so they are carried on to the ligands
            if newmolec.iscomplex: 
                if debug > 0: print(f"CELL.MOLECLIST: splitting complex")
                newmolec.split_complex(debug=debug)
            # Not needed here, as the reconstruction will take care of it
            self.moleclist.append(newmolec)

        return self.moleclist
   
    def arrange_cell_coord(self): 
        ## Updates the cell coordinates preserving the original atom ordering
        ## Do do so, it uses the variable atlist stored in each molecule
        new_coord = np.zeros((self.natoms,3))
        done = np.zeros((self.natoms)) ## Security check so that all coordinates are modified just once
        for mol in self.moleclist:
            for z in zip(mol.indices, mol.coord):
                for i in range(0,3):
                    new_coord[z[0]][i] = z[1][i]
                done[z[0]] += 1
        if np.any(done != 1): 
            print("CELL.ARRANGE_CELL_COORDINATES. WARNING!!! not all coordinates have been modified once")
            print(f"CELL.ARRANGE_CELL_COORDINATES. {done=}")
            return None
        self.coord = np.ndarray.tolist(new_coord)
        return self.coord

    def get_occurrence(self, substructure: object) -> int:
        occurrence = 0
        ## Molecules in Cell
        if hasattr(substructure,"subtype") and hasattr(self,"moleclist"): 
            if substructure.subtype == 'molecule':
                for m in self.moleclist:
                    issame = compare_species(substructure, m)
                    if issame: occurrence += 1
        return occurrence

#######################################################
    def reconstruct(self, cov_factor: float=None, metal_factor: float=None, debug: int=0):
        from Scope.Reconstruct import classify_fragments, fragments_reconstruct

        if not hasattr(self,"fragmented"): self.check_fragmentation()
        if not self.fragmented:
            print("CELL.RECONSTRUCT. Cell is not fragmented")
            return self.moleclist

        if not hasattr(self,"refmoleclist"): print("CELL.RECONSTRUCT. CELL missing list of reference molecules"); return
        if cov_factor is None:   cov_factor   = self.refmoleclist[0].cov_factor
        if metal_factor is None: metal_factor = self.refmoleclist[0].metal_factor

        ## Get the fragments, which is the moleclist of a fragmented cell
        fragments = self.get_moleclist(cov_factor=cov_factor, metal_factor=metal_factor, debug=debug)
        if fragments is None: self.error_get_fragments = True; return  
        else:                 self.error_get_fragments = False

        ## Classifies fragments
        molecules, fragments, hydrogens = classify_fragments(fragments, self.refmoleclist, debug=debug)
        if debug > 0: print(f"CELL.RECONSTRUCT: {len(molecules)} {molecules=}")
        if debug > 0: print(f"CELL.RECONSTRUCT: {len(fragments)} {fragments=}")
        if debug > 0: print(f"CELL.RECONSTRUCT: {len(hydrogens)} {hydrogens=}")

        ## Determines if Reconstruction is necessary
        if len(fragments) > 0 or len(hydrogens) > 0: self.is_fragmented = True
        else:                                        self.is_fragmented = False
        
        self.moleclist = []
        if not self.is_fragmented: 
            for mol in molecules:
                self.moleclist.append(mol)
            return self.moleclist     
        else :
            reconstructed_molecules, Warning = fragments_reconstruct(molecules, fragments, hydrogens, self.refmoleclist, self.cell_vector, cov_factor, metal_factor)
            
            if Warning:
                self.is_fragmented = True
                self.error_reconstruction = True 

            else :
                self.is_fragmented = False
                self.error_reconstruction = False 

            ## For consistency, we create the molecules once again, even if mol is already a molecule-class object.
            ## One must follow the same structure as in self.get_moleclist()
            if debug > 0: print(f"CELL.RECONSTRUCT: Creating and preparing molecules")
            for idx, mol in enumerate(reconstructed_molecules):
                if debug > 0: print(f"CELL.RECONSTRUCT: Doing molecule {idx} with {mol.formula=}")
                newmolec = molecule(mol.labels, mol.coord)
                newmolec.origin = "cell.reconstruct"
                newmolec.set_factors(cov_factor, metal_factor)
                newmolec.set_atoms(create_adjacencies=True, debug=debug)
                newmolec.add_parent(self, mol.cell_indices, debug=debug) 
                if debug > 0: print(f"CELL.RECONSTRUCT: Setting fractional coordinates")
                newmolec.set_fractional_coord(mol.frac_coord)
                if newmolec.iscomplex: newmolec.split_complex()
                self.moleclist.append(newmolec)         
            return self.moleclist

    #######################################################
    def reset_charge_assignment(self, debug: int=0):
        if not hasattr(self,"moleclist"): return None
        for mol in self.moleclist:
            mol.reset_charge()

    #######################################################
    def view(self, size: str='default'):
        import plotly.graph_objects as go
        from Scope.Read_Write import set_scene
        from Scope.Elementdata import ElementData  
        elemdatabase = ElementData()

        ### Adjusts size
        if   size.lower() == 'default': width=600; height=600;   marker_size=8;  text_size=9
        elif size.lower() == 'small':   width=400; height=400;   marker_size=6;  text_size=7 
        elif size.lower() == 'large':   width=800; height=800;   marker_size=10; text_size=12
        elif size.lower() == 'ultra':   width=1000; height=1000; marker_size=11; text_size=13

        if not hasattr(self,"adjmat"): self.get_adjmatrix()
        fig             = go.Figure()

        # Gather Data
        positions       = np.array(self.coord)
        symbols         = self.labels
        adjacencies     = self.adjmat

        # Gets bonds from adjacency matrix
        indices = np.argwhere(self.adjmat >0)
        unique_bonds    = set()

        for i in indices:
            unique_bonds.add(tuple(i))

        # Plot atoms as markers
        fig.add_trace(go.Scatter3d(
            x           = positions[:, 0],
            y           = positions[:, 1],
            z           = positions[:, 2],
            mode        ='markers',
            marker      = dict(
                size        = 10,
                color       = [elemdatabase.cpk_colors[l] for l in symbols],
                line        = dict(color='black', width=1),
            ),
            hoverinfo   = 'text',
            text        = symbols,
            showlegend  = False
        ))

        # Label atom + indices
        fig.add_trace(go.Scatter3d(
            x           = positions[:, 0],
            y           = positions[:, 1],
            z           = positions[:, 2],
            mode        = 'text',
            text        = self.labels,
            #text        = [str(i) for i in range(len(positions))],
            textfont    = dict(color='black', size=12),
            hoverinfo   = 'none',
            showlegend  = False
        ))

        ## Plot bonds as lines and calculate midpoints
        #midpoints   = []
        bond_pairs  = []

        for i, j in unique_bonds:
            # Add bond trace
            fig.add_trace(go.Scatter3d(
                x           = [positions[i, 0], positions[j, 0]],
                y           = [positions[i, 1], positions[j, 1]],
                z           = [positions[i, 2], positions[j, 2]],
                mode        = 'lines',
                line        = dict(color='gray', width=5),
                hoverinfo   = 'none',
                showlegend  = False
            ))

        #    # Calculate midpoints
        #    midpoint = (positions[i] + positions[j]) / 2
        #    midpoints.append(midpoint)
            bond_pairs.append((i, j))

        #midpoints = np.array(midpoints)

        set_scene(fig, np.array(self.coord), width=width, height=height)
        fig.show()

    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> CELL >>>                                    \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        to_print += f' Refcode               = {self.refcode}\n'
        to_print += f' Num Atoms             = {self.natoms}\n'
        to_print += f' Cell Parameters a:c   = {self.cell_param[0:3]}\n'
        to_print += f' Cell Parameters al:ga = {self.cell_param[3:6]}\n'
        if hasattr(self,"moleclist"):  
            to_print += f' # Molecules:          = {len(self.moleclist)}\n'
            to_print += f' With Formulae:                               \n'
            for idx, m in enumerate(self.moleclist):
                to_print += f'    {idx}: {m.formula} \n'
        to_print += '---------------------------------------------------\n'
        if hasattr(self,"refmoleclist"):
            to_print += f' # of Ref Molecules:   = {len(self.refmoleclist)}\n'
            to_print += f' With Formulae:                                  \n'
            for idx, ref in enumerate(self.refmoleclist):
                to_print += f'    {idx}: {ref.formula} \n'
        return to_print


#####################
### SPLIT SPECIES ### 
#####################
#def split_species_from_object(obj: object, debug: int=0):
#
#    if not hasattr(obj,"adjmat"): obj.get_adjmatrix()
#    if obj.adjmat is None: return None
#
#    degree = np.diag(self.adjnum)  # creates a matrix with adjnum as diagonal values. Needed for the laplacian
#    lap = obj.adjmat - degree     # computes laplacian
#
#    # creates block matrix
#    graph = csr_matrix(lap)
#    perm = reverse_cuthill_mckee(graph)
#    gp1 = graph[perm, :]
#    gp2 = gp1[:, perm]
#    dense = gp2.toarray()
#
#    # detects blocks in the block diagonal matrix called "dense"
#    startlist, endlist = get_blocks(dense)
#
#    nblocks = len(startlist)
#    # keeps track of the atom movement within the matrix. Needed later
#    atomlist = np.zeros((len(dense)))
#    for b in range(0, nblocks):
#        for i in range(0, len(dense)):
#            if (i >= startlist[b]) and (i <= endlist[b]):
#                atomlist[i] = b + 1
#    invperm = inv(perm)
#    atomlistperm = [int(atomlist[i]) for i in invperm]
#
#    # assigns atoms to molecules
#    blocklist = []
#    for b in range(0, nblocks):
#        atlist = []    # atom indices in the original ordering
#        for i in range(0, len(atomlistperm)):
#            if atomlistperm[i] == b + 1:
#                atlist.append(indices[i])
#        blocklist.append(atlist)
#    return blocklist

######################
####    IMPORT    ####
######################
def import_cell(old_cell: object, debug: int=0) -> object:
    assert hasattr(old_cell,"labels") 
    assert hasattr(old_cell,"coord") or hasattr(old_cell,"pos")
    assert hasattr(old_cell,"cellvec") or hasattr(old_cell,"cell_vector")
    assert hasattr(old_cell,"refcode")

    labels     = old_cell.labels
    refcode    = old_cell.refcode
    if   hasattr(old_cell,"coord"):       coord      = old_cell.coord
    elif hasattr(old_cell,"pos"):         coord      = old_cell.pos

    if   hasattr(old_cell,"cellvec"):     cell_vector = old_cell.cellvec
    elif hasattr(old_cell,"cell_vector"): cell_vector = old_cell.cell_vector

    if   hasattr(old_cell,"cell_param"):  cell_param  = old_cell.cell_param
    elif hasattr(old_cell,"cellparam"):   cell_param  = old_cell.cellparam
    else:                                 cell_param  = cellvec_2_cellparam(cell_vector)

    if   hasattr(old_cell,"frac_coord"):  frac_coord = old_cell.frac_coord
    else:                                 frac_coord = cart2frac(coord, cell_vector)

    new_cell = cell(refcode, labels, coord, frac_coord, cell_vector, cell_param)
    new_cell.subtype = "cell"
    new_cell.origin  = "import_cell"
    if debug > 0: print(f"IMPORT CELL: importing cell {new_cell}")

    ## Moleclist
    if debug > 0: print(f"IMPORT CELL: creating moleclist")
    moleclist = []
    for mol in old_cell.moleclist: 
        new_mol = import_molecule(mol, parent=new_cell, debug=debug)
        moleclist.append(new_mol)
    new_cell.set_moleclist(moleclist)

    ## Refmoleclist
    if debug > 0: print(f"IMPORT CELL: creating refmoleclist")
    new_cell.refmoleclist = []
    if hasattr(old_cell,"refmoleclist"):
        for rmol in old_cell.refmoleclist:
            new_cell.refmoleclist.append(import_molecule(rmol, parent=new_cell))
    elif hasattr(new_cell,"moleclist"):
        for mol in new_cell.moleclist:
            found = False
            for rmol in new_cell.refmoleclist:
                issame = compare_species(rmol, mol)
                if issame: found = True 
            if not found: new_cell.refmoleclist.append(mol)

    ## Temporary things that I'd like to remove from the import once sorted 
    if hasattr(old_cell,"warning_list"): new_cell.warning_list = old_cell.warning_list 

    return new_cell

################################
def import_gmol(gmol: object, debug: int=0) -> object:
    assert hasattr(gmol,"labels") 
    assert hasattr(gmol,"coord") or hasattr(gmol,"pos")

    labels     = gmol.labels

    if   hasattr(gmol,"coord"):             coord      = gmol.coord
    elif hasattr(gmol,"pos"):               coord      = gmol.pos

    if   hasattr(gmol,"parent_indices"):    indices    = gmol.parent_indices
    elif hasattr(gmol,"atlist"):            indices    = gmol.atlist
    else:                                   indices    = None

    if   hasattr(gmol,"radii"):             radii      = gmol.radii
    else:                                   radii      = None          

    if hasattr(gmol,"factor") and hasattr(gmol,"metal_factor"): 
        cov_factor   = gmol.factor
        metal_factor = gmol.metal_factor
    elif hasattr(gmol,"cov_factor") and hasattr(gmol,"metal_factor"): 
        cov_factor   = gmol.cov_factor
        metal_factor = gmol.metal_factor
    else:
        cov_factor   = 1.3
        metal_factor = 1.0

    new_gmol = specie(labels, coord, radii)
    new_gmol.origin = "import_gmol"
    new_gmol.set_factors(cov_factor, metal_factor)
    new_gmol.get_adjmatrix()           ## Necessary when importing atoms
    new_gmol.get_metal_adjmatrix()     ## Necessary when importing atoms
    if debug > 0: print(f"IMPORT GMOL: importing gmol {new_gmol.formula}")

    ## Smiles
    if   hasattr(gmol,"Smiles"): new_gmol.smiles = gmol.Smiles
    elif hasattr(gmol,"smiles"): new_gmol.smiles = gmol.smiles        
    elif debug > 0: print(f"IMPORT GMOL: SMILES could not be imported")

    ## Atoms
    if not hasattr(gmol,"atoms"):  
        if debug > 0: print(f"IMPORT GMOL: Creating Atoms")
        new_gmol.set_atoms(create_adjacencies=True, debug=debug)
    else:
        if debug > 0: print(f"IMPORT GMOL: importing atoms from old_molecule")
        atoms = []
        for kdx, at in enumerate(gmol.atoms):
            new_atom = import_atom(at, parent=new_gmol, index=kdx, debug=debug)
            atoms.append(new_atom)
        new_molec.set_atoms(atomlist=atoms, debug=debug)
    if debug > 0: print(f"IMPORT GMOL: Example of imported atom")
    if debug > 0: print(new_gmol.atoms[0])

    ## Charges
    if hasattr(gmol,"totcharge") and hasattr(gmol,"atcharge"):
        if debug > 0: print(f"IMPORT gmolEC: imported total charge and atomic charges")
        new_gmol.set_charges(gmol.totcharge, gmol.atcharge)
    elif hasattr(gmol,"totcharge") and not hasattr(gmol,"atcharge"):
        if debug > 0: print(f"IMPORT gmolEC: imported total charge but no atomic charges")
        new_gmol.set_charges(gmol.totcharge)

    ## Fractional coordinates
    if debug > 0: print(f"IMPORT gmolEC: trying to import fractional coordinates")
    if   hasattr(gmol,"frac_coord"):   new_gmol.set_fractional_coord(gmol.frac_coord, debug=debug)
    elif hasattr(gmol,"frac"):         new_gmol.set_fractional_coord(gmol.frac, debug=debug)
    else:                              new_gmol.get_fractional_coord(debug=debug)

    return new_gmol

################################
def import_molecule(mol: object, parent: object=None, debug: int=0) -> object:
    assert hasattr(mol,"labels") 
    assert hasattr(mol,"coord") or hasattr(mol,"pos")

    labels     = mol.labels

    if   hasattr(mol,"coord"):             coord      = mol.coord
    elif hasattr(mol,"pos"):               coord      = mol.pos

    if   hasattr(mol,"parent_indices"):    indices    = mol.parent_indices
    elif hasattr(mol,"atlist"):            indices    = mol.atlist
    else:                                  indices    = None

    if   hasattr(mol,"radii"):             radii      = mol.radii
    else:                                  radii      = None          

    if hasattr(mol,"factor") and hasattr(mol,"metal_factor"): 
        cov_factor   = mol.factor
        metal_factor = mol.metal_factor
    elif hasattr(mol,"cov_factor") and hasattr(mol,"metal_factor"): 
        cov_factor   = mol.cov_factor
        metal_factor = mol.metal_factor
    else:
        cov_factor   = 1.3
        metal_factor = 1.0

    new_molec = molecule(labels, coord, radii)
    new_molec.origin = "import_molecule"
    new_molec.set_factors(cov_factor, metal_factor)
    new_molec.get_adjmatrix()           ## Necessary when importing atoms
    new_molec.get_metal_adjmatrix()     ## Necessary when importing atoms
    if debug > 0: print(f"IMPORT MOLEC: importing molecule {new_molec.formula}")

    ## Parents
    if parent is not None:
        new_molec.add_parent(parent, indices, overwrite=False, debug=debug) 
        if debug > 0: print(f"IMPORT MOLEC: parent {parent.subtype=} added with {indices=}")
    else:
        if debug > 0: print(f"IMPORT MOLEC: parent is None") 

    ## Smiles
    if   hasattr(mol,"Smiles"): new_molec.smiles = mol.Smiles
    elif hasattr(mol,"smiles"): new_molec.smiles = mol.smiles        
    elif debug > 0: print(f"IMPORT MOLEC: SMILES could not be imported")

    ## Atoms
    if not hasattr(mol,"atoms"):  
        if debug > 0: print(f"IMPORT MOLEC: Creating Atoms")
        new_molec.set_atoms(create_adjacencies=True, debug=debug)
    else: 
        if debug > 0: print(f"IMPORT MOLEC: importing atoms from old_molecule")
        atoms = []
        for kdx, at in enumerate(mol.atoms): 
            new_atom = import_atom(at, parent=new_molec, index=kdx, debug=debug)
            atoms.append(new_atom)
        new_molec.set_atoms(atomlist=atoms, debug=debug)
    if debug > 0: print(f"IMPORT MOLEC: Example of imported atom")
    if debug > 0: print(new_molec.atoms[0])

    ## Substructures
    if debug > 0: print(f"IMPORT MOLEC: Importing Substructures")
    if not hasattr(mol,"ligandlist") or not hasattr(mol,"metalist"):
        if debug > 0: print(f"IMPORT MOLEC: splitting complex")
        new_molec.split_complex()
    else:
        # Ligands
        new_molec.ligands = []
        for lig in mol.ligandlist: 
            new_lig = import_ligand(lig, parent=new_molec, debug=debug)
            new_molec.ligands.append(new_lig)
            if debug > 0: print(f"-------")
        # Metals         !! now, metals are taken from molecule.atoms. Otherwise it wasn't working
        new_molec.metals = []
        for at in new_molec.atoms:
            if at.subtype == "metal": new_molec.metals.append(at)
        #for met in mol.metalist: 
        #    new_atom = import_atom(met, parent=new_molec, debug=debug)
        #    new_atom.origin = "import_molecule"
        #    new_molec.metals.append(new_atom)
        #    if debug > 0: print(f"-------")

    ## Charges
    if hasattr(mol,"totcharge") and hasattr(mol,"atcharge"):
        if debug > 0: print(f"IMPORT MOLEC: imported total charge and atomic charges")
        new_molec.set_charges(mol.totcharge, mol.atcharge)
    elif hasattr(mol,"totcharge") and not hasattr(mol,"atcharge"):
        if debug > 0: print(f"IMPORT MOLEC: imported total charge but no atomic charges")
        new_molec.set_charges(mol.totcharge)

    ## Fractional coordinates
    if debug > 0: print(f"IMPORT MOLEC: trying to import fractional coordinates")
    if   hasattr(mol,"frac_coord"):   new_molec.set_fractional_coord(mol.frac_coord, debug=debug)
    elif hasattr(mol,"frac"):         new_molec.set_fractional_coord(mol.frac, debug=debug)
    else:                             new_molec.get_fractional_coord(debug=debug)

    return new_molec

################################
def import_ligand(lig: object, parent: object=None, debug: int=0) -> object:
    assert hasattr(lig,"labels") and (hasattr(lig,"coord") or hasattr(lig,"pos"))
    labels     = lig.labels
    if   hasattr(lig,"coord"):       coord      = lig.coord
    elif hasattr(lig,"pos"):         coord      = lig.pos

    if   hasattr(lig,"parent_indices"):  indices    = lig.parent_indices
    elif hasattr(lig,"atlist"):          indices    = lig.atlist
    else:                                indices    = None

    if   hasattr(lig,"radii"):       radii      = lig.radii
    else:                            radii      = None          

    if debug > 0 and parent is None: print("IMPORT LIGAND: parent is NONE")

    new_ligand = ligand(labels, coord, radii)
    new_ligand.origin = "import_ligand"
    if debug > 0: print(f"IMPORT LIGAND: importing ligand with {new_ligand.formula}")

    ## Parents
    if parent is not None:
        new_ligand.add_parent(parent, indices, overwrite=False, debug=debug)
        if debug > 0: print(f"IMPORT LIGAND: parent {parent.subtype} added with {indices=}") 
    else:
        if debug > 0: print(f"IMPORT LIGAND: parent is None") 

    ## Charges
    if hasattr(lig,"totcharge") and hasattr(lig,"atcharge"):        
        new_ligand.set_charges(lig.totcharge, lig.atcharge)
        if debug > 0: print(f"IMPORT LIGAND: imported total charge and atomic charges")
    elif hasattr(lig,"totcharge") and not hasattr(lig,"atcharge"):  
        new_ligand.set_charges(lig.totcharge)
        if debug > 0: print(f"IMPORT LIGAND: imported total charge but no atomic charges")

    ## Smiles
    if   hasattr(lig,"Smiles"): new_ligand.smiles = lig.Smiles
    elif hasattr(lig,"smiles"): new_ligand.smiles = lig.smiles     
    elif debug > 0: print(f"IMPORT LIGAND: SMILES could not be imported")

    ## Rdkit Object
    if hasattr(lig,"object"): 
        new_ligand.rdkit_obj = lig.object
    else:
        if debug > 0: print(f"IMPORT LIGAND: RDKIT OBJECT could not be imported")
    
    ## Substructures
    if not hasattr(lig,"grouplist"): new_ligand.split_ligand()
    else: 
        new_ligand.groups = []
        for gr in lig.grouplist: 
            new_group = import_group(gr, parent=new_ligand, debug=debug)
            new_ligand.groups.append(new_group)
    
    ## Atoms
    #if not hasattr(lig,"atoms"):  
    if new_ligand.check_parent("molecule"):
        if debug > 0: print(f"IMPORT LIGAND: importing atoms from molecule")
        ## Tries to get them from parent molecule if exists:
        lig_idx_in_mol = new_ligand.get_parent_indices("molecule")
        parent_mol = new_ligand.get_parent("molecule")
        atoms = extract_from_list(lig_idx_in_mol, parent_mol.atoms, dimension=1) 
        new_ligand.set_atoms(atomlist=atoms, debug=debug)
    else:
        ## Otherwise creates them
        if debug > 0: print(f"IMPORT LIGAND: old ligand has no atoms and no parent molecule. Creating new atoms from scratch")
        new_ligand.set_atoms(debug=debug)
        if debug > 0: print(f"IMPORT LIGAND: {len(new_ligand.atoms)} atoms created. Printed below")
        if debug > 0: print(new_ligand.atoms)
    #else: 
    #    if debug > 0: print(f"IMPORT LIGAND: importing atoms from old ligand")
    #    atoms = []
    #    for at in lig.atoms: 
    #        new_atom = import_atom(at, parent=new_ligand, debug=debug)
    #        atoms.append(new_atom)
    #    new_ligand.set_atoms(atomlist=atoms, debug=debug)
            
    return new_ligand

################################
def import_group(old_group: object, parent: object=None, debug: int=0) -> object:

    ## In cell2mol_v1 cells, groups didn't have coordinates nor labels nor nothing
    assert hasattr(old_group,"atlist") or hasattr(old_group,"indices")
    
    ## Labels
    if   hasattr(old_group,"labels"):      labels     = old_group.labels
    elif parent is not None:
        if hasattr(parent,"labels"):       labels     = extract_from_list(old_group.atlist, parent.labels, dimension=1)
        
    ## Coordinates
    if   hasattr(old_group,"coord"):       coord      = old_group.coord
    elif hasattr(old_group,"pos"):         coord      = old_group.pos
    elif parent is not None:
        if hasattr(parent,"coord"):        coord     = extract_from_list(old_group.atlist, parent.coord, dimension=1)
        elif hasattr(parent,"coord"):      coord     = extract_from_list(old_group.atlist, parent.pos, dimension=1)

    if   hasattr(old_group,"parent_indices"):     indices    = old_group.parent_indices
    elif hasattr(old_group,"atlist"):             indices    = old_group.atlist
    else:                                         indices    = None

    ## Radii
    if   hasattr(old_group,"radii"):       radii  = old_group.radii
    else:                                  radii  = None          

    new_group = group(labels, coord, radii)
    new_group.origin = "import_group"
    if debug > 0: print(f"IMPORT GROUP: importing group with {new_group.formula}")

    ## Parents
    if parent is not None:
        new_group.add_parent(parent, indices, overwrite=False, debug=debug)
        if debug > 0: print(f"IMPORT GROUP: parent {parent.subtype} added with {indices=}") 
    else:
        if debug > 0: print(f"IMPORT GROUP: parent is None") 

    ## Charges
    if hasattr(old_group,"totcharge") and hasattr(old_group,"atcharge"):
        new_group.set_charges(old_group.totcharge, old_group.atcharge)
        if debug > 0: print(f"IMPORT GROUP: imported total charge and atomic charges")
    elif hasattr(old_group,"totcharge") and not hasattr(old_group,"atcharge"):
        new_group.set_charges(old_group.totcharge)
        if debug > 0: print(f"IMPORT GROUP: imported total charge but no atomic charges")

    ## Atoms
    #if not hasattr(old_group,"atoms"):  
    if new_group.check_parent("molecule"):
        if debug > 0: print(f"IMPORT GROUP: importing atoms from molecule")
        ## Tries to get them from parent molecule if exists:
        group_idx_in_mol = new_group.get_parent_indices("molecule")
        parent_mol = new_group.get_parent("molecule")
        atoms = extract_from_list(group_idx_in_mol, parent_mol.atoms, dimension=1) 
        new_group.set_atoms(atomlist=atoms, debug=debug)
    else:
        ## Otherwise creates them
        if debug > 0: print(f"IMPORT GROUP: old group has no atoms and no parent molecule. Creating new atoms from scratch")
        new_group.set_atoms(debug=debug)
        if debug > 0: print(f"IMPORT GROUP: {len(new_group.atoms)} atoms created. Printed below")
        if debug > 0: print(new_group.atoms)
            
    return new_group

################################
def import_atom(old_atom: object, parent: object=None, index: int=None, debug: int=0) -> object:
    assert hasattr(old_atom,"label") and (hasattr(old_atom,"coord") or hasattr(old_atom,"pos"))
    label     = old_atom.label
    
    if   hasattr(old_atom,"coord"):      coord      = old_atom.coord
    elif hasattr(old_atom,"pos"):        coord      = old_atom.pos

    if index is None:
        if   hasattr(old_atom,"parent_index"): index    = old_atom.parent_index
        elif hasattr(old_atom,"atlist"):       index    = old_atom.atlist
        elif hasattr(old_atom,"index"):        index    = old_atom.index
        else:                                  index    = None

    if   hasattr(old_atom,"radii"):      radii      = old_atom.radii
    else:                                radii      = None          

    if   hasattr(old_atom,"block"):      block      = old_atom.block
    else:                                block      = elemdatabase.elementblock[label]    

    if block == 'd' or block == 'f': 
        if debug > 0: print(f"IMPORT ATOM: importing metal {old_atom.label}")
        new_atom = metal(label, coord, radii=radii)
    else:                            
        if debug > 0: print(f"IMPORT ATOM: importing atom {old_atom.label}")
        new_atom = atom(label, coord, radii=radii)
    new_atom.origin = "import_atom"
    
    ## Parents
    if parent is not None:
        new_atom.add_parent(parent, index, overwrite=False, debug=debug)
        if debug > 0: print(f"IMPORT ATOM: parent {parent.subtype=} added with {index=}") 
    else:
        if debug > 0: print(f"IMPORT ATOM: parent is None") 
    
    ## Charge. For some weird reason, charge is "charge" in regular atoms, and "totcharge" in metals
    ## Also. In cell2mol version 1, atoms do not carry charge. Metals in metalist yes
    if debug > 0: print(f"IMPORT ATOM: now trying to import charge") 
    if hasattr(old_atom,"charge"):       
        if type(old_atom.charge) == int:    
            new_atom.set_charge(old_atom.charge)
            if new_atom.subtype == "atom"  and debug > 0: print(f"IMPORT ATOM: atom charge set to {new_atom.charge}") 
            if new_atom.subtype == "metal" and debug > 0: print(f"IMPORT ATOM: metal charge set to {new_atom.charge}") 

    elif hasattr(old_atom,"totcharge"):     
        if type(old_atom.totcharge) == int: 
            new_atom.set_charge(old_atom.totcharge)
            if new_atom.subtype == "atom"  and debug > 0: print(f"IMPORT ATOM: atom charge set to {new_atom.charge}") 
            if new_atom.subtype == "metal" and debug > 0: print(f"IMPORT ATOM: metal charge set to {new_atom.charge}") 
    else:
        if debug > 0: print(f"IMPORT ATOM: no charge found for old atom {old_atom.label}")

    ## Connectivity
    if debug > 0: print(f"IMPORT ATOM: inheriting connectivity for {new_atom.subtype=}")
    if index is not None: new_atom.inherit_connectivity("molecule", debug=debug)
    else:         print(f"IMPORT ATOM: connectivity could not be imported since Index=None")

    ## Factors for connectivity calculation
    ## In cell2mol version 1, atoms do not carry factors. Metals do
    if   hasattr(old_atom,"metal_factor") and hasattr(old_atom,"factor"): 
        cov_factor   = old_atom.factor
        metal_factor = old_atom.metal_factor
    elif hasattr(old_atom,"metal_factor") and hasattr(old_atom,"cov_factor"): 
        cov_factor   = old_atom.cov_factor
        metal_factor = old_atom.metal_factor
    elif new_atom.check_parent("molecule"):
        par = new_atom.get_parent("molecule")
        if   hasattr(par,"metal_factor") and hasattr(par,"factor"):
            cov_factor   = par.factor
            metal_factor = par.metal_factor
        elif hasattr(par,"metal_factor") and hasattr(par,"cov_factor"):
            cov_factor   = par.cov_factor
            metal_factor = par.metal_factor
        else:
            cov_factor   = 1.3 
            metal_factor = 1.0 
    else:
        cov_factor   = 1.3 
        metal_factor = 1.0 
    new_atom.set_factors(cov_factor, metal_factor)
            
    return new_atom
