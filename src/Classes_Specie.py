####################################################
####  Contains the SPECIE Class and Sub-Classes ####
####################################################

import numpy as np
from .Connectivity import * 
from .Classes_Atom import *
from .Other import get_metal_idxs, get_dist
from .Unit_cell_tools import * 
from .Elementdata import ElementData
elemdatabase = ElementData()

##############
### SPECIE ###
##############
class specie(object):
    def __init__(self, labels: list, coord: list, frac_coord: list=None, radii: list=None) -> None:

       # Sanity Checks
        assert len(labels) == len(coord)
        if frac_coord is not None:        
            assert len(coord) == len(frac_coord)
            self.frac_coord = frac_coord
            
        # Optional Information
        if radii is not None:   self.radii   = radii                  ## Radii are used to obtain the adjacency matrix 
        else:                   self.radii   = get_radii(labels)

        self.version              = "1.0"
        self.type                 = "specie"
        self.subtype              = "specie"
        self.origin               = "created"
        self.labels               = labels
        self.coord                = coord
        self.formula              = labels2formula(labels)
        self.eleccount            = labels2electrons(labels)
        self.natoms               = len(labels)
        self.iscomplex            = any((elemdatabase.elementblock[l] == "d") or (elemdatabase.elementblock[l] == "f") for l in self.labels)
        self.indices              = [*range(0,self.natoms,1)]

        ## Specie can be associated with "parents", which are other species related to it.
        ## For instance, a ligand can be related to the transition metal complex it belongs to.
        self.parents              = []
        self.parents_indices      = []

        ## Defaults. Were used in cell2mol version 1, kept for simplicity
        self.cov_factor           = 1.3
        self.metal_factor         = 1.0
        
        ## Bonds
        self.has_bonds            = False

    ######
    def __repr__(self, indirect: bool=False):
        to_print                   = ''
        if not indirect: to_print += '------------------------------------------------\n'
        if not indirect: to_print += '------------- SCOPE SPECIE Object --------------\n'
        if not indirect: to_print += '------------------------------------------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        if hasattr(self,'subtype'):        to_print += f' Sub-Type              = {self.subtype}\n'
        to_print += f' Number of Atoms       = {self.natoms}\n'
        to_print += f' Formula               = {self.formula}\n'
        if hasattr(self,"parents"):        to_print += f' Number of Parents     = {len(self.parents)}\n'
        if hasattr(self,"totcharge"):      to_print += f' Total Charge          = {self.totcharge}\n'
        if hasattr(self,"spin"):           to_print += f' Spin                  = {self.spin}\n'
        if hasattr(self,"smiles"):         to_print += f' SMILES                = {self.smiles}\n'
        if hasattr(self,"adjmat"):         to_print += f' Has Adjacency Matrix  = YES\n'
        else:                              to_print += f' Has Adjacency Matrix  = NO \n'
        if not hasattr(self,"has_bonds"):  to_print += f' Has Bonds             = NO\n'
        else:
            if self.has_bonds:             to_print += f' Has Bonds             = YES\n'
            else:                          to_print += f' Has Bonds             = NO\n'
        if not indirect: to_print += '\n'
        return to_print

    ######
    def save(self, filepath):
        from .Read_Write import save_binary
        save_binary(self, filepath)

    ######
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
                            self.parents_indices[idx] = new_indices
                        append = False
                if append: 
                    self.parents.append(p2)
                    self.parents_indices.append(new_indices)
                    if debug > 0: print(f"SPECIE.ADD_PARENT: added parent of parent with subtype {p2.subtype}. Indices are {new_indices}")

    ######
    def check_parent(self, subtype: str):
        ## checks if parent of a given subtype exists
        for p in self.parents:
            if p.subtype == subtype: return True
        return False

    ######
    def get_parent(self, subtype: str):
        ## retrieves parent of a given subtype 
        for p in self.parents:
            if hasattr(p, "subtype"):  
                if p.subtype.lower() == subtype.lower(): return p
            else:  
                print(f"Warning. Parent with type: {p.type} does not have subtype")
        return None

    ######
    def get_parent_indices(self, subtype: str):
        ## retrieves parent indices of a given subtype 
        for idx, p in enumerate(self.parents):
            if hasattr(p, "subtype"):  
                if p.subtype.lower() == subtype.lower(): return self.parents_indices[idx]
        return None

    ######
    def get_centroid(self):
        self.centroid = compute_centroid(self.coord)
        if hasattr(self,"frac_coord"): self.frac_centroid = compute_centroid(self.frac_coord)
        return self.centroid

    ######
    def set_fractional_coord(self, frac_coord: list, debug: int=0) -> None:
        if debug > 0: print(f"SPECIE.SET_FRACTIONAL_COORD: length comparison: {len(frac_coord)} vs {len(self.coord)}")
        if debug > 0: print(f"SPECIE.SET_FRACTIONAL_COORD: length comparison: {np.shape(frac_coord)} vs {np.shape(self.coord)}")
        if debug > 0: print(f"SPECIE.SET_FRACTIONAL_COORD: 1st items: {frac_coord[0]} vs {self.coord[0]}")
        #assert len(frac_coord) == len(self.coord)  ## Removed 'cos it was complaining despite both vectors having the same length
        self.frac_coord = frac_coord 

    ######
    def get_fractional_coord(self, cell_vector=None, debug: int=0) -> None:
        from .Reconstruct import cart2frac

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

    ######
    def get_atomic_numbers(self):
        if not hasattr(self,"atoms"): self.set_atoms()
        self.atnums = []
        for at in self.atoms:
            self.atnums.append(at.atnum)
        return self.atnums

    ######
    def set_element_count(self, heavy_only: bool=False):
        self.element_count = get_element_count(self.labels, heavy_only=heavy_only)
        return self.element_count

    ######
    def set_adj_types(self):
        if not hasattr(self,"adjmat"): self.get_adjmatrix()
        self.adj_types = get_adjacency_types(self.labels, self.adjmat)
        return self.adj_types

    ######
    def reset_charge(self):
        if hasattr(self,"totcharge"):      delattr(self,"totcharge")
        if hasattr(self,"atomic_charges"): delattr(self,"atomic")
        if hasattr(self,"smiles"):         delattr(self,"smiles")
        if hasattr(self,"rdkit_obj"):      delattr(self,"rdkit_obj")
        if hasattr(self,"poscharges"):     delattr(self,"poscharges") 
        for a in self.atoms:               
            a.reset_charge() 

    ######
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

    ######
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

                #################################
                # Prepares fractional coordinates
                #################################
                # First tries to get it from parent cell, if exists
                if not hasattr(self,"frac_coord"): 
                    if self.check_parent("cell"):
                       par = self.get_parent("cell")
                       if hasattr(par,"cell_vector"): self.get_fractional_coord(debug=debug)
                # If it managed, then it established the frac_coord to atoms 
                if hasattr(self,"frac_coord"):
                    if ismetal: newatom = metal(l, self.coord[idx], self.frac_coord[idx], radii=self.radii[idx])
                    else:       newatom =  atom(l, self.coord[idx], self.frac_coord[idx],radii=self.radii[idx])
                # Otherwise, frac_coord is empty
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

    ######
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

    ######
    def set_bonds(self, debug: int=0):
        ## Creats bond objects using the information contained in the RDKit object
        ## The RDKit object is necessary, as it is the only one that contains the bond order information (lewis structure)
        if not hasattr(self,"rdkit_obj"): 
            if debug > 0: print(f"SPECIE.SET_BONDS: Can't set bonds, specie lacks an rdkit_object: {self.formula}")
            return False
        natoms_rdkit = self.rdkit_obj.GetNumAtoms() 
        if debug >= 1: print(f"SPECIE.SET_BONDS: {self.formula=}, {self.subtype=} {self.smiles=}")

        if self.natoms == natoms_rdkit:  
            if debug >= 2: print(f"\tNumber of atoms in {self.subtype} object and RDKit object are equal: {self.natoms} {natoms_rdkit}")

            for idx, rdkit_atom in enumerate(self.rdkit_obj.GetAtoms()): # e.g. idx 0, 1, 2, 3, 4, 5, 6, 7, 8
                if debug >= 2: print(f"\t{idx=}", rdkit_atom.GetSymbol(), "Number of bonds in rdkit_obj:", len(rdkit_atom.GetBonds()))
                if len(rdkit_atom.GetBonds()) == 0:
                    if debug >= 1: print(f"\tNO BONDS CREATED for {self.atoms[idx].label} due to no bonds in {self.subtype} RDKit object")
                else:
                    for b in rdkit_atom.GetBonds():
                        bond_startatom = b.GetBeginAtomIdx()
                        bond_endatom   = b.GetEndAtomIdx()
                        bond_order     = b.GetBondTypeAsDouble()

                        ## Checks the labels involved in the bond
                        if self.atoms[bond_endatom].label != self.rdkit_obj.GetAtomWithIdx(bond_endatom).GetSymbol():
                            if debug >= 1: print(f"\tError with Bond EndAtom", self.atoms[bond_endatom].label, self.rdkit_obj.GetAtomWithIdx(bond_endatom).GetSymbol())
                        else:
                            if bond_endatom == idx:
                                start = bond_endatom
                                end   = bond_startatom
                            elif bond_startatom == idx:
                                start = bond_startatom
                                end   = bond_endatom      

                            # create new bond object and add it to both atoms
                            if end > start:
                                if debug >=2: print(f"\tBOND CREATED", idx, start, end, bond_order, self.atoms[start].label, self.atoms[end].label)
                                new_bond = bond(self.atoms[start], self.atoms[end], bond_order)
                                self.atoms[start].add_bond(new_bond)
                                self.atoms[end].add_bond(new_bond)
                
                    if hasattr(self.atoms[idx], "bonds"):
                        if debug >=1 : print(f"\tBONDS", [(bd.atom1.label, bd.atom2.label, bd.order, round(bd.distance,3)) for bd in self.atoms[idx].bonds])
                    else:
                        if self.natoms == 1:
                            if debug >=1: print(f"\tNO BONDS CREATED for {self.atoms[idx].label} because it is the only atom in {self.subtype} object")
                            pass
                        else:
                            if debug >=1: print(f"\tNO BONDS for {self.atoms[idx].label} with {self.subtype} RDKit object index {idx}. Please check the RDKit object.")
                            return False # return False if no bonds are created
        else:
            if debug >= 1: print(f"\tNumber of atoms in {self.subtype} object and RDKit object are different: {self.natoms} {natoms_rdkit}")
            if debug >= 2: print(f"\t{[(i, atom.label) for i, atom in enumerate(self.atoms)]}")
            if debug >= 2: print(f"\t{[(i, atom.GetSymbol()) for i, atom in enumerate(self.rdkit_obj.GetAtoms())]}")       
            non_bonded_atoms = list(range(0, natoms_rdkit))[self.natoms:]
            if debug >= 2: print(f"\tNON_BONDED_ATOMS", non_bonded_atoms)

            for idx, rdkit_atom in enumerate(self.rdkit_obj.GetAtoms()): # e.g. idx 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
                if debug >= 2: print(f"\t{idx=}", rdkit_atom.GetSymbol(), "Number of bonds :", len(rdkit_atom.GetBonds()))
                if len(rdkit_atom.GetBonds()) == 0:
                    if debug >= 1: print(f"\tNO BONDS CREATED for {rdkit_atom.GetSymbol()} due to no bonds in {self.subtype} RDKit object")
                else:
                    for b in rdkit_atom.GetBonds():
                        bond_startatom = b.GetBeginAtomIdx()
                        bond_endatom   = b.GetEndAtomIdx()
                        bond_order     = b.GetBondTypeAsDouble()
  
                        if bond_startatom in non_bonded_atoms or bond_endatom in non_bonded_atoms:
                            if debug >= 2: print(f"\tNO BOND CREATED {bond_startatom=} or {bond_endatom=} is not in the self.atoms. It belongs to {non_bonded_atoms=}.")
                        else :
                            if bond_endatom == idx:
                                start = bond_endatom
                                end   = bond_startatom
                            elif bond_startatom == idx:
                                start = bond_startatom
                                end   = bond_endatom   

                            # create new bond object and add it to both atoms
                            if end > start:
                                if debug >=2: print(f"\tBOND CREATED", idx, start, end, bond_order, self.atoms[start].label, self.atoms[end].label)
                                new_bond = bond(self.atoms[start], self.atoms[end], bond_order)
                                self.atoms[start].add_bond(new_bond)
                                self.atoms[end].add_bond(new_bond)
                
                    if idx not in non_bonded_atoms:
                        if hasattr(self.atoms[idx], "bonds"):
                            if debug >=2: 
                                print(f"\tBONDS", [(bd.atom1.label, bd.atom2.label, bd.order, round(bd.distance,3)) for bd in self.atoms[idx].bonds])
                        else:
                            if self.natoms == 1:
                                if debug >=1: print(f"\tNO BONDS CREATED for {self.atoms[idx].label} because it is the only atom in {self.subtype} object")
                                pass
                            else:
                                if debug >=1: print(f"\tNO BONDS for {self.atoms[idx].label} with {self.subtype} RDKit object index {idx}. Please check the RDKit object.")
                                return False # return False if no bonds are created
                    else :
                        if debug >=1: print(f"\tNO BONDS for {rdkit_atom.GetSymbol()} with {self.subtype} RDKit object index {idx} because it is an added atom")
        self.has_bonds = True
        return True 

    ######
    def set_factors(self, cov_factor: float=1.3, metal_factor: float=1.0) -> None:
        self.cov_factor   = cov_factor
        self.metal_factor = metal_factor
        if hasattr(self,"atoms"):
            for at in self.atoms:
                at.set_factors(cov_factor=self.cov_factor, metal_factor=self.metal_factor)

    ######
    def get_adjmatrix(self):
        isgood, adjmat, adjnum = get_adjmatrix(self.labels, self.coord, self.cov_factor, self.metal_factor, radii=self.radii)
        if isgood:
            self.adjmat = adjmat
            self.adjnum = adjnum
        else:
            self.adjmat = None
            self.adjnum = None
        return self.adjmat, self.adjnum

    ######
    def get_metal_adjmatrix(self):
        isgood, madjmat, madjnum = get_adjmatrix(self.labels, self.coord, self.cov_factor, self.metal_factor, radii=self.radii, metal_only=True)
        if isgood:
            self.madjmat = madjmat
            self.madjnum = madjnum
        else:
            self.madjmat = None
            self.madjnum = None
        return self.madjmat, self.madjnum

    ######
    def get_occurrence(self, substructure: object, debug: int=0) -> int:
        """
        Counts the number of times a given substructure appears within the current object.
        The method supports different types of substructures:
        - Ligands within Complexes (when substructure.subtype == 'ligand' and self.subtype == 'molecule')
        - Groups within Ligands    (when substructure.subtype == 'group' and self.subtype == 'ligand')
        - Atoms within Species     (when substructure.type == 'atom' and self.type == 'specie')
        The method uses helper functions such as `compare_species` and `compare_atoms` to determine equivalence,
        and may call `split_complex`, `split_ligand`, or `set_atoms` to ensure the relevant attributes are initialized.
        Args:
            substructure (object): The substructure to search for within the current object.
            debug (int, optional): Debug level for comparison functions. Defaults to 0.
        Returns:
            int: The number of times the substructure appears within the current object.
        """
    
        ## Finds how many times a substructure appears in self
        occurrence = 0

        ## Case of Ligands in Complexes or Groups in Ligands
        done = False
        if hasattr(substructure,"subtype") and hasattr(self,"subtype"):
            if substructure.subtype == 'ligand' and self.subtype == 'molecule':
                if not hasattr(self,"ligands"): self.split_complex()
                if self.ligands is not None:
                    for l in self.ligands:
                        issame = compare_species(substructure, l, debug=debug)
                        if issame: occurrence += 1
                    done = True 
            elif substructure.subtype == 'group' and self.subtype == 'ligand':
                if not hasattr(self,"ligands"): self.split_complex()
                if self.ligands is not None:
                    for l in self.ligands:
                        if not hasattr(l,"groups"): self.split_ligand()
                        for g in l.groups:
                            issame = compare_species(substructure, g, debug=debug)
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

    ######
    def set_spin(self, spin: int) -> None:
        self.spin = spin

    ######
    def check_fragmentation(self, cov_factor: float=1.3, metal_factor: float=None, debug: int=0):
        blocklist = split_species(self.labels, self.coord, cov_factor=cov_factor)
        if len(blocklist) > 1: self.isfragmented = True
        else:                  self.isfragmented = False
        return self.isfragmented

    #########################################
    ### Functions to Interact with States ###
    #########################################
    def add_state(self, name: object, debug: int=0):
        from .Classes_State import state
        if not hasattr(self,"states"): setattr(self,"states",list([]))
        exists, new_state = self.find_state(name)
        if exists:  
            if debug > 0: print(f"SPECIE.ADD_STATE. State with same {name=} found, returning it")
            return new_state
        else:
            if debug > 0: print("SPECIE.ADD_STATE. Creating new state, returning it")
            new_state = state(self, name, debug=debug)
            self.states.append(new_state)
        return new_state

    ######
    def find_state(self, search_name: str, debug: int=0):
        from .Classes_State import state
        if not hasattr(self,"states"): setattr(self,"states",list([]))
        if debug > 0: print(f"SPECIE.FIND_STATE: Searching {search_name} in SPECIE object with {len(self.states)} states")
        for sta in self.states:
            if debug > 0: print(f"SPECIE.FIND_STATE: Comparing {search_name} with {sta.name}")
            if sta.name == search_name: 
                if debug > 0: print(f"SPECIE.FIND_STATE: state {search_name} found")
                return True, sta
        if debug > 0: print(f"SPECIE.FIND_STATE: state {search_name} not found")
        return False, None

    ######
    def rmsd(self, other, reorder=True, center_method='centroid', debug: int=0):
        from .Other import rmsd
        value = rmsd(self.labels, self.coord, other.labels, other.coord, reorder=reorder, center_method=center_method, debug=debug)   
        return value

    ###################################
    ### Functions to print/visualize ##
    ###################################
    def print_xyz(self):
        print(self.natoms)
        print("")
        for idx, l in enumerate(self.labels):
            print("%s  %.6f  %.6f  %.6f" % (l, self.coord[idx][0], self.coord[idx][1], self.coord[idx][2]))

    ######
    def view(self, show_indices: bool=False, size: str='default'):
        import plotly.graph_objects as go
        from .Read_Write import set_scene
        from .Elementdata import ElementData  
        elemdatabase = ElementData()

        size_map = {'default': (600, 600, 8, 9), 'small': (400, 400, 6, 7), 
                    'large': (800, 800, 10, 12), 'ultra': (1000, 1000, 11, 13)}
        width, height, marker_size, text_size = size_map.get(size.lower(), size_map['default'])

        if not hasattr(self, "adjmat"): self.get_adjmatrix()
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
            mode='text', text=[str(i) for i in range(len(positions))] if show_indices else symbols,
            textfont=dict(color='black', size=text_size), hoverinfo='none', showlegend=False))

        for i, j in unique_bonds:
            fig.add_trace(go.Scatter3d(
                x=[positions[i, 0], positions[j, 0]], y=[positions[i, 1], positions[j, 1]],
                z=[positions[i, 2], positions[j, 2]], mode='lines',
                line=dict(color='gray', width=5), hoverinfo='none', showlegend=False))

        set_scene(fig, positions, width=width, height=height)
        fig.show()

    ###########################
    ### Other Dunder Methods ##
    ###########################
    def __add__(self, other):
        ## To be implemented
        if not isinstance(other, type(self)): return self
        return self

    def __len__(self):
        return self.natoms

    def __eq__(self, other):
        elems = elemdatabase.elementnr.keys()
        
        # a pair of species is compared on the basis of:
        # 1) the total number of atoms
        if (self.natoms != other.natoms): 
            #print(f"COMPARE_SPECIES. FALSE: {self.natoms=} vs. {other.natoms=}")
            return False

        # 2) the total number of electrons (as sum of atomic number)
        if (self.eleccount != other.eleccount): 
            #print(f"COMPARE_SPECIES. FALSE: {self.eleccount=} vs. {other.eleccount=}")
            return False

        # 3) the number of atoms of each type
        if not hasattr(self,"element_count"): self.set_element_count()
        if not hasattr(other,"element_count"): other.set_element_count()
        for kdx, elem in enumerate(self.element_count):
            if elem != other.element_count[kdx]: 
                #print(f"COMPARE_SPECIES. FALSE, different {elem} count:")
                return False       
        # 4) the number of adjacencies between each pair of element types
        if not hasattr(self,"adj_types"):     self.set_adj_types()
        if not hasattr(other,"adj_types"):     other.set_adj_types()

        count = 0
        for kdx, (elem, row1) in enumerate(zip(elems, self.adj_types)):
            for ldx, (elem2, val1) in enumerate(zip(elems, row1)):
                val2 = other.adj_types[kdx, ldx]
                if val1 != val2: 
                    count += 1
        if count > 0 : return False
        return True

###############
### MOLECULE ##
###############
class molecule(specie):
    def __init__(self, labels: list, coord: list, frac_coord: list=None, radii: list=None) -> None:
        specie.__init__(self, labels, coord, frac_coord, radii)
        self.subtype = "molecule"

    def __repr__(self, indirect: bool=False):
        to_print                   = ''
        if not indirect: to_print += '--------------------------------------------------\n'
        if not indirect: to_print += '------------- SCOPE MOLECULE Object --------------\n'
        if not indirect: to_print += '--------------------------------------------------\n'
        to_print += specie.__repr__(self, indirect=True)
        if hasattr(self,"ligands"):  
            if self.ligands is not None: to_print += f' # Ligands             = {len(self.ligands)}\n'
        if hasattr(self,"metals"):   
            if self.metals is not None:  to_print += f' # Metals              = {len(self.metals)}\n'
        if not indirect: to_print += '\n'
        return to_print

    ######
    def reset_charge(self):
        specie.reset_charge(self)      ## First, uses the generic specie-class function for itself and its atoms
        if hasattr(self,"ligands"):    ## Second, resets the child classes 
            for lig in self.ligands:
                lig.reset_charge()
        if hasattr(self,"metals"):    
            for met in self.metals:
                met.reset_charge()

    ######
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
                self.metals.append(self.atoms[m])                            
        return self.ligands, self.metals
        
    ######
    def fix_ligands_rdkit_obj(self, debug: int=0):
        if debug > 0 and self.iscomplex:
            print(f"SPECIE.FIX_LIGANDS_RDKIT_OBJ. Fixing specie {self.formula}")
            print(f"    With ligand formula and smiles:")
            for lig in self.ligands:
                print(f"        {lig.formula} {lig.smiles}")
        if self.iscomplex and not hasattr(self,"ligands"): self.split_complex(debug=debug)
        if not self.iscomplex: return self.smiles
        self.smiles = []
        for lig in self.ligands:
            if not hasattr(lig,"smiles"): 
                raise ValueError(self, lig) 
            lig.fix_rdkit_obj(debug=debug)
            lig.set_smiles_from_rdkit_obj(debug=debug) 
            self.smiles.append(lig.smiles)
        return self.smiles

    ######
    def get_cshm(self, ref_shape: str='OC-6', overwrite: bool=False, debug: int=0):
        if not hasattr(self,"metals"): self.split_complex(debug=debug)
        self.cshm = []
        for met in self.metals:
            self.cshm.append(met.get_cshm(ref_shape=ref_shape, overwrite=overwrite, debug=debug))
        return self.cshm

    ######
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

    ######
    def set_bonds(self, debug: int=0):
        # If it is a regular molecule, just use the specie function
        if not self.iscomplex: return super().set_bonds(debug=debug)
        # If it is a transition metal complex, then:
        if not hasattr(self,"metals") or not hasattr(self,"ligands"): self.split_complex(debug=debug)
        for lig in self.ligands:            ## Runs for each ligand individually
            lig.set_bonds(debug=debug)      
        self.set_metal_ligand_bonds(debug=debug) ## Creates bond between metals and ligands...
        self.set_metal_metal_bonds(debug=debug)  ## ... and between metals 
        self.has_bonds = True
        return True

    ######
    def set_metal_ligand_bonds(self, debug: int=0):
        if not self.iscomplex: return False
        if not hasattr(self,"metals") or not hasattr(self,"ligands"): self.split_complex(debug=debug)
        for lig in self.ligands:
            for at in lig.atoms:
                count = 0
                for met in self.metals: 
                    isconnected = at.check_connectivity(met, debug=debug)
                    if isconnected:
                        index_1 = at.get_parent_index("molecule")
                        index_2 = met.get_parent_index("molecule")
                        if index_1 < index_2 : 
                            bond_startatom = at
                            bond_endatom   = met
                        else:
                            bond_startatom = met
                            bond_endatom   = at
                        newbond = bond(bond_startatom, bond_endatom, 0.5, subtype="metal-ligand")

                        at.add_bond(newbond)
                        met.add_bond(newbond)
                        count += 1 
                if hasattr(at, "madjnum"): 
                    if count != at.madjnum: 
                        if debug >= 1: print(f"MOL.CREATE_BONDS: error creating bonds for atom: \n{at}\n of ligand: \n{lig}\n")
                        if debug >= 1: print(f"MOL.CREATE_BONDS: count differs from atom.mconnec: {count}, {at.madjnum}")
        return True

    ######
    def set_metal_metal_bonds(self, debug: int=0):
        if not self.iscomplex: return False
        if not hasattr(self,"metals") or not hasattr(self,"ligands"): self.split_complex(debug=debug)
        if len(self.metals) < 2: return False
        if debug >= 1: print(f"MOL.CREATE_BONDS: Creating Metal-Metal Bonds for selfecule {self.formula}")
        if debug >= 2: print(f"MOL.CREATE_BONDS: Metals: {self.metals}")
        for idx, met1 in enumerate(self.metals):
            for jdx, met2 in enumerate(self.metals):
                if jdx > idx: 
                    if met1.check_connectivity(met2, debug=debug):
                        index_1 = met1.get_parent_index("molecule")
                        index_2 = met2.get_parent_index("molecule")
                        if index_1 < index_2 : 
                            bond_startatom = met1
                            bond_endatom   = met2
                        else:
                            bond_startatom = met2
                            bond_endatom   = met1
                        newbond = bond(bond_startatom, bond_endatom, 0, subtype="metal-metal")
                        met1.add_bond(newbond) 
                        met2.add_bond(newbond) 
        return True

###############
### LIGAND ####
###############
class ligand(specie):
    def __init__(self, labels: list, coord: list, frac_coord: list=None, radii: list=None) -> None:
        if frac_coord is not None:        
            assert len(frac_coord) == len(coord)
            self.frac_coord = frac_coord
        specie.__init__(self, labels, coord, frac_coord, radii)
        self.subtype  = "ligand"

    ######
    def __repr__(self, indirect: bool=False):
        to_print                   = ''
        if not indirect: to_print += '------------------------------------------------\n'
        if not indirect: to_print += '------------- SCOPE LIGAND Object --------------\n'
        if not indirect: to_print += '------------------------------------------------\n'
        to_print += specie.__repr__(self, indirect=True)
        if hasattr(self,"rdkit_obj"):  to_print += f' # HAS RDKIT OBJECT    = YES\n'
        else:                          to_print += f' # HAS RDKIT OBJECT    = NO\n'
        if not indirect: to_print += '\n'
        return to_print

    ######
    def set_hapticity(self, hapttype):
        self.hapticity = True 
        self.hapttype  = hapttype 

    ######
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
    
    ######
    def get_connected_idx(self, debug: int=0):
        ## Remember madjmat should not be computed at the ligand level. Since the metal is not there.
        ## Now we operate at the molecular level. We get the parent molecule, and the indices of the ligand atoms in the molecule
        self.connected_idx = [] 
        if not hasattr(self,"madjnum"): self.inherit_adjmatrix("molecule")
        for idx, con in enumerate(self.madjnum):
            if con > 0: self.connected_idx.append(idx)
        return self.connected_idx 

    ######
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

    ######
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
        connatoms      = extract_from_list(connected_idx, self.atoms, dimension=1)
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
            gr_atoms        = extract_from_list(b, connatoms, dimension=1)
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

    ######
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

    ######
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

    ######
    def fix_rdkit_obj(self, debug: int=0):
        """
        Fixes common issues in the RDKit molecule object (`rdkit_obj`) associated with the ligand, particularly those arising from conversion tools like cell2mol.
        The function applies corrections to specific bonding and charge patterns that are known to be problematic, including:
            - N2-CSe+ ([N--]C#[Se+])
            - N-CSe+  ([N-]C#[Se+])
            - N2-CS+  ([N--]C#[S+])
            - N-CS+   ([N-]C#[S+])
        For each pattern, it adjusts bond orders and formal charges in both the RDKit molecule and the ligand's atom objects.
        Additional corrections include:
            - Removal of extra atoms added during conversion, with charge compensation for neighboring atoms.
            - Correction of adjacent atoms with opposite formal charges, except for nitrosyl groups.
            - Attempts to sanitize the molecule using RDKit's `Chem.SanitizeMol`. If sanitization fails (e.g., due to nitrosyl groups), the original molecule is preserved and the failed attempt is stored in `self.failed_rdkit_obj`.
        Args:
            debug (int, optional): If greater than 0, prints debug information during the fixing process. Default is 0.
        Returns:
            rdkit.Chem.Mol: The fixed RDKit molecule object. If fixing or sanitization fails, returns the original `rdkit_obj`.
        """

        from rdkit import Chem
        if not hasattr(self,"rdkit_obj"): return None

        ## Initializes the rdkit object as mol
        mol = self.rdkit_obj

        #########################
        # Fixes N2-CSe+ pattern # 
        #########################
        pattern = Chem.MolFromSmarts("[N--]C#[Se+]")
        # Check if the substructure exists in the molecule
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                atom_indices = {mol.GetAtomWithIdx(idx).GetAtomicNum(): idx for idx in match}
                n_idx = atom_indices[7]    # Atomic number of nitrogen
                c_idx = atom_indices[6]    # Atomic number of carbon
                se_idx = atom_indices[34]  # Atomic number of selenium

                # Get the editable version of the molecule
                rw_mol = Chem.RWMol(mol)

                # Replace the triple bond between C and Se with a double bond
                rw_mol.RemoveBond(c_idx, se_idx)
                rw_mol.AddBond(c_idx, se_idx, Chem.BondType.DOUBLE)

                # Replace the single bond between N and C with a double bond
                rw_mol.RemoveBond(n_idx, c_idx)
                rw_mol.AddBond(n_idx, c_idx, Chem.BondType.DOUBLE)

                if hasattr(self,"bonds"):
                    # Also in the atoms object, using ligand indices
                    bond = self.atoms[c_idx].find_bond(se_idx, "ligand")
                    if bond is not None: bond.order = 2.0
                    bond = self.atoms[n_idx].find_bond(c_idx, "ligand")
                    if bond is not None: bond.order = 2.0

                # Get Current Charges
                fcharge_N  = rw_mol.GetAtomWithIdx(n_idx).GetFormalCharge()
                fcharge_Se = rw_mol.GetAtomWithIdx(se_idx).GetFormalCharge()

                # Update formal charges
                rw_mol.GetAtomWithIdx(n_idx).SetFormalCharge(fcharge_N+1)   
                rw_mol.GetAtomWithIdx(c_idx).SetFormalCharge(0)             
                rw_mol.GetAtomWithIdx(se_idx).SetFormalCharge(fcharge_Se-1) 

                # Also in ligand object
                self.atoms[n_idx].charge  = fcharge_N+1
                self.atoms[se_idx].charge = fcharge_Se-1

                # Update the molecule
                mol = rw_mol.GetMol()

        #########################
        # Fixes N-CSe+ pattern # 
        #########################
        pattern = Chem.MolFromSmarts("[N-]C#[Se+]")

        # Check if the substructure exists in the molecule
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                atom_indices = {mol.GetAtomWithIdx(idx).GetAtomicNum(): idx for idx in match}
                n_idx = atom_indices[7]    # Atomic number of nitrogen
                c_idx = atom_indices[6]    # Atomic number of carbon
                se_idx = atom_indices[34]  # Atomic number of selenium

                # Get the editable version of the molecule
                rw_mol = Chem.RWMol(mol)

                # Replace the triple bond between C and Se with a double bond
                rw_mol.RemoveBond(c_idx, se_idx)
                rw_mol.AddBond(c_idx, se_idx, Chem.BondType.DOUBLE)

                # Replace the single bond between N and C with a double bond
                rw_mol.RemoveBond(n_idx, c_idx)
                rw_mol.AddBond(n_idx, c_idx, Chem.BondType.DOUBLE)

                if hasattr(self,"bonds"):
                    # Also in the atoms object, using ligand indices
                    bond = self.atoms[c_idx].find_bond(se_idx, "ligand")
                    bond.order = 2.0
                    bond = self.atoms[n_idx].find_bond(c_idx, "ligand")
                    bond.order = 2.0

                # Get Current Charges
                fcharge_N  = rw_mol.GetAtomWithIdx(n_idx).GetFormalCharge()
                fcharge_Se = rw_mol.GetAtomWithIdx(se_idx).GetFormalCharge()

                # Update formal charges
                rw_mol.GetAtomWithIdx(n_idx).SetFormalCharge(fcharge_N+1)   
                rw_mol.GetAtomWithIdx(c_idx).SetFormalCharge(0)             
                rw_mol.GetAtomWithIdx(se_idx).SetFormalCharge(fcharge_Se-1) 

                # Also in ligand object
                self.atoms[n_idx].charge  = fcharge_N+1
                self.atoms[se_idx].charge = fcharge_Se-1

                # Update the molecule
                mol = rw_mol.GetMol()

        ########################
        # Fixes N2-CS+ pattern # 
        ########################
        pattern = Chem.MolFromSmarts("[N--]C#[S+]")

        # Check if the substructure exists in the molecule
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                atom_indices = {mol.GetAtomWithIdx(idx).GetAtomicNum(): idx for idx in match}
                n_idx = atom_indices[7]   # Atomic number of nitrogen
                c_idx = atom_indices[6]   # Atomic number of carbon
                s_idx = atom_indices[16]  # Atomic number of sulfur

                # Get the editable version of the molecule
                rw_mol = Chem.RWMol(mol)

                # Replace the triple bond between C and Se with a double bond
                rw_mol.RemoveBond(c_idx, s_idx)
                rw_mol.AddBond(c_idx, s_idx, Chem.BondType.DOUBLE)

                # Replace the single bond between N and C with a double bond
                rw_mol.RemoveBond(n_idx, c_idx)
                rw_mol.AddBond(n_idx, c_idx, Chem.BondType.DOUBLE)

                if hasattr(self,"bonds"):
                    # Also in the atoms object, using ligand indices
                    bond = self.atoms[c_idx].find_bond(s_idx, "ligand")
                    bond.order = 2.0
                    bond = self.atoms[n_idx].find_bond(c_idx, "ligand")
                    bond.order = 2.0

                # Get Current Charges
                fcharge_N  = rw_mol.GetAtomWithIdx(n_idx).GetFormalCharge()
                fcharge_S  = rw_mol.GetAtomWithIdx(s_idx).GetFormalCharge()

                # Update formal charges
                rw_mol.GetAtomWithIdx(n_idx).SetFormalCharge(fcharge_N+1)   
                rw_mol.GetAtomWithIdx(c_idx).SetFormalCharge(0)             
                rw_mol.GetAtomWithIdx(s_idx).SetFormalCharge(fcharge_S-1) 

                # Also in ligand object
                self.atoms[n_idx].charge  = fcharge_N+1
                self.atoms[s_idx].charge  = fcharge_S-1

                # Update the molecule
                mol = rw_mol.GetMol()

        #######################
        # Fixes N-CS+ pattern # 
        #######################
        pattern = Chem.MolFromSmarts("[N-]C#[S+]")

        # Check if the substructure exists in the molecule
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                atom_indices = {mol.GetAtomWithIdx(idx).GetAtomicNum(): idx for idx in match}
                n_idx = atom_indices[7]   # Atomic number of nitrogen
                c_idx = atom_indices[6]   # Atomic number of carbon
                s_idx = atom_indices[16]  # Atomic number of sulfur

                # Get the editable version of the molecule
                rw_mol = Chem.RWMol(mol)

                # Replace the triple bond between C and Se with a double bond
                rw_mol.RemoveBond(c_idx, s_idx)
                rw_mol.AddBond(c_idx, s_idx, Chem.BondType.DOUBLE)

                # Replace the single bond between N and C with a double bond
                rw_mol.RemoveBond(n_idx, c_idx)
                rw_mol.AddBond(n_idx, c_idx, Chem.BondType.DOUBLE)

                if hasattr(self,"bonds"):
                    # Also in the atoms object, using ligand indices
                    bond = self.atoms[c_idx].find_bond(s_idx, "ligand")
                    if bond is not None: bond.order = 2.0
                    bond = self.atoms[n_idx].find_bond(c_idx, "ligand")
                    if bond is not None: bond.order = 2.0

                # Get Current Charges
                fcharge_N  = rw_mol.GetAtomWithIdx(n_idx).GetFormalCharge()
                fcharge_S  = rw_mol.GetAtomWithIdx(s_idx).GetFormalCharge()

                # Update formal charges
                rw_mol.GetAtomWithIdx(n_idx).SetFormalCharge(fcharge_N+1)   
                rw_mol.GetAtomWithIdx(c_idx).SetFormalCharge(0)             
                rw_mol.GetAtomWithIdx(s_idx).SetFormalCharge(fcharge_S-1) 

                # Also in ligand object
                self.atoms[n_idx].charge  = fcharge_N+1
                self.atoms[s_idx].charge  = fcharge_S-1

                # Update the molecule
                mol = rw_mol.GetMol()

        #####################
        # Fixes Added atoms # 
        #####################
        natoms_rdkit = mol.GetNumAtoms() 
        if natoms_rdkit > self.natoms:
            if debug > 0: print(f"Fixing_Added_Atoms: {self.natoms=} {natoms_rdkit=}")
            rw_mol = Chem.RWMol(mol)
            to_remove = []
            for a in rw_mol.GetAtoms():
                i    = a.GetIdx()
                neig = a.GetNeighbors()
                if i > self.natoms-1:
                    label = a.GetSymbol()
                    if debug > 0: print(f"FIX_RDKIT: Atom {i=} {label=} flagged for removal")
                    to_remove.append(i)
                    for n in neig:
                        fcharge = n.GetFormalCharge()
                        if   label == 'H' : n.SetFormalCharge(fcharge-1) 
                        elif label == 'O' : n.SetFormalCharge(fcharge-2)
                        elif label == 'Cl': n.SetFormalCharge(fcharge+1)
                        else: 
                            if debug > 0: print(f"FIX_RDKIT: Unknown atom type as added atom {label=}")

            if len(to_remove) > 0:
                to_remove.sort(reverse=True)
                for idx, r in enumerate(to_remove):
                    if debug > 0: print(f"FIX_RDKIT: removing atom {r=}")
                    rw_mol.RemoveAtom(r)
            
            mol = rw_mol.GetMol()

        #######################################
        # Fixes Zwitterions in Adjacent Atoms #
        #######################################
        rw_mol = Chem.RWMol(mol)
        isnitrosyl = False
        for a in rw_mol.GetAtoms():
            ## Searches for adjacent atoms with opposite formal charges
            fcharge = a.GetFormalCharge()
            if fcharge > 0:
                label = a.GetSymbol()
                neig = a.GetNeighbors()
                neig_labels = [n.GetSymbol() for n in neig]
                idx = a.GetIdx()

                if debug > 0: print(f"FIX_RDKIT:{neig_labels=} {label=} {neig_labels.count('O')=}")
                fix = True
                ## Except Nitrosyls
                if label == 'N' and len(neig) == 3 and neig_labels.count('O') == 2: 
                    fix = False
                    isnitrosyl = True
                
                ## If Must be fixed
                if fix: 
                    for n in neig:
                        if n.GetIdx() > idx:
                            compensated_charge = n.GetFormalCharge() + fcharge 
                            ## Initiates Correction
                            if compensated_charge == 0:
                                a.SetFormalCharge(0)
                                n.SetFormalCharge(0)
        mol = rw_mol.GetMol()

        ##################
        ## Sanitize Step. Rdkit it too restrictive sometimes, but I still prefer to use it
        ##################
        ## If nitrosyl has been detected. Sanitize will fail
        if isnitrosyl:
            if debug > 0: print("Found nitrosyl group. RDKIT will fail at sanitizing. Preserving old one")
            if debug > 0: print("The mol-object attempt that failed is stored as self.failed_rdkit_obj")
            self.failed_rdkit_obj = mol
            return self.rdkit_obj

        ## Otherwise we try to sanitize
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:
            if debug > 0: print("Error Sanitizing the 'fixed' rdkit object. Preserving old one")
            if debug > 0: print("The mol-object attempt that failed is stored as self.failed_rdkit_obj")
            self.failed_rdkit_obj = mol
            return self.rdkit_obj

        ## If everything worked, we have new rdkit_obj
        self.rdkit_obj = mol
        return self.rdkit_obj

    ######
    def set_smiles_from_rdkit_obj(self, debug: int=0):
        from rdkit import Chem
        if not hasattr(self,"rdkit_obj"): return None
        self.smiles = Chem.MolToSmiles(self.rdkit_obj)
        return self.smiles

    ######
    # Moved to specie class
    ######
    #def set_bonds(self, debug: int=0):
    #    ## Creats bond objects using the information contained in the RDKit object
    #    ## The RDKit object is necessary, as it is the only one that contains the bond order information (lewis structure)
    #    if not hasattr(self,"rdkit_obj"): 
    #        if debug >= 1: print(f"LIG.SET_BONDS: Can't set bonds, ligand has not rdkit_object")
    #        return False
    #    natoms_rdkit = self.rdkit_obj.GetNumAtoms() 
    #    if debug >= 1: print(f"LIG.SET_BONDS: {self.formula=}, {self.subtype=} {self.smiles=}")

    #    if self.natoms == natoms_rdkit:  
    #        if debug >= 2: print(f"\tNumber of atoms in {self.subtype} object and RDKit object are equal: {self.natoms} {natoms_rdkit}")

    #        for idx, rdkit_atom in enumerate(self.rdkit_obj.GetAtoms()): # e.g. idx 0, 1, 2, 3, 4, 5, 6, 7, 8
    #            if debug >= 2: print(f"\t{idx=}", rdkit_atom.GetSymbol(), "Number of bonds in rdkit_obj:", len(rdkit_atom.GetBonds()))
    #            if len(rdkit_atom.GetBonds()) == 0:
    #                if debug >= 1: print(f"\tNO BONDS CREATED for {self.atoms[idx].label} due to no bonds in {self.subtype} RDKit object")
    #            else:
    #                for b in rdkit_atom.GetBonds():
    #                    bond_startatom = b.GetBeginAtomIdx()
    #                    bond_endatom   = b.GetEndAtomIdx()
    #                    bond_order     = b.GetBondTypeAsDouble()

    #                    ## Checks the labels involved in the bond
    #                    if self.atoms[bond_endatom].label != self.rdkit_obj.GetAtomWithIdx(bond_endatom).GetSymbol():
    #                        if debug >= 1: print(f"\tError with Bond EndAtom", self.atoms[bond_endatom].label, self.rdkit_obj.GetAtomWithIdx(bond_endatom).GetSymbol())
    #                    else:
    #                        if bond_endatom == idx:
    #                            start = bond_endatom
    #                            end   = bond_startatom
    #                        elif bond_startatom == idx:
    #                            start = bond_startatom
    #                            end   = bond_endatom      

    #                        # create new bond object and add it to both atoms
    #                        if end > start:
    #                            if debug >=2: print(f"\tBOND CREATED", idx, start, end, bond_order, self.atoms[start].label, self.atoms[end].label)
    #                            new_bond = bond(self.atoms[start], self.atoms[end], bond_order)
    #                            self.atoms[start].add_bond(new_bond)
    #                            self.atoms[end].add_bond(new_bond)
    #            
    #                if hasattr(self.atoms[idx], "bonds"):
    #                    if debug >=1 : print(f"\tBONDS", [(bd.atom1.label, bd.atom2.label, bd.order, round(bd.distance,3)) for bd in self.atoms[idx].bonds])
    #                else:
    #                    if self.natoms == 1:
    #                        if debug >=1: print(f"\tNO BONDS CREATED for {self.atoms[idx].label} because it is the only atom in {self.subtype} object")
    #                        pass
    #                    else:
    #                        if debug >=1: print(f"\tNO BONDS for {self.atoms[idx].label} with {self.subtype} RDKit object index {idx}. Please check the RDKit object.")
    #                        return False # return False if no bonds are created
    #    else:
    #        if debug >= 1: print(f"\tNumber of atoms in {self.subtype} object and RDKit object are different: {self.natoms} {natoms_rdkit}")
    #        if debug >= 2: print(f"\t{[(i, atom.label) for i, atom in enumerate(self.atoms)]}")
    #        if debug >= 2: print(f"\t{[(i, atom.GetSymbol()) for i, atom in enumerate(self.rdkit_obj.GetAtoms())]}")       
    #        non_bonded_atoms = list(range(0, natoms_rdkit))[self.natoms:]
    #        if debug >= 2: print(f"\tNON_BONDED_ATOMS", non_bonded_atoms)

    #        for idx, rdkit_atom in enumerate(self.rdkit_obj.GetAtoms()): # e.g. idx 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    #            if debug >= 2: print(f"\t{idx=}", rdkit_atom.GetSymbol(), "Number of bonds :", len(rdkit_atom.GetBonds()))
    #            if len(rdkit_atom.GetBonds()) == 0:
    #                if debug >= 1: print(f"\tNO BONDS CREATED for {rdkit_atom.GetSymbol()} due to no bonds in {self.subtype} RDKit object")
    #            else:
    #                for b in rdkit_atom.GetBonds():
    #                    bond_startatom = b.GetBeginAtomIdx()
    #                    bond_endatom   = b.GetEndAtomIdx()
    #                    bond_order     = b.GetBondTypeAsDouble()
  
    #                    if bond_startatom in non_bonded_atoms or bond_endatom in non_bonded_atoms:
    #                        if debug >= 2: print(f"\tNO BOND CREATED {bond_startatom=} or {bond_endatom=} is not in the self.atoms. It belongs to {non_bonded_atoms=}.")
    #                    else :
    #                        if bond_endatom == idx:
    #                            start = bond_endatom
    #                            end   = bond_startatom
    #                        elif bond_startatom == idx:
    #                            start = bond_startatom
    #                            end   = bond_endatom   

    #                        # create new bond object and add it to both atoms
    #                        if end > start:
    #                            if debug >=2: print(f"\tBOND CREATED", idx, start, end, bond_order, self.atoms[start].label, self.atoms[end].label)
    #                            new_bond = bond(self.atoms[start], self.atoms[end], bond_order)
    #                            self.atoms[start].add_bond(new_bond)
    #                            self.atoms[end].add_bond(new_bond)
    #            
    #                if idx not in non_bonded_atoms:
    #                    if hasattr(self.atoms[idx], "bonds"):
    #                        if debug >=2: 
    #                            print(f"\tBONDS", [(bd.atom1.label, bd.atom2.label, bd.order, round(bd.distance,3)) for bd in self.atoms[idx].bonds])
    #                    else:
    #                        if self.natoms == 1:
    #                            if debug >=1: print(f"\tNO BONDS CREATED for {self.atoms[idx].label} because it is the only atom in {self.subtype} object")
    #                            pass
    #                        else:
    #                            if debug >=1: print(f"\tNO BONDS for {self.atoms[idx].label} with {self.subtype} RDKit object index {idx}. Please check the RDKit object.")
    #                            return False # return False if no bonds are created
    #                else :
    #                    if debug >=1: print(f"\tNO BONDS for {rdkit_atom.GetSymbol()} with {self.subtype} RDKit object index {idx} because it is an added atom")
    #    self.has_bonds = True
    #    return True 

###############
#### GROUP ####
###############
class group(specie):
    def __init__(self, labels: list, coord: list, frac_coord: list=None, radii: list=None) -> None:
        if frac_coord is not None:        
            assert len(frac_coord) == len(coord)
            self.frac_coord = frac_coord
        specie.__init__(self, labels, coord, frac_coord, radii)
        self.subtype = "group"

    ######
    def __repr__(self, indirect: bool=False) -> str:
        to_print                   = ''
        if not indirect: to_print += '------------------------------------------------\n'
        if not indirect: to_print += '------------- SCOPE GROUP Object ---------------\n'
        if not indirect: to_print += '------------------------------------------------\n'
        to_print += specie.__repr__(self, indirect=True)
        if not indirect: to_print += '\n'
        return to_print

    ######
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

    ######
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

    ######
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

    ######
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

    ######
    def get_denticity(self, debug: int=0):
        self.denticity = 0
        for a in self.atoms: 
            if debug > 0: print(f"GROUP.GET_DENTICITY. Evaluating Atom with {a.madjnum=} and so far {self.denticity=}")
            self.denticity += a.madjnum      
        return self.denticity

###############
### IMPORTS ###
###############
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

######
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

######
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
    else:
        if debug > 0: print(f"IMPORT LIGAND: SMILES could not be imported")

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
            
    return new_ligand

######
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
