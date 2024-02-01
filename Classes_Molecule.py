import numpy as np
from Scope.Adapted_from_cell2mol import * 
from Scope.Other import get_metal_idxs

################################
####  BASIS FOR CELL2MOL 2  ####
################################
class specie(object):
    def __init__(self, labels: list, coord: list, indices: list=None, radii: list=None, parent: object=None) -> None:

       # Sanity Checks
        assert len(labels) == len(coord)
        if indices is not None: assert len(labels) == len(indices)
        if radii   is not None: assert len(labels) == len(radii)

        # Optional Information
        if radii   is not None: self.radii   = radii
        else:                   self.radii   = get_radii(labels)
        if parent  is not None: self._parent = parent
#        if parent  is not None: self.occurence = parent.get_occurrence(self)

        self.type      = "specie"
        self.version   = "0.1"
        self.labels    = labels
        self.coord     = coord
        self.formula   = labels2formula(labels)
        self.eleccount = labels2electrons(labels)
        self.natoms    = len(labels)
        self.iscomplex = any((elemdatabase.elementblock[l] == "d") or (elemdatabase.elementblock[l] == "f") for l in self.labels)

        if indices is not None: self.indices = indices
        else:                   self.indices = [*range(0,self.natoms,1)]

        self.cov_factor   = 1.3
        self.metal_factor = 1.0

    def get_centroid(self):
        self.centroid = get_centroid(self.coord)
        return self.centroid

    def set_fractional_coord(self, frac_coord: list) -> None:
        self.frac_coord = frac_coord 

    def get_atomic_numbers(self):
        if not hasattr(self,"atoms"): self.set_atoms()
        self.atnums = []
        for at in self.atoms:
            self.atnums.append(at.atnum)
        return self.atnums

    def set_atoms(self, atomlist: object=None):
        if atomlist is None:
            self.atoms = []
            for idx, l in enumerate(self.labels):
                newatom = atom(l, self.coord[idx], parent=self, index=idx, radii=self.radii[idx])
                self.atoms.append(newatom)
        else:
            self.atoms = atomlist

    def set_element_count(self, heavy_only: bool=False):
        self.element_count = get_element_count(self.labels, heavy_only=heavy_only)
        return self.element_count

    def set_adj_types(self):
        if not hasattr(self,"adjmat"): self.get_adjmatrix()
        self.adj_types = get_adjacency_types(self.labels, self.adjmat)
        return self.adj_types

    def set_adjacency_parameters(self, cov_factor: float, metal_factor: float) -> None:
        # Stores the covalentradii factor and metal factor that were used to generate the molecule
        self.cov_factor   = cov_factor
        self.metal_factor = metal_factor

    def set_charges(self, totcharge: int, atomic_charges: list=None) -> None:
        self.totcharge = totcharge
        if atomic_charges is not None:
            self.atomic_charges = atomic_charges
            if not hasattr(self,"atoms"): self.set_atoms()
            for idx, a in enumerate(self.atoms):
                a.set_charge(self.atomic_charges[idx])

    def get_adjmatrix(self):
        isgood, adjmat, adjnum = get_adjmatrix(self.labels, self.coord, self.cov_factor, self.radii)
        if isgood:
            self.adjmat = adjmat
            self.adjnum = adjnum
        else:
            self.adjmat = None
            self.adjnum = None
        return self.adjmat, self.adjnum

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

    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> SPECIE >>>                                  \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        if hasattr(self,'subtype'): to_print += f' Sub-Type              = {self.subtype}\n'
        to_print += f' Number of Atoms       = {self.natoms}\n'
        to_print += f' Formula               = {self.formula}\n'
        if hasattr(self,"adjmat"):     to_print += f' Has Adjacency Matrix  = YES\n'
        else:                          to_print += f' Has Adjacency Matrix  = NO \n'
        if hasattr(self,"occurrence"): to_print += f' Occurrence in Parent  = {self.occurrence}\n'
        if hasattr(self,"totcharge"):  to_print += f' Total Charge          = {self.totcharge}\n'
        if hasattr(self,"spin"):       to_print += f' Spin                  = {self.spin}\n'
        if hasattr(self,"smiles"):     to_print += f' SMILES                = {self.smiles}\n'
        to_print += '---------------------------------------------------\n'
        return to_print

###############
### MOLECULE ##
###############
class molecule(specie):
    def __init__(self, labels: list, coord: list, indices: list=None, radii: list=None, parent: object=None) -> None:
        self.subtype = "molecule"
        specie.__init__(self, labels, coord, indices, radii, parent)

############
    def split_complex(self, debug: int=0):
        if not self.iscomplex: 
            self.ligands = None
            self.metals  = None
            #print("MOLECULE.SPLIT_COMPLEX: This molecule is not a transition metal complex");
    
        else: 
            self.ligands = []
            self.metals  = []
            # Identify Metals and the rest
            metal_idx = get_metal_idxs(self.labels, debug=debug)
            rest_idx  = list(idx for idx in self.indices if idx not in metal_idx)
            # Split the "rest" to obtain the ligands
            rest_labels  = extract_from_list(rest_idx, self.labels, dimension=1)
            rest_coords  = extract_from_list(rest_idx, self.coord, dimension=1)
            rest_indices = extract_from_list(rest_idx, self.indices, dimension=1)
            if hasattr(self,"cov_factor"): blocklist = split_species(rest_labels, rest_coords, factor=self.cov_factor)
            else:                          blocklist = split_species(rest_labels, rest_coords)
            ## Arranges Ligands
            for b in blocklist:
                lig_labels  = extract_from_list(b, self.labels, dimension=1) 
                lig_coord   = extract_from_list(b, self.coord, dimension=1) 
                lig_indices = extract_from_list(b, self.indices, dimension=1) 
                lig_radii   = extract_from_list(b, self.radii, dimension=1) 
                newligand   = ligand(lig_labels, lig_coord, indices=lig_indices, radii=lig_radii, parent=self)
                # If fractional coordinates are available...
                if hasattr(self,"frac_coord"): 
                    lig_frac_coords  = extract_from_list(b, self.frac_coord, dimension=1)
                    newligand.set_fractional_coord(lig_frac_coords)
                self.ligands.append(newligand)
            ## Arranges Metals
            for m in metal_idx:
                newmetal    = metal(self.labels[m], self.coord[m], self.indices[m], self.radii[m], parent=self)
                self.metals.append(newmetal)
        return self.ligands, self.metals
        
    def get_metal_adjmatrix(self):
        if not hasattr(self,"adjmat"): self.get_adjmatrix() 
        if self.adjmat is None: return None, None
        madjmat = np.zeros((self.natoms,self.natoms))
        madjnum = np.zeros((self.natoms)) 
        metal_idx = get_connected_idx(self, debug=debug)
        for idx, row in enumerate(self.adjmat):
            for jdx, col in enumerate(row):
                if not idx in metal_idx: madjmat[idx,jdx] = int(0) 
                else:                    madjmat[idx,jdx] = adjmat[idx,jdx]
        for i in range(0, len(adjmat)):
            madjnum[i] = np.sum(madjmat[i, :])
        return self.madjmat, self.madjnum

############

###############
### LIGAND ####
###############
class ligand(specie):
    def __init__(self, labels: list, coord: list, indices: list=None, radii: list=None, parent: object=None) -> None:
        self.subtype  = "ligand"
        specie.__init__(self, labels, coord, indices, radii, parent)

    def set_hapticity(self, hapttype):
        self.hapticity = True 
        self.hapttype  = hapttype 

    def get_connected_idx(self, debug: int=0):
       self.connected_idx = [] 
       if not hasattr(self._parent,"adjmat"): madjmat, madjnum = self._parent.get_metal_adjmatrix()
       self.madjmat = extract_from_list(self.indices, madjmat, dimension=2)
       self.madjnum = extract_from_list(self.indices, madjnum, dimension=1)
       for idx, con in enumerate(self.madjnum):
           if con > 0: self.connected_idx.append(idx) 
       return self.connected_idx 
        
    def split_ligand(self, debug: int=0):
        self.groups = []
        # Identify Connected and Unconnected atoms (to the metal)
        if not hasattr(self,"connected_idx"): self.get_connected_idx()
        conn_idx = self.connected_idx
        rest_idx = list(idx for idx in self.indices if idx not in conn_idx)

        # Split the "ligand to obtain the groups
        conn_labels  = extract_from_list(conn_idx, self.labels, dimension=1)
        conn_coords  = extract_from_list(conn_idx, self.coord, dimension=1)
        conn_indices = extract_from_list(conn_idx, self.indices, dimension=1)

        if hasattr(self,"cov_factor"): blocklist = split_species(conn_labels, conn_coords, factor=self.cov_factor)
        else:                          blocklist = split_species(conn_labels, conn_coords)
        ## Arranges Groups 
        for b in blocklist:
            gr_labels  = extract_from_list(b, self.labels, dimension=1)
            gr_coord   = extract_from_list(b, self.coord, dimension=1)
            gr_indices = extract_from_list(b, self.indices, dimension=1)
            gr_radii   = extract_from_list(b, self.radii, dimension=1)
            newgroup   = group(gr_labels, gr_coord, gr_indices, gr_radii, parent=self)
            self.groups.append(newgroup)
            
        return self.groups

###############
#### GROUP ####
###############
class group(specie):
    def __init__(self, labels: list, coord: list, indices: list=None, radii: list=None, parent: object=None) -> None:
        self.subtype = "group"
        specie.__init__(self, labels, coord, indices, radii, parent)

###############
### ATOM ######
###############
class atom(object):
    def __init__(self, label: str, coord: list, index: int=None, radii: float=None, frac_coord: list=None, parent: object=None) -> None:
        self.type    = "atom"
        self.version = "0.1"
        self.label = label
        self.coord = coord
        self.atnum = elemdatabase.elementnr[label]
        self.block = elemdatabase.elementblock[label]

        if index is not None:      self.index = index
        if radii is None:          self.radii = getradii(label)
        else:                      self.radii = radii
        if frac_coord is not None: self.frac_coord = frac_coord
        if parent is not None:     self._parent  = parent
        if parent is not None:     self.occurence = parent.get_occurrence(self)

    def set_adjacency_parameters(self, cov_factor: float, metal_factor: float) -> None:
        self.cov_factor   = cov_factor
        self.metal_factor = metal_factor

    def set_charge(self, charge: int) -> None:
        self.charge = charge 

    def information(self, cov_factor: float, metal_factor: float) -> None:
        self.cov_factor   = cov_factor
        self.metal_factor = metal_factor

    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> ATOM >>>                                    \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        if hasattr(self,'subtype'): to_print += f' Sub-Type              = {self.subtype}\n'
        to_print += f' Label                 = {self.label}\n'
        to_print += f' Atomic Number         = {self.atnum}\n'
        if hasattr(self,"occurrence"): to_print += f' Occurrence in Parent  = {self.occurrence}\n'
        if hasattr(self,"charge"):     to_print += f' Atom Charge           = {self.charge}\n'
        to_print += '----------------------------------------------------\n'
        return to_print

###############
#### METAL ####
###############
class metal(atom):
    def __init__(self, label: str, coord: list, index: int=None, radii: float=None, frac_coord: list=None, parent: object=None) -> None:
        self.subtype      = "metal"
        atom.__init__(self, label, coord, index, radii, frac_coord, parent)

##############
#### CELL ####
##############
class cell(object):
    def __init__(self, refcode: str, labels: list, pos: list, cellvec: list, cellparam: list) -> None:
        self.version    = "0.1"
        self.type       = "cell"
        self.refcode    = refcode
        self.labels     = labels 
        self.coord      = pos
        self.pos        = pos
        self.cellvec    = cellvec
        self.cellparam  = cellparam
        self.natoms     = len(labels)
 
    def set_fractional_coord(self, frac_coord: list) -> None:
        self.frac_coord = frac_coord 

    def set_moleclist(self, moleclist: list) -> None:
        self.moleclist = moleclist

    def get_moleclist(self):
        if not hasattr(self,"labels") or not hasattr(self,"pos"): return None
        if len(self.labels) == 0 or len(self.pos) == 0: return None
        factor = 1.3

        #try:
        blocklist = split_species(self.labels, self.pos, factor=factor)
        self.moleclist = []
        for b in blocklist:
            mol_labels  = extract_from_list(b, self.labels, dimension=1)
            mol_coords  = extract_from_list(b, self.coord, dimension=1)
            newmolec    = molecule(mol_labels, mol_coords)
            # If fractional coordinates are available...
            if hasattr(self,"frac_coord"): 
                mol_frac_coords  = extract_from_list(b, self.frac_coord, dimension=1)
                newmolec.set_fractional_coord(mol_frac_coords)
            # This must be below the frac_coord, so they are carried on to the ligands
            if newmolec.iscomplex: 
                newmolec.split_complex()
            self.moleclist.append(newmolec)
        return self.moleclist
   
    def arrange_cell_coord(self): 
        ## Updates the cell coordinates preserving the original atom ordering
        ## Do do so, it uses the variable atlist stored in each molecule
        self.coord = np.zeros((self.natoms,3))
        for mol in self.moleclist:
            for z in zip(mol.indices, mol.coord):
                for i in range(0,3):
                    self.coord[z[0]][i] = z[1][i]
        self.coord = np.ndarray.tolist(self.coord)

    def get_occurrence(self, substructure: object) -> int:
        occurrence = 0
        ## Molecules in Cell
        if hasattr(substructure,"subtype") and hasattr(self,"moleclist"): 
            if substructure.subtype == 'molecule':
                for m in self.moleclist:
                    issame = compare_species(substructure, m)
                    if issame: occurrence += 1
        return occurrence

    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> CELL >>>                                    \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        to_print += f' Refcode               = {self.refcode}\n'
        to_print += f' Num Atoms             = {self.natoms}\n'
        to_print += f' Cell Parameters a:c   = {self.cellparam[0:3]}\n'
        to_print += f' Cell Parameters al:ga = {self.cellparam[3:6]}\n'
        if hasattr(self,"moleclist"):  to_print += f' List of Molecules     = \n{self.moleclist}\n'
        to_print += '---------------------------------------------------\n'
        return to_print

######################
####    IMPORT    ####
######################
def import_cell(cell: object) -> object:
    assert hasattr(cell,"labels") and (hasattr(cell,"coord") or hasattr(cell,"pos"))
    assert hasattr(cell,"cellvec")
    assert hasattr(cell,"refcode")

    labels     = cell.labels
    refcode    = cell.refcode
    if   hasattr(cell,"coord"):       coord      = cell.coord
    elif hasattr(cell,"pos"):         coord      = cell.pos

    cellvec    = cell.cellvec
    if   hasattr(cell,"cellparam"):   cellparam  = cell.cellparam
    else:                             cellparam  = cellvec_2_cellparam(cell.cellvec)

    if   hasattr(cell,"frac_coord"):  frac_coord = cell.frac_coord
    else:                             frac_coord = cart2frac(coord, cell.cellvec)

    new_cell = Scope.Classes_Molecule.cell(refcode, labels, coord, cellvec, cellparam)
    new_cell.set_fractional_coord(frac_coord)

    if not hasattr(cell,"moleclist"): new_cell.get_moleclist()
    else: 
        moleclist = []
        for mol in cell.moleclist: 
            new_mol = import_molecule(mol, parent=new_cell)
            moleclist.append(new_mol)
        new_cell.set_moleclist(moleclist)

    return new_cell

################################
def import_molecule(mol: object, parent: object=None) -> object:
    assert hasattr(mol,"labels") and (hasattr(mol,"coord") or hasattr(mol,"pos"))
    labels     = mol.labels
    if   hasattr(mol,"coord"):       coord      = mol.coord
    elif hasattr(mol,"pos"):         coord      = mol.pos

    if   hasattr(mol,"indices"):     indices    = mol.indices
    elif hasattr(mol,"atlist"):      indices    = mol.atlist
    else:                            indices    = None

    if   hasattr(mol,"radii"):       radii      = mol.radii
    else:                            radii      = None          

    new_molec = Scope.Classes_Molecule.molecule(labels, coord, indices, radii, parent)

    ## Charges
    if hasattr(mol,"totcharge") and hasattr(mol,"atcharge"):
        new_molec.set_charges(mol.totcharge, mol.atcharge)
    elif hasattr(mol,"totcharge") and not hasattr(mol,"atcharge"):
        new_molec.set_charges(mol.totcharge)

    ## Smiles
    if   hasattr(mol,"Smiles"): new_molec.smiles = mol.Smiles
    elif hasattr(mol,"smiles"): new_molec.smiles = mol.smiles        

    ## Substructures
    if not hasattr(mol,"ligandlist") or not hasattr(mol,"metalist"):  new_molec.split_complex()
    if hasattr(mol,"ligandlist"): 
        ligands = []
        for lig in mol.ligandlist: 
            new_lig = import_ligand(lig, parent=new_molec)
            ligands.append(new_lig)
        new_molec.ligands = ligands
    if hasattr(mol,"metalist"): 
        metals = []
        for met in mol.metalist: 
            new_atom = import_atom(met, parent=new_molec)
            metals.append(new_atom)
        new_molec.metals = metals

    ## Atoms
    if not hasattr(mol,"atoms"):  new_molec.set_atoms()
    else: 
        atoms = []
        for at in mol.atoms: 
            new_atom = import_atom(at, parent=new_molec)
            atoms.append(new_atom)
        new_molec.set_atoms(atoms)

    return new_molec

################################
def import_ligand(lig: object, parent: object=None) -> object:
    assert hasattr(lig,"labels") and (hasattr(lig,"coord") or hasattr(lig,"pos"))
    labels     = lig.labels
    if   hasattr(lig,"coord"):       coord      = lig.coord
    elif hasattr(lig,"pos"):         coord      = lig.pos

    if   hasattr(lig,"indices"):     indices    = lig.indices
    elif hasattr(lig,"atlist"):      indices    = lig.atlist
    else:                            indices    = None

    if   hasattr(lig,"radii"):       radii      = lig.radii
    else:                            radii      = None          

    new_ligand = Scope.Classes_Molecule.ligand(labels, coord, indices, radii, parent)
    
    ## Charges
    if hasattr(lig,"totcharge") and hasattr(lig,"atcharge"):
        new_ligand.set_charges(lig.totcharge, lig.atcharge)
    elif hasattr(lig,"totcharge") and not hasattr(lig,"atcharge"):
        new_ligand.set_charges(lig.totcharge)

    ## Smiles
    if   hasattr(lig,"Smiles"): new_ligand.smiles = lig.Smiles
    elif hasattr(lig,"smiles"): new_ligand.smiles = lig.smiles     
    
    ## Substructures
    if not hasattr(lig,"grouplist"):  new_ligand.split_ligand()
    else: 
        groups = []
        for gr in lig.grouplist: 
            new_group = import_group(gr, parent=new_ligand)
            groups.append(new_group)
    new_ligand.groups = groups
    
    ## Atoms
    if not hasattr(lig,"atoms"):  new_ligand.set_atoms()
    else: 
        atoms = []
        for at in lig.atoms: 
            new_atom = import_atom(at, parent=new_ligand)
            atoms.append(new_atom)
        new_ligand.set_atoms(atoms)
            
    return new_ligand

################################
def import_group(group: object, parent: object=None) -> object:
    assert hasattr(group,"atlist") or hasattr(group,"indices")
    
    if   hasattr(group,"labels"):      labels     = group.labels
    elif parent is not None:
        if hasattr(parent,"labels"):   labels     = extract_from_list(group.atlist, parent.labels, dimension=1)
        
    if   hasattr(group,"coord"):       coord      = group.coord
    elif hasattr(group,"pos"):         coord      = group.pos
    elif parent is not None:
        if hasattr(parent,"coord"):    coord     = extract_from_list(group.atlist, parent.coord, dimension=1)
        elif hasattr(parent,"coord"):  coord     = extract_from_list(group.atlist, parent.pos, dimension=1)

    if   hasattr(group,"indices"):     indices    = group.indices
    elif hasattr(group,"atlist"):      indices    = group.atlist
    else:                              indices    = None

    if   hasattr(group,"radii"):       radii      = group.radii
    else:                              radii      = None          

    new_group = Scope.Classes_Molecule.group(labels, coord, indices, radii, parent)
    
    ## Charges
    if hasattr(group,"totcharge") and hasattr(group,"atcharge"):
        new_group.set_charges(group.totcharge, group.atcharge)
    elif hasattr(group,"totcharge") and not hasattr(group,"atcharge"):
        new_group.set_charges(group.totcharge)

    ## Smiles
    if   hasattr(group,"Smiles"): new_group.smiles = group.Smiles
    elif hasattr(group,"smiles"): new_group.smiles = group.smiles     
    
    ## Atoms
    if not hasattr(group,"atoms"):  new_group.set_atoms()
    else: 
        atoms = []
        for at in group.atoms: 
            new_atom = import_atom(at, parent=new_group)
            atoms.append(new_atom)
        new_group.set_atoms(atoms)
            
    return new_group

################################
def import_atom(atom: object, parent: object=None) -> object:
    assert hasattr(atom,"label") and (hasattr(atom,"coord") or hasattr(atom,"pos"))
    label     = atom.label
    
    if   hasattr(atom,"coord"):      coord      = atom.coord
    elif hasattr(atom,"pos"):        coord      = atom.pos

    if   hasattr(atom,"index"):      indices    = atom.index
    else:                            indices    = None

    if   hasattr(atom,"radii"):      radii      = atom.radii
    else:                            radii      = None          

    if   hasattr(atom,"block"):      block      = atom.block
    else:                            block      = elemdatabase.elementblock[label]    
    
    if block == 'd' or block == 'f': new_atom = Scope.Classes_Molecule.metal(label, coord, indices, radii, parent)
    else:                            new_atom = Scope.Classes_Molecule.atom(label, coord, indices, radii, parent)
    
    ## Charge
    if hasattr(atom,"charge"):       new_atom.set_charge(atom.charge)
            
    return new_atom
