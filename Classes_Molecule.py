import numpy as np
from Scope.Adapted_from_cell2mol import * 
from Scope.Other import get_metal_idxs
from Scope.Unit_cell_tools import * 
#from Scope.Reconstruct import * #error when loading here 

################################
####  BASIS FOR CELL2MOL 2  ####
################################
class specie(object):
    def __init__(self, labels: list, coord: list, parent_indices: list=None, radii: list=None, parent: object=None) -> None:

       # Sanity Checks
        assert len(labels) == len(coord)
        if parent_indices is not None: assert len(labels) == len(parent_indices)
        if radii   is not None: assert len(labels) == len(radii)

        # Optional Information
        if radii   is not None: self.radii   = radii
        else:                   self.radii   = get_radii(labels)
        if parent  is not None: self.parent = parent

        self.type      = "specie"
        self.version   = "0.1"
        self.labels    = labels
        self.coord     = coord
        self.formula   = labels2formula(labels)
        self.eleccount = labels2electrons(labels)
        self.natoms    = len(labels)
        self.iscomplex = any((elemdatabase.elementblock[l] == "d") or (elemdatabase.elementblock[l] == "f") for l in self.labels)

        if parent_indices is not None: self.parent_indices = parent_indices
        else:                          self.parent_indices = [*range(0,self.natoms,1)]
        self.indices = [*range(0,self.natoms,1)]

        self.cov_factor   = 1.3
        self.metal_factor = 1.0

    def get_centroid(self):
        self.centroid = compute_centroid(self.coord)
        if hasattr(self,"frac_coord"): self.frac_centroid = compute_centroid(self.frac_coord)
        return self.centroid

    def set_fractional_coord(self, frac_coord: list) -> None:
        self.frac_coord = frac_coord 

    def get_fractional_coord(self, cell_vector=None) -> None:
        from Scope.Reconstruct import cart2frac
        if cell_vector is None:
            if hasattr(self,"parent"):
                if hasattr(self.parent,"cellvec"): cell_vector = self.parent.cellvec.copy()
            else:     print("CLASS_SPECIE: get_fractional coordinates. Missing cell vector. Please provide it"); return None
        else:         
            self.frac_coord = cart2frac(self.coord, cell_vector)
        return self.frac_coord

    def get_atomic_numbers(self):
        if not hasattr(self,"atoms"): self.set_atoms()
        self.atnums = []
        for at in self.atoms:
            self.atnums.append(at.atnum)
        return self.atnums

    ############
    def set_atoms(self, atomlist=None, overwrite_parent: bool=False, create_adjacencies: bool=False, debug: int=0):
        ## If the atom objects already exist, and you want to set them in self from a different specie
        if atomlist is not None:
            if debug > 0: print(f"SPECIE.SET_ATOMS: received {atomlist=}")
            self.atoms = atomlist.copy()
            if overwrite_parent:
                for at in self.atoms:
                    setattr(at,"parent",self) #at.parent = self

        ## If not, that is, if the atom objects must be created from scratch....
        else:
            self.atoms = []
            for idx, l in enumerate(self.labels):
                if debug > 0: print(f"SPECIE.SET_ATOMS: creating atom for label {l}")
                ## For each l in labels, create an atom class object.
                ismetal = elemdatabase.elementblock[l] == "d" or elemdatabase.elementblock[l] == "f"
                if debug > 0: print(f"SPECIE.SET_ATOMS: {ismetal=}")
                if ismetal: newatom = metal(l, self.coord[idx], parent_index=idx, radii=self.radii[idx], parent=self)
                else:       newatom = atom(l, self.coord[idx], parent_index=idx, radii=self.radii[idx], parent=self)
                if debug > 0: print(f"SPECIE.SET_ATOMS: added atom to specie: \n{newatom}")
                self.atoms.append(newatom)

        if create_adjacencies:
            if not hasattr(self,"adjmat"):  self.get_adjmatrix()
            if not hasattr(self,"madjmat"): self.get_metal_adjmatrix()
            if self.adjmat is not None and self.madjmat is not None:
                for idx, at in enumerate(self.atoms):
                    at.set_adjacencies(self.adjmat[idx],self.madjmat[idx],self.adjnum[idx],self.madjnum[idx])

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

    def set_charges(self, totcharge: int=None, atomic_charges: list=None) -> None:
        ## Sets total charge  
        if totcharge is not None:                              self.totcharge = totcharge
        elif totcharge is None and atomic_charges is not None: self.totcharge = np.sum(atomic_charges)
        elif totcharge is None and atomic_charges is None:     self.totcharge = "Unknown" 
        ## Sets atomic charges
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

    def get_metal_adjmatrix(self):
        isgood, madjmat, madjnum = get_adjmatrix(self.labels, self.coord, self.cov_factor, self.radii, metal_only=True)
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

    def view(self):
        import plotly.graph_objects as go
        from Scope.Read_Write import set_scene
        from Scope.Elementdata import ElementData  
        elemdatabase = ElementData()

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

        set_scene(fig, np.array(self.coord))
        fig.show()


    ## To be implemented
    def __add__(self, other):
        if not isinstance(other, type(self)): return self
        return self

    ############
    def __repr__(self, indirect: bool=False):
        to_print = ""
        if not indirect: to_print  += f'------------- Cell2mol SPECIE Object --------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        if hasattr(self,'subtype'): to_print += f' Sub-Type               = {self.subtype}\n'
        to_print += f' Number of Atoms       = {self.natoms}\n'
        to_print += f' Formula               = {self.formula}\n'
        if hasattr(self,"adjmat"):     to_print += f' Has Adjacency Matrix  = YES\n'
        else:                          to_print += f' Has Adjacency Matrix  = NO \n'
        if hasattr(self,"parent"):
            if self.parent is not None:    to_print += f' Has Parent            = YES\n'
            else:                          to_print += f' Has Parent            = NO \n'
        if hasattr(self,"occurrence"): to_print += f' Occurrence in Parent  = {self.occurrence}\n'
        if hasattr(self,"totcharge"):  to_print += f' Total Charge          = {self.totcharge}\n'
        if hasattr(self,"spin"):       to_print += f' Spin                  = {self.spin}\n'
        if hasattr(self,"smiles"):     to_print += f' SMILES                = {self.smiles}\n'
        if not indirect: to_print += '---------------------------------------------------\n'
        return to_print

###############
### MOLECULE ##
###############
class molecule(specie):
    def __init__(self, labels: list, coord: list, parent_indices: list=None, radii: list=None, parent: object=None) -> None:
        self.subtype = "molecule"
        specie.__init__(self, labels, coord, parent_indices, radii, parent)

    def split_complex(self, debug: int=0):
         if not hasattr(self,"atoms"): self.set_atoms()
         if not self.iscomplex:        self.ligands = None; self.metals = None
         else:
             self.ligands = []
             self.metals  = []
             # Identify Metals and the rest
             metal_idx = list([self.indices[idx] for idx in get_metal_idxs(self.labels, debug=debug)])
             rest_idx  = list(idx for idx in self.indices if idx not in metal_idx)
             if debug > 0 :  print(f"MOLECULE.SPLIT COMPLEX: labels={self.labels}")
             if debug > 0 :  print(f"MOLECULE.SPLIT COMPLEX: parent_indices={self.parent_indices}")
             if debug > 0 :  print(f"MOLECULE.SPLIT COMPLEX: metal_idx={metal_idx}")
             if debug > 0 :  print(f"MOLECULE.SPLIT COMPLEX: rest_idx={rest_idx}")

             # Split the "rest" to obtain the ligands
             rest_labels  = extract_from_list(rest_idx, self.labels, dimension=1)
             rest_coord   = extract_from_list(rest_idx, self.coord, dimension=1)
             rest_indices = extract_from_list(rest_idx, self.indices, dimension=1)
             rest_radii   = extract_from_list(rest_idx, self.radii, dimension=1)
             rest_atoms   = extract_from_list(rest_idx, self.atoms, dimension=1)
             if debug > 0:
                 print(f"SPLIT COMPLEX: rest labels: {rest_labels}")
                 print(f"SPLIT COMPLEX: rest coord: {rest_coord}")
                 print(f"SPLIT COMPLEX: rest indices: {rest_indices}")
                 print(f"SPLIT COMPLEX: rest radii: {rest_radii}")

             if hasattr(self,"frac_coord"): rest_frac = extract_from_list(rest_idx, self.frac_coord, dimension=1)
             if debug > 0: print(f"SPLIT COMPLEX: rest labels: {rest_labels}")
             if debug > 0: print(f"SPLIT COMPLEX: splitting species with {len(rest_labels)} atoms in block")
             if hasattr(self,"cov_factor"): blocklist = split_species(rest_labels, rest_coord, radii=rest_radii, cov_factor=self.cov_factor, debug=debug)
             else:                          blocklist = split_species(rest_labels, rest_coord, radii=rest_radii, cov_factor=self.cov_factor, debug=debug)

             ## Arranges Ligands
             for b in blocklist:
                 if debug > 0: print(f"PREPARING BLOCK: {b}")
                 lig_indices = extract_from_list(b, rest_indices, dimension=1)
                 lig_labels  = extract_from_list(b, rest_labels, dimension=1)
                 lig_coord   = extract_from_list(b, rest_coord, dimension=1)
                 lig_radii   = extract_from_list(b, rest_radii, dimension=1)
                 lig_atoms   = extract_from_list(b, rest_atoms, dimension=1)
                 if debug > 0: print(f"CREATING LIGAND: {labels2formula(lig_labels)}")
                 if debug > 0: print(f"CREATING LIGAND with atoms: {lig_atoms}")
                 newligand   = ligand(lig_labels, lig_coord, parent_indices=lig_indices, radii=lig_radii, parent=self)
                 # We pass the molecule atoms to the ligand, and create an intermediate parent (ligand)
                 newligand.set_atoms(atomlist=lig_atoms, overwrite_parent=True)
                 # If fractional coordinates are available...
                 if hasattr(self,"frac_coord"):
                     lig_frac_coord = extract_from_list(b, rest_frac, dimension=1)
                     newligand.set_fractional_coord(lig_frac_coord)
                 self.ligands.append(newligand)

             ## Arranges Metals
             for m in metal_idx:
                 newmetal    = metal(self.labels[m], self.coord[m], self.indices[m], self.radii[m], parent=self)
                 self.metals.append(newmetal)
         return self.ligands, self.metals
        
    def __repr__(self):
        to_print = ""
        to_print += f'------------- Cell2mol MOLECULE Object --------------\n'
        to_print += specie.__repr__(self, indirect=True)
        if hasattr(self,"ligands"):  
            if self.ligands is not None: to_print += f' # Ligands             = {len(self.ligands)}\n'
        if hasattr(self,"metals"):   
            if self.metals is not None:  to_print += f' # Metals              = {len(self.metals)}\n'
        to_print += '---------------------------------------------------\n'
        return to_print

############

###############
### LIGAND ####
###############
class ligand(specie):
    def __init__(self, labels: list, coord: list, parent_indices: list=None, radii: list=None, parent: object=None) -> None:
        self.subtype  = "ligand"
        specie.__init__(self, labels, coord, parent_indices, radii, parent)

    def set_hapticity(self, hapttype):
        self.hapticity = True 
        self.hapttype  = hapttype 

    #######################################################
    def get_connected_idx(self, debug: int=0):
        self.connected_idx = []
        if hasattr(self,"madjnum"):
            for idx, con in enumerate(self.madjnum):
                if con > 0: self.connected_idx.append(self.indices[idx])
        else:
            ## Remember madjmat should not be computed at the ligand level. Since the metal is not there.
            ## This is why we compute it at the parent (molecule) level
            if not hasattr(self.parent,"madjnum"): self.parent.get_metal_adjmatrix()
            self.madjmat = extract_from_list(self.parent_indices, self.parent.madjmat, dimension=2)    ## Here we use parent_indices because we're operating with a molecule variable (madjmat and madjnum)
            self.madjnum = extract_from_list(self.parent_indices, self.parent.madjnum, dimension=1)
            for idx, con in enumerate(self.madjnum):
                if con > 0: self.connected_idx.append(self.indices[idx])
        return self.connected_idx

    #######################################################
    def get_connected_atoms(self, debug: int=0):
        if not hasattr(self,"atoms"):         self.set_atoms()
        if not hasattr(self,"connected_idx"): self.get_connected_idx()
        self.connected_atoms = []
        for idx, at in enumerate(self.atoms):
            if idx in self.connected_idx: # and at.mconnec > 0:
                self.connected_atoms.append(at)
            elif idx in self.connected_idx: # and at.mconnec == 0:
                print("WARNING: Atom appears in connected_idx, but has mconnec=0")
        return self.connected_atoms
        
    #######################################################
    def split_ligand(self, debug: int=0):
        # Split the "ligand to obtain the groups
        self.groups = []
        # Identify Connected and Unconnected atoms (to the metal)
        if not hasattr(self,"connected_idx"): self.get_connected_idx()

        ## Creates the list of variables
        conn_idx     = self.connected_idx
        if debug > 0: print(f"LIGAND.SPLIT_LIGAND: {self.indices=}")
        if debug > 0: print(f"LIGAND.SPLIT_LIGAND: {conn_idx=}")
        conn_labels  = extract_from_list(conn_idx, self.labels, dimension=1)
        conn_coord   = extract_from_list(conn_idx, self.coord, dimension=1)
        conn_radii   = extract_from_list(conn_idx, self.radii, dimension=1)
        conn_atoms   = extract_from_list(conn_idx, self.atoms, dimension=1)
        if hasattr(self,"cov_factor"): blocklist = split_species(conn_labels, conn_coord, cov_factor=self.cov_factor, debug=debug)
        else:                          blocklist = split_species(conn_labels, conn_coord, debug=debug)

        ## Arranges Groups
        for b in blocklist:
            if debug > 0: print(f"LIGAND.SPLIT_LIGAND: block={b}")
            gr_indices = extract_from_list(b, conn_idx, dimension=1)
            if debug > 0: print(f"LIGAND.SPLIT_LIGAND: {gr_indices=}")
            gr_labels  = extract_from_list(gr_indices, conn_labels, dimension=1)
            gr_coord   = extract_from_list(gr_indices, conn_coord, dimension=1)
            gr_radii   = extract_from_list(gr_indices, conn_radii, dimension=1)
            gr_atoms   = extract_from_list(gr_indices, conn_atoms, dimension=1)
            newgroup = group(gr_labels, gr_coord, parent_indices=gr_indices, radii=gr_radii, parent=self)
            ## For group, we do not update parent. Atoms remain at the ligand level
            newgroup.set_atoms(atomlist=gr_atoms, overwrite_parent=False)
            newgroup.get_closest_metal()
            self.groups.append(newgroup)
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
        if not hasattr(self,"groups"):      self.split_ligand(debug=2)
        if debug > 0: print(f"LIGAND.Get_denticity: checking connectivity of ligand {self.formula}")
        if debug > 0: print(f"LIGAND.Get_denticity: initial connectivity is {len(self.connected_idx)}")
        self.denticity = 0
        for g in self.groups:
            if debug > 0: print(f"LIGAND.Get_denticity: checking denticity of group \n{g}")
            self.denticity += g.get_denticity(debug=debug)      ## A check is also performed at the group level
        if debug > 0: print(f"LIGAND.Get_denticity: final connectivity is {self.denticity}")
        return self.denticity

###############
#### GROUP ####
###############
class group(specie):
    def __init__(self, labels: list, coord: list, parent_indices: list=None, radii: list=None, parent: object=None) -> None:
        self.subtype = "group"
        specie.__init__(self, labels, coord, parent_indices, radii, parent)

   #######################################################
    def remove_atom(self, index: int, debug: int=0):
        if debug > 0: print(f"GROUP.REMOVE_ATOM: deleting atom {index=} from group")
        if index > self.natoms: return None
        if not hasattr(self,"atoms"): self.set_atoms()
        self.atoms.pop(index)
        self.labels.pop(index)
        self.coord.pop(index)
        self.parent_indices.pop(index)
        self.indices.pop(index)
        self.radii.pop(index)
        self.formula   = labels2formula(self.labels)
        self.eleccount = labels2electrons(self.labels)   ### Assuming neutral specie (so basically this is the sum of atomic numbers)
        self.natoms    = len(self.labels)
        self.iscomplex = any((elemdatabase.elementblock[l] == "d") or (elemdatabase.elementblock[l] == "f") for l in self.labels)
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
        for met in self.parent.parent.metals:
            bpos = np.array(met.coord)
            dist.append(np.linalg.norm(apos - bpos))
        # finds the closest Metal Atom (tgt)
        self.closest_metal = self.parent.parent.metals[np.argmin(dist)]
        return self.closest_metal

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
        return self.haptic_type

    #######################################################
    def get_denticity(self, debug: int=0):
        self.denticity = 0
        for a in self.atoms:
            self.denticity += a.mconnec
        return self.denticity

    #######################################################
    def __repr__(self):
        to_print = ""
        to_print += f'------------- Cell2mol GROUP Object --------------\n'
        to_print += specie.__repr__(self, indirect=True)
        to_print += '---------------------------------------------------\n'
        return to_print

###############
### BOND ######
###############
class bond(object):
    def __init__(self, atom1: object, atom2: object, bond_order: int=1):
        self.type       = "bond"
        self.version    = "0.1"
        self.atom1      = atom1
        self.atom2      = atom2
        self.order      = bond_order

    def __repr__(self):
        to_print += f'------------- Cell2mol BOND Object --------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        to_print += f' Atom 1                = {self.atom1.parent_index}\n'
        to_print += f' Atom 2                = {self.atom2.parent_index}\n'
        to_print += f' Bond Order            = {self.order}\n'
        to_print += '----------------------------------------------------\n'
        return to_print

###############
### ATOM ######
###############
class atom(object):
    def __init__(self, label: str, coord: list, parent_index: int=None, radii: float=None, parent: object=None, frac_coord: list=None) -> None:
        self.type          = "atom"
        self.version       = "0.1"
        self.label         = label
        self.coord         = coord
        self.atnum         = elemdatabase.elementnr[label]
        self.block         = elemdatabase.elementblock[label]

        if parent_index is not None:   self.parent_index = parent_index
        if radii is None:          self.radii = get_radii(label)
        else:                      self.radii = radii
        if parent is not None:     self.parent  = parent
        if parent is not None:     self.occurence = parent.get_occurrence(self)
        if frac_coord is not None: self.frac_coord = frac_coord

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
    def set_adjacency_parameters(self, cov_factor: float, metal_factor: float) -> None:
        self.cov_factor   = cov_factor
        self.metal_factor = metal_factor

    #######################################################
    def reset_charge(self) -> None:
        if hasattr(self,"charge"):    delattr(self,"charge")
        if hasattr(self,"poscharges"): delattr(self,"charge")

    #######################################################
    def set_charge(self, charge: int) -> None:
        self.charge = charge

    #######################################################
    def set_adjacencies(self, adjmat, madjmat, connectivity: int, metal_connectivity: int=0):
        self.connec  = int(connectivity)
        self.mconnec = int(metal_connectivity)
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
        for met in self.parent.parent.metals:
            bpos = np.array(met.coord)
            dist.append(np.linalg.norm(apos - bpos))
        self.closest_metal = self.parent.parent.metals[np.argmin(dist)]
        return self.closest_metal

    #######################################################
    def information(self, cov_factor: float, metal_factor: float) -> None:
        self.cov_factor   = cov_factor
        self.metal_factor = metal_factor

    #######################################################
    def __repr__(self, indirect: bool=False):
        to_print = ""
        if not indirect: to_print += f'------------- Cell2mol ATOM Object ----------------\n'
        to_print += f' Version                      = {self.version}\n'
        to_print += f' Type                         = {self.type}\n'
        if hasattr(self,'subtype'): to_print += f' Sub-Type                    = {self.subtype}\n'
        to_print += f' Label                        = {self.label}\n'
        to_print += f' Atomic Number                = {self.atnum}\n'
        #if self.parent_index is not None: to_print += f' Index in Parent ({self.parent.subtype})     = {self.parent_index}\n'
        if hasattr(self,"occurrence"): to_print += f' Occurrence in Parent         = {self.occurrence}\n'
        if hasattr(self,"mconnec"):    to_print += f' Metal Adjacency (mconnec)    = {self.mconnec}\n'
        if hasattr(self,"connec"):     to_print += f' Regular Adjacencies (connec) = {self.mconnec}\n'
        if hasattr(self,"charge"):     to_print += f' Atom Charge                  = {self.charge}\n'
        if not indirect: to_print += '----------------------------------------------------\n'
        return to_print

    #######################################################
    def reset_mconnec(self, index: int, value: int=0, debug: int=0):
        ## Careful. In group.atoms, the parent_index is that of the molecule
        if debug > 0: print(f"ATOM.RESET_MCONN: resetting mconnec for atom {self.label=} with {self.parent_index=}")
        self.mconnec = value                                                                             # Corrects data of atom object in self
        if hasattr(self,"parent"):   ## Normally, ligand level
            if debug > 0: print(f"ATOM.RESET_MCONN: updating ligand atoms and madjnum")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {self.parent.madjnum=}")
            self.parent.atoms[index].mconnec = value                                    # Corrects data of atom object in ligand class
            self.parent.madjnum[index] = value                                          # Corrects data in metal_adjacency number of the ligand class
            if debug > 0: print(f"ATOM.RESET_MCONN: final {self.parent.madjnum=}")
            if self.parent.subtype == "ligand": self.parent.get_connected_idx(debug=debug)
            if self.parent.subtype == "ligand": self.parent.get_connected_atoms(debug=debug)
            if hasattr(self.parent,"parent"):
                if debug > 0: print(f"ATOM.RESET_MCONN: updating molecule atoms and madjnum")
                self.parent.parent.atoms[self.parent.parent_indices[index]] = value   # Corrects data of atom object in molecule class
                self.parent.parent.madjnum[self.parent.parent_indices[index]] = value # Corrects data in metal_adjacency number of the molecule class

###############
#### METAL ####
###############
class metal(atom):
    def __init__(self, label: str, coord: list, parent_index: int=None, radii: float=None, parent: object=None, frac_coord: list=None) -> None:
        self.subtype      = "metal"
        atom.__init__(self, label, coord, parent_index=parent_index, radii=radii, parent=parent, frac_coord=frac_coord)

    #######################################################
    def get_valence_elec (self, m_ox: int):
        """ Count valence electrons for a given transition metal and metal oxidation state """
        v_elec = elemdatabase.valenceelectrons[self.label] - m_ox
        if v_elec >= 0 :  self.valence_elec = v_elec
        else :            self.valence_elec = elemdatabase.elementgroup[self.label] - m_ox
        return self.valence_elec

    #######################################################
    def get_coord_sphere(self, debug: int=0):
        if not hasattr(self,"parent"): 
            print("GET_COORD_SPHERE: no parent, returning None")
            return None
        if not hasattr(self.parent,"adjmat"): self.parent.get_adjmatrix()
        adjmat = self.parent.adjmat.copy()

        ## Cordination sphere defined as a collection of atoms
        self.coord_sphere = []
        for idx, at in enumerate(adjmat[self.parent_index]):
            if at >= 1: self.coord_sphere.append(self.parent.atoms[idx])
        return self.coord_sphere

    #######################################################
    def get_coord_sphere_formula(self, debug: int=0):
        if not hasattr(self,"coord_sphere"): self.get_coord_sphere(debug=1)
        self.coord_sphere_formula = labels2formula(list(at.label for at in self.coord_sphere))
        return self.coord_sphere_formula

    #######################################################
    def get_connected_groups(self, debug: int=0):
        # metal.groups will be used for the calculation of the relative metal radius
        # and define the coordination geometry of the metal /hapicitiy/ hapttype
        if not hasattr(self,"parent"): return None
        self.groups = []
        for group in self.parent.ligand.groups:
            tmplabels = self.label.copy()
            tmpcoord  = self.coord.copy()
            tmplabels.append(group.labels)
            tmpcoord.append(group.coord)
            isgood, tmpadjmat, tmpadjnum = get_adjmatrix(tmplabels, tmpcoord, metal_only=True)
            if isgood and any(tmpadjnum) > 0: self.groups.append(group)
        return self.groups

    ############
    def reset_charge(self):
        atom.reset_charge(self)     ## First uses the generic atom class function for itself
        if hasattr(self,"poscharges"):   delattr(self,"poscharge")

    def __repr__(self):
        to_print = ""
        to_print += f'------------- Cell2mol METAL Object --------------\n'
        to_print += atom.__repr__(self, indirect=True)
        to_print += '----------------------------------------------------\n'
        return to_print


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

    def get_adjmatrix(self):
        isgood, adjmat, adjnum = get_adjmatrix(self.labels, self.coord, 1.3, get_radii(self.labels))
        if isgood:
            self.adjmat = adjmat
            self.adjnum = adjnum
        else:
            self.adjmat = None
            self.adjnum = None
        return self.adjmat, self.adjnum

    def set_moleclist(self, moleclist: list) -> None:
        self.moleclist = moleclist

    def get_moleclist(self, blocklist=None):
        if not hasattr(self,"labels") or not hasattr(self,"pos"): return None
        if len(self.labels) == 0 or len(self.pos) == 0: return None
        cov_factor = 1.3

        if blocklist is None: blocklist = split_species(self.labels, self.pos, cov_factor=cov_factor)
        self.moleclist = []
        for b in blocklist:
            mol_labels  = extract_from_list(b, self.labels, dimension=1)
            mol_coords  = extract_from_list(b, self.coord, dimension=1)
            mol_indices = extract_from_list(b, range(self.natoms), dimension=1)
            newmolec    = molecule(mol_labels, mol_coords, parent_indices=mol_indices)
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

    #######################################################
    def reconstruct(self, cov_factor: float=None, metal_factor: float=None, debug: int=0):
        if not hasattr(self,"refmoleclist"): print("CELL.RECONSTRUCT. CELL missing list of reference molecules"); return None
        if not hasattr(self,"moleclist"): self.get_moleclist()
        blocklist    = self.moleclist.copy() # In principle, in moleclist now there are both fragments and molecules
        if cov_factor is None:   cov_factor   = self.refmoleclist[0].cov_factor
        if metal_factor is None: metal_factor = self.refmoleclist[0].metal_factor
        ## Classifies fragments
        for b in blocklist:
            if not hasattr(b,"frac_coord"):       b.get_fractional_coord(self.cellvec)
        moleclist, fraglist, Hlist = classify_fragments(blocklist, self.refmoleclist, debug=debug)

        ## Determines if Reconstruction is necessary
        if len(fraglist) > 0 or len(Hlist) > 0: self.is_fragmented = True
        else:                                   self.is_fragmented = False

        if not self.is_fragmented: return self.moleclist
        self.moleclist, Warning = fragments_reconstruct(moleclist,fraglist,Hlist,self.refmoleclist,self.cellvec,cov_factor,metal_factor)
        if Warning:      self.is_fragmented = True;  self.error_reconstruction = True
        else:            self.is_fragmented = False; self.error_reconstruction = False
        return self.moleclist

    def view(self):
        import plotly.graph_objects as go
        from Scope.Read_Write import set_scene
        from Scope.Elementdata import ElementData  
        elemdatabase = ElementData()

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

        set_scene(fig, np.array(self.coord))
        fig.show()

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
    assert hasattr(old_cell,"labels") and (hasattr(old_cell,"coord") or hasattr(old_cell,"pos"))
    assert hasattr(old_cell,"cellvec")
    assert hasattr(old_cell,"refcode")


    labels     = old_cell.labels
    refcode    = old_cell.refcode
    if   hasattr(old_cell,"coord"):       coord      = old_cell.coord
    elif hasattr(old_cell,"pos"):         coord      = old_cell.pos

    cellvec    = old_cell.cellvec
    if   hasattr(old_cell,"cellparam"):   cellparam  = old_cell.cellparam
    else:                                 cellparam  = cellvec_2_cellparam(old_cell.cellvec)

    if   hasattr(old_cell,"frac_coord"):  frac_coord = old_cell.frac_coord
    else:                                 frac_coord = cart2frac(coord, old_cell.cellvec)

    new_cell = cell(refcode, labels, coord, cellvec, cellparam)
    new_cell.set_fractional_coord(frac_coord)
    if debug > 0: print(f"IMPORT CELL: created new_cell {new_cell}")

    ## Moleclist
    moleclist = []
    for mol in old_cell.moleclist: 
        if debug > 0: print(f"IMPORT CELL: importing molecule {mol.formula} with charge {mol.totcharge}")
        new_mol = import_molecule(mol, parent=new_cell, debug=debug)
        moleclist.append(new_mol)
    new_cell.set_moleclist(moleclist)

    ## Refmoleclist
    new_cell.refmoleclist = []
    if hasattr(old_cell,"refmoleclist"):
        for rmol in old_cell.refmoleclist:
            new_cell.refmoleclist.append(import_molecule(rmol))
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
def import_molecule(mol: object, parent: object=None, debug: int=0) -> object:
    assert hasattr(mol,"labels") and (hasattr(mol,"coord") or hasattr(mol,"pos"))
    labels     = mol.labels
    if   hasattr(mol,"coord"):       coord      = mol.coord
    elif hasattr(mol,"pos"):         coord      = mol.pos

    if   hasattr(mol,"parent_indices"):     indices    = mol.parent_indices
    elif hasattr(mol,"atlist"):             indices    = mol.atlist
    else:                                   indices    = None

    if   hasattr(mol,"radii"):       radii      = mol.radii
    else:                            radii      = None          
 
    if parent is None: print(f"IMPORT MOLECULE {mol.formula}: parent is NONE")

    new_molec = molecule(labels, coord, indices, radii, parent)
    if debug > 0: print(f"IMPORT MOLEC: created new_molec with {new_molec.formula}")

    ## Charges
    if hasattr(mol,"totcharge") and hasattr(mol,"atcharge"):
        if debug > 0: print(f"IMPORT MOLEC: old molecule has total charge and atomic charges")
        new_molec.set_charges(mol.totcharge, mol.atcharge)
    elif hasattr(mol,"totcharge") and not hasattr(mol,"atcharge"):
        if debug > 0: print(f"IMPORT MOLEC: old molecule has total charge but no atomic charges")
        new_molec.set_charges(mol.totcharge)
    else:
        if debug > 0: print(f"IMPORT MOLEC: old molecule has no total charge nor atomic charges")

    ## Smiles
    if   hasattr(mol,"Smiles"): new_molec.smiles = mol.Smiles
    elif hasattr(mol,"smiles"): new_molec.smiles = mol.smiles        

    ## Substructures
    if not hasattr(mol,"ligandlist") or not hasattr(mol,"metalist"):  new_molec.split_complex()
    if hasattr(mol,"ligandlist"): 
        ligands = []
        for lig in mol.ligandlist: 
            if debug > 0: print(f"IMPORT MOLEC: old molecule has ligand {lig.formula}")
            new_lig = import_ligand(lig, parent=new_molec, debug=debug)
            ligands.append(new_lig)
        new_molec.ligands = ligands
    if hasattr(mol,"metalist"): 
        metals = []
        for met in mol.metalist: 
            if debug > 0: print(f"IMPORT MOLEC: old molecule has metal {met.label}")
            new_atom = import_atom(met, parent=new_molec, debug=debug)
            metals.append(new_atom)
        new_molec.metals = metals

    ## Atoms
    if not hasattr(mol,"atoms"):  
        if debug > 0: print(f"IMPORT MOLEC: old molecule has no atoms")
        new_molec.set_atoms()
    else: 
        atoms = []
        for at in mol.atoms: 
            if debug > 0: print(f"IMPORT MOLEC: old molecule has atom {at.label}")
            new_atom = import_atom(at, parent=new_molec, debug=debug)
            atoms.append(new_atom)
        new_molec.set_atoms(atoms)

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

    new_ligand = ligand(labels, coord, indices, radii, parent)
    if debug > 0: print(f"IMPORT LIGAND: created new_ligand with {new_ligand.formula}")
    
    ## Charges
    if hasattr(lig,"totcharge") and hasattr(lig,"atcharge"):
        new_ligand.set_charges(lig.totcharge, lig.atcharge)
    elif hasattr(lig,"totcharge") and not hasattr(lig,"atcharge"):
        new_ligand.set_charges(lig.totcharge)

    ## Smiles
    if   hasattr(lig,"Smiles"): new_ligand.smiles = lig.Smiles
    elif hasattr(lig,"smiles"): new_ligand.smiles = lig.smiles     

    ## Rdkit Object
    if   hasattr(lig,"object"): new_ligand.rdkit_obj = lig.object
    
    ## Substructures
    if not hasattr(lig,"grouplist"):  new_ligand.split_ligand()
    else: 
        groups = []
        for gr in lig.grouplist: 
            new_group = import_group(gr, parent=new_ligand, debug=debug)
            groups.append(new_group)
    new_ligand.groups = groups
    
    ## Atoms
    if not hasattr(lig,"atoms"):  new_ligand.set_atoms()
    else: 
        atoms = []
        for at in lig.atoms: 
            new_atom = import_atom(at, parent=new_ligand, debug=debug)
            atoms.append(new_atom)
        new_ligand.set_atoms(atoms)
            
    return new_ligand

################################
def import_group(old_group: object, parent: object=None, debug: int=0) -> object:
    assert hasattr(old_group,"atlist") or hasattr(old_group,"indices")
    
    if   hasattr(old_group,"labels"):      labels     = old_group.labels
    elif parent is not None:
        if hasattr(parent,"labels"):       labels     = extract_from_list(old_group.atlist, parent.labels, dimension=1)
        
    if   hasattr(old_group,"coord"):       coord      = old_group.coord
    elif hasattr(old_group,"pos"):         coord      = old_group.pos
    elif parent is not None:
        if hasattr(parent,"coord"):        coord     = extract_from_list(old_group.atlist, parent.coord, dimension=1)
        elif hasattr(parent,"coord"):      coord     = extract_from_list(old_group.atlist, parent.pos, dimension=1)

    if   hasattr(old_group,"parent_indices"):     indices    = old_group.parent_indices
    elif hasattr(old_group,"atlist"):             indices    = old_group.atlist
    else:                                         indices    = None

    if   hasattr(old_group,"radii"):       radii      = old_group.radii
    else:                                  radii      = None          

    if parent is None: print("IMPORT GROUP: parent is NONE")

    new_group = group(labels, coord, indices, radii, parent)
    if debug > 0: print(f"IMPORT GROUP: created new_group with {new_group.formula}")
    
    ## Charges
    if hasattr(old_group,"totcharge") and hasattr(old_group,"atcharge"):
        new_group.set_charges(old_group.totcharge, old_group.atcharge)
    elif hasattr(old_group,"totcharge") and not hasattr(old_group,"atcharge"):
        new_group.set_charges(old_group.totcharge)

    ## Smiles
    if   hasattr(old_group,"Smiles"): new_group.smiles = old_group.Smiles
    elif hasattr(old_group,"smiles"): new_group.smiles = old_group.smiles     
    
    ## Atoms
    if not hasattr(old_group,"atoms"):  new_group.set_atoms()
    else: 
        atoms = []
        for at in old_group.atoms: 
            new_atom = import_atom(at, parent=new_group, debug=debug)
            atoms.append(new_atom)
        new_group.set_atoms(atoms)
            
    return new_group

################################
def import_atom(old_atom: object, parent: object=None, debug: int=0) -> object:
    assert hasattr(old_atom,"label") and (hasattr(old_atom,"coord") or hasattr(old_atom,"pos"))
    label     = old_atom.label
    
    if   hasattr(old_atom,"coord"):      coord      = old_atom.coord
    elif hasattr(old_atom,"pos"):        coord      = old_atom.pos

    if   hasattr(old_atom,"parent_index"):      index      = old_atom.parent_index
    elif hasattr(old_atom,"atlist"):            index      = old_atom.atlist
    else:                                       index      = None

    if   hasattr(old_atom,"radii"):      radii      = old_atom.radii
    else:                                radii      = None          

    if   hasattr(old_atom,"block"):      block      = old_atom.block
    else:                                block      = elemdatabase.elementblock[label]    

    if parent is None: print("IMPORT ATOM: parent is NONE")
    
    if block == 'd' or block == 'f': 
        new_atom = metal(label, coord, parent_index=index, radii=radii, parent=parent)
        if debug > 0: print(f"IMPORT ATOM: created new_metal {new_atom.label}")
    else:                            
        new_atom = atom(label, coord, parent_index=index, radii=radii, parent=parent)
        if debug > 0: print(f"IMPORT ATOM: created new_atom {new_atom.label}")
    
    ## Charge. For some weird reason, charge is "charge" in metals, and "totcharge" in regular atoms
    if hasattr(old_atom,"charge"):       
        if type(old_atom.charge) == (int):      new_atom.set_charge(old_atom.charge)
        elif hasattr(old_atom,"totcharge"):     
            if type(old_atom.totcharge) == int: new_atom.set_charge(old_atom.totcharge)
            
    return new_atom
