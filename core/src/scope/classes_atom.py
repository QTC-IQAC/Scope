###########################################################
####  Contains the Atom & Bond Classes and Sub-Classes ####
###########################################################

import numpy as np
from scope.connectivity               import * 
from scope.geometry                   import get_dist
from scope.elementdata                import ElementData
from scope.operations.dicts_and_lists import extract_from_list 
elemdatabase = ElementData()

############
### ATOM ###
############
class Atom(object):
    """
    Represent a single atom in SCOPE.

    Attributes:
        object_type (str):              Object category (`"atom"`).
        object_subtype (str):           Atom subtype.
        label (str):                    Atomic symbol.
        coord (list):                   Cartesian coordinates.
        parents (list):                 Parent species that contain the atom.
        bonds (list):                   Bonds attached to the atom.

    Methods:
        set_charge():                   Store an atomic charge.
        set_spin():                     Store an atomic spin.
        add_parent():                   Link the atom to a parent object.
        set_adjacencies():              Store covalent and metal adjacency data.
        check_connectivity():           Test connectivity against another atom.
    """
    def __init__(self, label: str, coord: list, frac_coord: list=None, radii: float=None) -> None:
        self.object_type          = "atom"
        self.object_subtype       = "atom"
        self.version              = "1.0"
        self.origin               = "created"
        self.label                = label
        self.coord                = coord
        self.atnum                = elemdatabase.elementnr[label]
        self.block                = elemdatabase.elementblock[label]
        self.parents              = []
        self.parents_index        = []
        self.bonds                = []
        self.formula              = label

        if frac_coord is not None:        self.frac_coord = frac_coord
        if radii is None:                 self.radii = get_radii(label)
        else:                             self.radii = radii

    #####################
    ## Spin and Charge ##
    #####################
    def reset_spin(self) -> None:
        self.spin = int(0)

    def reset_charge(self) -> None:
        self.charge = int(0)

    def set_charge(self, charge: int) -> None:
        self.charge = int(charge)

    def set_spin(self, spin: int) -> None:
        self.spin           = int(spin)
        self.multiplicity   = int(spin + 1)
        self.ms             = float(spin/2)

    ###########
    ## Other ##
    ###########
    def __eq__(self, other, check_coordinates: bool=False, debug: int=0):
        if not isinstance(other, type(self)): return False
        if debug > 0:
            print("Comparing Atoms")
            print(self)
            print(other)
        # Compares Species, Coordinates, Charge and Spin
        if (self.label != other.label): return False
        if hasattr(self,"charge") and hasattr(other,"charge"):
            if (self.charge != other.charge): return False
        if hasattr(self,"spin") and hasattr(other,"spin"):
            if (self.spin != other.spin): return False
        if hasattr(self,"adjnum") and hasattr(other,"adjnum"):
            if (self.adjnum != other.adjnum): return False
        if hasattr(self,"madjnum") and hasattr(other,"madjnum"):
            if (self.madjnum != other.madjnum): return False
        if check_coordinates:
            if (self.coord[0] != other.coord[0]): return False
            if (self.coord[1] != other.coord[1]): return False
            if (self.coord[2] != other.coord[2]): return False
        return True

    def get_decorated_label(self, typ: str="spin"):
        ## Function for Quantum Espresso Inputs
        if   typ.lower() == "spin" and hasattr(self,"spin"):                   return self.label + str(self.spin)
        elif typ.lower() == "multiplicity" and hasattr(self,"multiplicity"):   return self.label + str(self.multiplicity)
        else: print("ATOM.GET_MOD_LABEL: {typ=} is not implemented")

    #############
    ## Parents ##
    #############
    def add_parent(self, parent: object, index: int, overwrite: bool=True, debug: int=0):
        ## associates a parent specie to self. The Atom indices of self in parent are given in "indices"
        ## if parent of the same subtype already in self.parent then it is overwritten
        ## this is to avoid having a substructure (e.g. a ligand) in more than one superstructure (e.g. a molecule) 
        append = True
        for idx, p in enumerate(self.parents):
            if p.object_subtype == parent.object_subtype:
                if overwrite: 
                    if debug > 0: print(f"ATOM.ADD_PARENT: Parent with same subtype {parent.object_subtype} exists, but overwrite={overwrite}. Overwriting")
                    self.parents[idx]         = parent
                    self.parents_index[idx]   = index
                else:
                    if debug > 0: print(f"ATOM.ADD_PARENT: Parent with same subtype {parent.object_subtype} exists, so new won't be added")
                append = False
        if append: 
            self.parents.append(parent)
            self.parents_index.append(index)
            if debug > 0: print(f"ATOM.ADD_PARENT: added parent with subtype {parent.object_subtype}. Atom index is {index=}")

    ######
    def check_parent(self, subtype: str, search_by: str="subtype"):
        ## checks if parent of a given subtype exists
        for p in self.parents:
            if search_by.lower() == "subtype":
                if p.object_subtype == subtype: return True
            elif search_by.lower() == "type":
                if p.object_type == subtype: return True
            else:
                raise ValueError(f"ATOM.CHECK_PARENT: unknown {search_by=}")
        return False

    ######
    def get_parent(self, subtype: str, search_by: str="subtype"):
        ## retrieves parent of a given subtype 
        for p in self.parents:
            if search_by.lower() == "subtype":
                if p.object_subtype == subtype: return p
            elif search_by.lower() == "type":
                if p.object_type == subtype: return p
            else:
                raise ValueError(f"ATOM.GET_PARENT: unknown {search_by=}")
        return None
    
    ######
    def get_parent_index(self, subtype: str, search_by: str="subtype"):
        ## retrieves parent of a given subtype 
        for idx, p in enumerate(self.parents):
            if search_by.lower() == "subtype":
                if p.object_subtype == subtype: return self.parents_index[idx]
            elif search_by.lower() == "type":
                if p.object_type == subtype: return self.parents_index[idx]
            else:
                raise ValueError(f"ATOM.GET_PARENT_INDEX: unknown {search_by=}")
        return None

    ##################
    ## Connectivity ##
    ##################
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

    ######
    def check_connectivity(self, other: object, debug: int=0):
        ## Checks whether two Atoms are connected (through the adjacency)
        if not isinstance(other, type(self)): return False
        labels = list([self.label,other.label]) 
        coords = list([self.coord,other.coord]) 
        isgood, adjmat, adjnum = get_adjmatrix(labels, coords)
        if isgood and adjnum[0] > 0: return True
        else:                        return False

    ######
    def set_adjacencies(self, adjmat, madjmat, connectivity: int, metal_connectivity: int=0):
        self.adjnum  = int(connectivity)
        self.madjnum = int(metal_connectivity)
        self.adjacency       = []
        self.metal_adjacency = []
        for idx, c in enumerate(adjmat):   ## The Atom only receives one row of adjmat, so this is not a matrix anymore. Keep in mind that the idx are the indices of parent
            if c > 0: self.adjacency.append(idx)
        for idx, c in enumerate(madjmat):  ## The Atom only receives one row of madjmat, so this is not a matrix anymore
            if c > 0: self.metal_adjacency.append(idx)

    ######
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

    ######
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

    ######
    def set_factors(self, cov_factor: float=1.3, metal_factor: float=1.0) -> None:
        self.cov_factor   = cov_factor
        self.metal_factor = metal_factor

    ######
    def reset_madjnum(self, met, diff: int=-1, debug: int=0):
        if debug > 0: print(f"ATOM.RESET_MCONN: resetting madjnum (and connec) for atom {self.label=}")
        if debug > 0: print(f"ATOM.RESET_MCONN: initial {self.adjnum=} {self.madjnum=}")
        if debug > 0: print(f"ATOM.RESET_MCONN: initial = {self.adjacency=} {self.metal_adjacency=}")
        self.madjnum += diff
        self.adjnum  += diff

        if debug > 0:  print(f"ATOM.RESET_MCONN: initial {met.adjnum=} {met.madjnum=}")
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
            if debug > 0: print(f"ATOM.RESET_MCONN: updating ligand atoms and adjnum")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.adjnum=} {len(lig.adjnum)}") 
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.atoms[lig_idx].adjnum=} {lig.atoms[lig_idx].madjnum=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.madjnum[lig_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.adjnum[lig_idx]=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {met.adjnum=} {met.madjnum=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {lig.madjmat[lig_idx]=} {lig.adjmat[lig_idx]=}")
            # Correct Ligand Data
            lig.madjnum[lig_idx] += diff                    # Corrects data in metal_adjacency number of the ligand class
            lig.adjnum[lig_idx]  += diff                    # Corrects data in adjacency number of the ligand class

            lig.atoms[lig_idx].set_adjacencies(lig.adjmat[lig_idx], lig.madjmat[lig_idx], lig.adjnum[lig_idx], lig.madjnum[lig_idx])

            # we should delete the adjacencies, but not a priority 
            if debug > 0: print(f"ATOM.RESET_MCONN: final {lig.madjnum=}")
            if debug > 0: print(f"ATOM.RESET_MCONN: final {lig.adjnum=}")  
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
            if debug > 0: print(f"ATOM.RESET_MCONN: updating molecule atoms and adjnum")
            if debug > 0: print(f"ATOM.RESET_MCONN: initial {mol.adjnum=} {len(mol.adjnum)}") 
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

    ###########
    ## Bonds ##
    ###########
    def find_bond(self, end_idx: int, parent_str: str, debug: int=0):
        ## Finds bond object in self.bonds. It needs an index, and the parent string to identify the parent.
        ## Remember indices change between parents
        if not hasattr(self,"bonds"): self.bonds = []
        for b in self.bonds:
            if   b.atom1.get_parent_index(parent_str) == end_idx and b.atom2.get_parent_index(parent_str) == self.get_parent_index(parent_str):
                return b
            elif b.atom2.get_parent_index(parent_str) == end_idx and b.atom1.get_parent_index(parent_str) == self.get_parent_index(parent_str):
                return b
        return None
            
    ######
    def add_bond(self, newbond: object, debug: int=0):
        if not hasattr(self,"bonds"): self.bonds = []
        at1 = newbond.atom1
        at2 = newbond.atom2
        found = False
        for b in self.bonds:
            issame1 = b.atom1.__eq__(at1, check_coordinates=True)
            issame2 = b.atom2.__eq__(at2, check_coordinates=True)
            issame3 = b.atom1.__eq__(at2, check_coordinates=True)
            issame4 = b.atom2.__eq__(at1, check_coordinates=True)
            if (issame1 and issame2) or (issame3 and issame4):
                if debug > 0: print(f"ATOM.ADD_BOND found the same bond with atoms:") 
                if debug > 0: print(f"atom1: {b.atom1}") 
                if debug > 0: print(f"atom2: {b.atom2}") 
                found = True    ### It means that the same bond has already been defined
        if not found: self.bonds.append(newbond)

    ###################
    ## Visualization ##
    ###################
    def __repr__(self, indirect: bool=False):
        to_print = ''
        if not indirect: to_print += '------------------------------------------------\n'
        if not indirect: to_print += '------------- SCOPE ATOM Object ----------------\n'
        if not indirect: to_print += '------------------------------------------------\n'
        to_print += f' Version                      = {self.version}\n'
        to_print += f' Type                         = {self.object_type}\n'
        if hasattr(self,'object_subtype'): to_print += f' Sub-Type                     = {self.object_subtype}\n'
        to_print += f' Label                        = {self.label}\n'
        to_print += f' Atomic Number                = {self.atnum}\n'
        if hasattr(self,"charge"):     to_print += f' Atom Charge                  = {self.charge}\n'
        if hasattr(self,"spin"):       to_print += f' Spin (alpha - beta)          = {self.spin}\n'
        # Adjacency and Metal Adjacency
        if hasattr(self,"mconnec"):    to_print += f' Metal Adjacencies (mconnec)  = {self.mconnec}\n'
        elif hasattr(self,"madjnum"):  to_print += f' Metal Adjacencies (madjnum)  = {self.madjnum}\n'
        else:                          to_print += f' No Metal Adjacency Info\n'
        if hasattr(self,"connec"):     to_print += f' Regular Adjacencies (connec) = {self.connec}\n'
        elif hasattr(self,"adjnum"):   to_print += f' Regular Adjacencies (adjnum) = {self.adjnum}\n'
        else:                          to_print += f' No Adjacency Info\n'
        if not indirect: to_print += '\n'
        return to_print

###############
#### METAL ####
###############
class Metal(Atom):
    """
    Represent a metal atom with coordination-sphere helpers.

    Attributes:
        object_subtype (str):           Atom subtype (`"metal"`).
        coord_sphere (list):            Neighboring atoms around the metal.
        coord_sphere_formula (str):     Formula of the coordination sphere.

    Methods:
        get_valence_elec():             Estimate valence electrons from oxidation state.
        get_coord_sphere():             Retrieve atoms bound to the metal.
        get_coord_sphere_formula():     Summarize the coordination sphere formula.
        get_connected_groups():         Find ligand groups attached to the metal.
        get_cshm():                     Compute a continuous shape measure.
    """
    def __init__(self, label: str, coord: list, frac_coord: list=None, radii: float=None) -> None:
        Atom.__init__(self, label, coord, frac_coord=frac_coord, radii=radii)
        self.object_subtype = "metal"

    ###########
    ## Other ##
    ###########
    def __eq__(self, other, check_coordinates: bool=False, debug: int=0):
        if not isinstance(other, type(self)): return False
        same_atom = super().__eq__(other, check_coordinates=check_coordinates)
        if not same_atom: return False
        if not hasattr(self,"coord_sphere_formula"): self.get_coord_sphere_formula()
        if not hasattr(other,"coord_sphere_formula"): other.get_coord_sphere_formula()
        if (self.coord_sphere_formula != other.coord_sphere_formula):
            if debug > 0: print("COMPARE_METALS. Different coordination sphere")
            if debug > 0: print(self.coord_sphere_formula)
            if debug > 0: print(other.coord_sphere_formula)
            return False
        return True 

    ###################
    ## Visualization ##
    ###################
    def __repr__(self):
        to_print =  '-----------------------------------------------\n'
        to_print += '------------- SCOPE METAL Object --------------\n'
        to_print += '-----------------------------------------------\n'
        to_print += Atom.__repr__(self, indirect=True)
        to_print += '\n'
        return to_print

    ######
    def get_valence_elec (self, m_ox: int):
        """ Count valence electrons for a given transition metal and metal oxidation state """
        v_elec = elemdatabase.valenceelectrons[self.label] - m_ox
        if v_elec >= 0 :  self.valence_elec = v_elec
        else :            self.valence_elec = elemdatabase.elementgroup[self.label] - m_ox
        return self.valence_elec

    #########################
    ## Coordination Sphere ##
    #########################
    def get_coord_sphere_idx(self, debug: int=0):
        if not self.check_parent("molecule"): 
            print(f"METAL.Get_coord_sphere_idx. Metal {self.label} does not have parent molecule")
            self.coord_sphere_idx = None
        else:
            mol = self.get_parent("molecule")
            pidx = self.get_parent_index("molecule")
            if not hasattr(mol,"adjmat"): mol.get_adjmatrix()
            adjmat = mol.adjmat.copy()
            ## Collect Cordination sphere idx as a collection of indice
            self.coord_sphere_idx = [idx for idx, at in enumerate(adjmat[pidx]) if at >= 1]
        return self.coord_sphere_idx

    ######
    def get_coord_sphere(self, debug: int=0):
        if not hasattr(self,"coord_sphere_idx"): self.get_coord_sphere_idx(debug=debug)
        if self.coord_sphere_idx is None:
            self.coord_sphere = None 
        else:
            mol = self.get_parent("molecule")
            ## Collect Cordination sphere defined as a collection of atoms
            self.coord_sphere = [at for idx, at in enumerate(mol.atoms) if idx in self.coord_sphere_idx]
        return self.coord_sphere

    ######
    def get_coord_sphere_formula(self, debug: int=0):
        if not hasattr(self,"coord_sphere"): self.get_coord_sphere(debug=debug)
        if self.coord_sphere is None:
            self.coord_sphere_formula = None 
        else:
            self.coord_sphere_formula = labels2formula(list([at.label for at in self.coord_sphere])) 
            if debug > 0: print(f"METAL.Get_coord_sphere_formula: {self.get_parent_index('molecule')} {self.label} {self.coord_sphere_formula}")
        return self.coord_sphere_formula 

    ##################
    ## Connectivity ##
    ##################
    def get_connected_groups(self, debug: int=0):
        from scope.connectivity import split_group
        # metal.groups will be used for the calculation of the relative metal radius 
        # and define the coordination geometry of the metal hapticitiy and hapttype    
        if not self.check_parent("molecule"): return None
        mol = self.get_parent("molecule")
        self.groups = []
        for lig in mol.ligands:
            for group in lig.groups:
                ligand_indices = [ a.get_parent_index("ligand") for a in group.atoms ]
                tmplabels = []
                tmpcoord  = []
                tmplabels.append(self.label)
                tmpcoord.append(self.coord)
                tmplabels.extend(group.labels)
                tmpcoord.extend(group.coord)
                isgood, tmpadjmat, tmpadjnum = get_adjmatrix(tmplabels, tmpcoord, metal_only=True)
                # if isgood and any(tmpadjnum) > 0: self.groups.append(group)
                if isgood:
                    if all(tmpadjnum[1:]): 
                        self.groups.append(group)
                    elif any(tmpadjnum[1:]): 
                        if debug > 0: print(f"METAL.GET_CONNECTED_GROUPS: Metal {self.label} is connected to {group.formula} but not all atoms are connected")
                        conn_idx = [ idx for idx, num in enumerate(tmpadjnum[1:]) if num == 1 ]
                        conn_ligand_indices = [ ligand_indices[idx] for idx, num in enumerate(tmpadjnum[1:]) if num == 1 ]
                        if debug > 0: print(f"METAL.GET_CONNECTED_GROUPS: get_connected_groups {tmpadjnum[1:]=} {conn_idx=} {conn_ligand_indices=} {ligand_indices=}")
                        splitted_groups = split_group(group, conn_idx, conn_ligand_indices, debug=debug)
                        for g in splitted_groups:
                            self.groups.append(g)
                            if debug > 0: print(f"METAL.GET_CONNECTED_GROUPS: Metal {self.label} is connected to {g.formula}")
                    else:
                        if debug > 0: print(f"METAL.GET_CONNECTED_GROUPS: Metal {self.label} is not connected to {group.formula}")
        return self.groups

    ######
    def get_connected_metals(self, debug: int=0):
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
            isgood, tmpadjmat, tmpadjnum = get_adjmatrix(tmplabels, tmpcoord, metal_only=True)
            if isgood:
                if all(tmpadjnum[1:]): 
                    self.metals.append(met)
                else:
                    if debug > 0: print(f"METAL.GET_CONNECTED_METALS: Metal {self.label} is not connected to {met.label}")
        if debug >= 2 : print(f"METAL.GET_CONNECTED_METALS: {self.label} connected to {len(self.metals)} metals {[m.label for m in self.metals]}")
        return self.metals

    ##########################
    ## Geometric Parameters ##
    ##########################
    def get_cshm(self, ref_shape: str='OC-6', overwrite: bool=False, debug: int=0):
        try:
            import cosymlib as cml
        except ImportError: 
            raise ImportError("This function requires cosymlib, which is currently not installed. It was removed from the dependencies of scope as it restrained the versions of other software. Please install it manually if you need this functionality.")
        from scope.cshm import get_cshm_ref

        # Prepares coordiantes of the coordination sphere
        if hasattr(self, "cshm") and not overwrite: return self.cshm
        coord_sphere_coords = [self.coord]
        coord_sphere_coords += [at.coord for at in self.get_coord_sphere()]
        # Gets Coordinates of the reference shape
        ref_coords = get_cshm_ref(ref_shape)

        class CustomShape(cml.shape.Shape):
            def get_positions(self):
                return self._coordinates

        # Wrap in shape objects
        shape_current = CustomShape(np.array(coord_sphere_coords))
        shape_reference = CustomShape(np.array(ref_coords))
        self.cshm = shape_current.measure(shape_reference)
        return self.cshm

############
### BOND ###
############
class Bond(object):
    """
    Represent a bond between two atoms.

    Attributes:
        object_type (str):              Object category (`"bond"`).
        object_subtype (str):           Bond subtype.
        atom1 (object):                 First bonded atom.
        atom2 (object):                 Second bonded atom.
        order (int | float):            Formal bond order.
        distance (float):               Interatomic distance.

    Methods:
        __repr__():                     Return a formatted description of the bond.
    """
    def __init__(self, atom1: object, atom2: object, bond_order: int=1, subtype: str="intraspecie"):
        self.object_type = "bond"
        self.object_subtype = subtype
        self.version    = "1.0"
        self.atom1      = atom1
        self.atom2      = atom2
        self.order      = bond_order
        self.distance   = np.linalg.norm(np.array(atom1.coord) - np.array(atom2.coord))

    def __repr__(self):
        to_print =  '----------------------------------------------\n'
        to_print += '------------- SCOPE BOND Object --------------\n'
        to_print += '----------------------------------------------\n'
        to_print += f' Version                 = {self.version}\n'
        to_print += f' Type                    = {self.object_type}\n'
        to_print += f' Atom1 Index in Molecule = {self.atom1.get_parent_index("molecule")}\n'
        to_print += f' Atom2 Index in Molecule = {self.atom2.get_parent_index("molecule")}\n'
        to_print += f' Label Atom1             = {self.atom1.label}\n'
        to_print += f' Label Atom2             = {self.atom2.label}\n'
        to_print += f' Formal Bond Order       = {self.order}\n'
        to_print += f' Distance                = {round(self.distance,3)} Angstrom\n'
        to_print += '\n'
        return to_print

###############
### IMPORTS ###
###############
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
        new_atom = Metal(label, coord, radii=radii)
    else:                            
        if debug > 0: print(f"IMPORT ATOM: importing atom {old_atom.label}")
        new_atom = Atom(label, coord, radii=radii)
    new_atom.origin = "import_atom"
    
    ## Parents
    if parent is not None:
        new_atom.add_parent(parent, index, overwrite=False, debug=debug)
    if debug > 0: print(f"IMPORT ATOM: parent {parent.object_subtype=} added with {index=}") 
    else:
        if debug > 0: print(f"IMPORT ATOM: parent is None") 
    
    ## Charge. For some weird reason, charge is "charge" in regular atoms, and "totcharge" in metals
    ## Also. In cell2mol version 1, atoms do not carry charge. Metals in metalist yes
    if debug > 0: print(f"IMPORT ATOM: now trying to import charge") 
    if hasattr(old_atom,"charge"):       
        if type(old_atom.charge) == int:    
            new_atom.set_charge(old_atom.charge)
            if new_atom.object_subtype == "atom"  and debug > 0: print(f"IMPORT ATOM: atom charge set to {new_atom.charge}") 
            if new_atom.object_subtype == "metal" and debug > 0: print(f"IMPORT ATOM: metal charge set to {new_atom.charge}") 

    elif hasattr(old_atom,"totcharge"):     
        if type(old_atom.totcharge) == int: 
            new_atom.set_charge(old_atom.totcharge)
            if new_atom.object_subtype == "atom"  and debug > 0: print(f"IMPORT ATOM: atom charge set to {new_atom.charge}") 
            if new_atom.object_subtype == "metal" and debug > 0: print(f"IMPORT ATOM: metal charge set to {new_atom.charge}") 
    else:
        if debug > 0: print(f"IMPORT ATOM: no charge found for old atom {old_atom.label}. Setting Default: charge=0")
        new_atom.set_charge(int(0))

    ## Spin.
    if debug > 0: print(f"IMPORT ATOM: now trying to import spin") 
    if hasattr(old_atom,"spin"):
        if type(old_atom.spin) == int:    
            new_atom.set_spin(old_atom.spin)
    else:
        if debug > 0: print(f"IMPORT ATOM: no spin found for old atom {old_atom.label}. Setting Default: spin=0")
        new_atom.set_spin(int(0))

    ## Connectivity
    if debug > 0: print(f"IMPORT ATOM: inheriting connectivity for {new_atom.object_subtype=}")
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
