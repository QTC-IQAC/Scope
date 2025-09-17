##############################
##### General CELL Class #####
##############################

import numpy as np
from .Connectivity import *
from .Geometry import cellparam_2_cellvec, cellvec_2_cellparam, cart2frac
from .Elementdata import ElementData
elemdatabase = ElementData()

##############
#### CELL ####
##############
class cell(object):
    def __init__(self, name: str, labels: list, coord: list, cell_vector: list=None, cell_param: list=None) -> None:
        self.version              = "1.0"
        self.type                 = "cell"
        self.subtype              = "cell"
        self.origin               = "created"
        self.name                 = name
        self.labels               = labels 
        self.coord                = coord
        self.formula              = labels2formula(labels)
        self.natoms               = len(labels)

        ## Gets Cell Parameters, Vectors and Fractional Coordinates
        if   cell_vector is None and cell_param is None:
            raise ValueError("Either cell_vector or cell_param must be provided")
        elif cell_vector is None and cell_param is not None:
            self.cell_param       = cell_param
            self.cell_vector      = cellparam_2_cellvec(cell_param)
        elif cell_vector is not None and cell_param is None:
            self.cell_vector      = cell_vector
            self.cell_param       = cellvec_2_cellparam(cell_vector)

        self.frac_coord           = cart2frac(self.coord, self.cell_vector)

    def __repr__(self, indirect: bool=False):
        to_print = ''   
        if not indirect: to_print += '-------------------------------\n'
        if not indirect: to_print += '   >>> SCOPE CELL Object >>>   \n'
        if not indirect: to_print += '-------------------------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        to_print += f' Name                  = {self.name}\n'
        to_print += f' Num Atoms             = {self.natoms}\n'
        to_print += f' Cell Parameters a:c   = {self.cell_param[0:3]}\n'
        to_print += f' Cell Parameters al:ga = {self.cell_param[3:6]}\n'
        if hasattr(self,"moleclist"):  
            to_print += f' # Molecules:          = {len(self.moleclist)}\n'
            to_print += f' With Formulae:                               \n'
            for idx, m in enumerate(self.moleclist):
                to_print += f'    {idx}: {m.formula} \n'
        to_print += '-------------------------------\n'
        if hasattr(self,"refmoleclist"):
            to_print += f' # of Ref Molecules:   = {len(self.refmoleclist)}\n'
            to_print += f' With Formulae:                                  \n'
            for idx, ref in enumerate(self.refmoleclist):
                to_print += f'    {idx}: {ref.formula} \n'
        if not indirect: to_print += '\n'
        return to_print

    ######
    def save(self, filepath):
        save_binary(self, filepath)

    ######
    def associate_cif(self, cif: object) -> None:
        self.cif       = cif 
        return self.cif

    ######
    def set_path(self, path: str) -> None:
        self.path = path

    ######
    def fix_cell_coord(self, debug: int=0) -> None:
        ## In cell2mol, the cell object does not have the coordinates of the reconstructed cell. 
        ## However, the molecule and atom objects are updated (i.e. reconstructed). We use this info to update the cell
        if not hasattr(self,"moleclist"): self.get_moleclist()

        self.labels = []
        self.coord  = []
        indices     = []
        for mol in self.moleclist:
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

    ## This function mimics the specie-class function with the same name
    def check_parent(self, subtype):
        if subtype == "cell": return True
        else:                 return False

    ## This function mimics the specie-class function with the same name
    def get_parent(self, subtype):
        if subtype == "cell": return self
        else:                 return None

    ## This function mimics the specie-class function with the same name
    def get_parent_indices(self, subtype: str):
        if subtype == "cell": return list(range(0,self.natoms))
        else:                 return None

    ######
    def get_adjmatrix(self):
        isgood, adjmat, adjnum = get_adjmatrix(self.labels, self.coord, self.cov_factor, self.metal_factor, get_radii(self.labels))
        if isgood:
            self.adjmat = adjmat
            self.adjnum = adjnum
        else:
            self.adjmat = None
            self.adjnum = None
        return self.adjmat, self.adjnum

    ######
    def set_moleclist(self, moleclist: list) -> None:
        self.moleclist = moleclist

    ######
    def check_fragmentation(self, reconstruct: bool = False, debug: int=0):

        if not hasattr(self,"moleclist"): self.get_moleclist()
        self.fragmented = False

        # First comparison with ref_moleclist
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

    ######
    def get_moleclist(self, overwrite: bool=False, cov_factor: float=1.3, metal_factor: float=1.0, debug: int=0):

        ## Overwrite and Warning
        if not overwrite and hasattr(self,"moleclist"): 
            if debug > 0: print(f"CELL.GET_MOLECLIST. Moleclist already exists and default is overwrite=False")
            return self.moleclist

        ## Security
        if not hasattr(self,"labels") or not hasattr(self,"coord"): 
            if debug > 0: print(f"CELL.GET_MOLECLIST. Labels or coordinates not found. Returning None")
            return None
        if len(self.labels) == 0 or len(self.coord) == 0:           
            if debug > 0: print(f"CELL.GET_MOLECLIST. Empty labels or coordinates. Returning None")
            return None
        if debug > 0: print(f"CELL.GET_MOLECLIST passed initial checks")

        ## If fragmented, then it reconstructs the unit cell
        if not hasattr(self,"fragmented"): self.check_fragmentation()
        if self.fragmented:
            print("CELL.GET_MOLECLIST. Reconstructing unit cell")
            self.reconstruct(debug=debug)
            print("CELL.GET_MOLECLIST. Arranging new coordinates from moleclist")
            new_coord = self.fix_cell_coord()
            if new_coord is None:
                print("CELL.GET_MOLECLIST. NONE received from FIX_CELL_COORD. Stopping")
                return 

        blocklist = split_species(self.labels, self.coord, cov_factor=cov_factor, debug=debug)
        
        if blocklist is None: 
            return None
        else :
            if debug > 0: print(f"CELL.GET_MOLECLIST: found {len(blocklist)} blocks")
            if debug > 0: print(f"CELL.GET_MOLECLIST: {blocklist=}")
        
        self.moleclist = []
        for b in blocklist:
            if debug > 0: print(f"CELL.GET_MOLECLIST: doing block={b}")
            mol_labels      = extract_from_list(b, self.labels, dimension=1)
            mol_coord       = extract_from_list(b, self.coord, dimension=1)
            mol_frac_coord  = extract_from_list(b, self.frac_coord, dimension=1)
            # Creates Molecule Object
            newmolec        = molecule(mol_labels, mol_coord, mol_frac_coord)
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
                if debug > 0: print(f"CELL.GET_MOLECLIST: splitting complex")
                newmolec.split_complex(debug=debug)
            # Not needed here, as the reconstruction will take care of it
            self.moleclist.append(newmolec)
        return self.moleclist

    ######
    def reconstruct(self, cov_factor: float=None, metal_factor: float=None, debug: int=0):
        from .Reconstruct import classify_fragments, fragments_reconstruct

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

    ######
    def reset_charge_assignment(self, debug: int=0):
        if not hasattr(self,"moleclist"): return None
        for mol in self.moleclist:
            mol.reset_charge()

    ######
    def view(self, size: str='default'):
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
            mode='text', text=symbols, textfont=dict(color='black', size=text_size), 
            hoverinfo='none', showlegend=False))

        for i, j in unique_bonds:
            fig.add_trace(go.Scatter3d(
                x=[positions[i, 0], positions[j, 0]], y=[positions[i, 1], positions[j, 1]],
                z=[positions[i, 2], positions[j, 2]], mode='lines',
                line=dict(color='gray', width=5), hoverinfo='none', showlegend=False))

        set_scene(fig, positions, width=width, height=height)
        fig.show()

    #########################################
    ### Functions to Interact with States ###
    #########################################
    def add_state(self, name: object, debug: int=0):
        from .Classes_State import state
        if not hasattr(self,"states"): setattr(self,"states",list([]))
        exists, new_state = self.find_state(name)
        if exists:  
            if debug > 0: print(f"CELL.ADD_STATE. State with same {name=} found, returning it")
            return new_state
        else:
            if debug > 0: print("CELL.ADD_STATE. Creating new state, returning it")
            new_state = state(self, name, debug=debug)
            self.states.append(new_state)
        return new_state

    ######
    def find_state(self, search_name: str, debug: int=0):
        from .Classes_State import state
        if not hasattr(self,"states"): setattr(self,"states",list([]))
        if debug > 0: print(f"CELL.FIND_STATE: Searching {search_name} in Cell object with {len(self.states)} states")
        for sta in self.states:
            if debug > 0: print(f"CELL.FIND_STATE: Comparing {search_name} with {sta.name}")
            if sta.name == search_name: 
                if debug > 0: print(f"CELL.FIND_STATE: state {search_name} found")
                return True, sta
        if debug > 0: print(f"CELL.FIND_STATE: state {search_name} not found")
        return False, None

######################
####    IMPORT    ####
######################
def import_cell(old_cell: object, debug: int=0) -> object:
    from .Classes_Specie  import import_molecule
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

    # If cell_param and cell_vector are None, the creation of the cell will fail
    new_cell            = cell(name, labels, coord, cell_vector, cell_param)
    new_cell.subtype    = "cell"
    new_cell.origin     = "import_cell"
    if debug > 0: print(f"IMPORT CELL: importing cell {new_cell}")

    ## Moleclist
    if debug > 0: print(f"IMPORT CELL: creating moleclist")
    moleclist = []
    for mol in old_cell.moleclist: 
        new_mol = import_molecule(mol, parent=new_cell, debug=debug)
        new_mol.set_bonds()
        new_mol.fix_ligands_rdkit_obj(debug=debug)
        moleclist.append(new_mol)
    new_cell.set_moleclist(moleclist)
    new_cell.fix_cell_coord()    ## In cell2mol, the cell object does not have the coordinates of the reconstructed cell. 
                                 ## However, the molecule and atom objects are updated (i.e. reconstructed). We use this info to update the cell

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
