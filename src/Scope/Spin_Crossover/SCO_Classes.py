################################################
##### SYSTEM & CELL Classes Adapted to SCO #####
################################################
import os
from copy import deepcopy

from Scope.Classes_Cell         import *
from Scope.Classes_Cif          import *
from Scope.Classes_State        import *
from Scope.Classes_Specie       import *
from Scope.Classes_System       import system
from Scope.Read_Write           import load_binary, print_xyz

from Scope.Spin_Crossover.SCO_Structure import *

########################################
##### SYSTEM Object Adapted to SCO #####
########################################
class sco_system(system):
    def __init__(self, refcode: str) -> None:
        system.__init__(self, refcode)
        self.subtype              = "sco_system" 
        self.refcode              = refcode
        self.refcode_wo_digits    = ''.join([i for i in refcode if not i.isdigit()])

    def __repr__(self):
        to_print  = ''
        to_print += '-------------------------------------\n'
        to_print += '-- >>> SCOPE SCO-System Object >>> --\n'
        to_print += '-------------------------------------\n'
        to_print += system.__repr__(self, indirect=True)
        to_print += '\n'
        return to_print

    ######
    ## Load Cells from paths
    ######
    def load_single_cell2mol_folder(self, folder: str, overwrite: bool=False, debug: int=0):
        """
        Imports a sco_cell-class object from a cell2mol folder containing both a Cell-class ".gmol" object and a ".cif" file, 
        and integrates the data into the current SCO system.

        Args:
            path (str): Path to the folder containing cell2mol files.
            overwrite (bool, optional): If True, overwrites existing cell data. Defaults to False.
            debug (int, optional): Debug level for verbose output. Defaults to 0.

        Returns:
            self: The updated SCO system instance with the loaded cell and associated CIF data.

        Notes:
            - The method searches for files in the specified folder that match both the Cell ".gmol" and ".cif" formats.
            - Both files must be present for the cell to be loaded and integrated.
            - The Cell and CIF objects are linked to each other after loading.
            - If either file is missing or cannot be loaded, the method returns without modifying the system.
        """
        from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal
        if folder[-1] != '/': folder += '/'
        if not os.path.isdir(folder): return None
        cell_loaded  = False
        cif_loaded  = False

        # A first loop over the list of files, to ensure that all necessary data is there
        for fil in sorted(os.listdir(folder)):
            if fil.endswith(".gmol") and fil.startswith("Cell"):
                #try: 
                if debug > 0: print("Trying to load cell2mol CELL file from", folder+fil)
                new_cell        = import_cell(load_binary(folder+fil))        ## Imports to the generic Cell class
                if new_cell is None: continue
                new_cell        = convert_to_sco_cell(new_cell)               ## Converts to the SCO adapted Cell class
                if new_cell is None: continue
                cell_loaded     = True
                cell_path       = folder+fil
                if debug > 0: print("File loaded successfully")
                #except Exception as exc: 
                    #if debug > 0: print("Cell could not be loaded from", folder+fil)
                    #if debug > 0: print(exc)
                    #else: pass
            elif fil.endswith(".cif"):
                #try:
                if debug > 0: print("Trying to load CIF file from", folder+fil)
                new_cif         = cif(fil, folder+fil)
                cif_loaded      = True
                if debug > 0: print("File loaded successfully")
                #except Exception as exc: 
                    #if debug > 0: print("Cif file could not be loaded from", folder+fil)
                    #if debug > 0: print(exc)
                    #else: pass
        if debug > 0: print(f"IMPORT_cell2mol_FOLDER. Path: {folder}, {cell_loaded=} {cif_loaded=}")

        ## In SCO projects reading from cell2mol, I set that both the Cell and the cif have to be there.
        ## If the information could be properly retrieved, then it keeps the cell
        if cif_loaded and cell_loaded: 
            ## Prepare the cell object
            if debug > 0: print(f"IMPORT_cell2mol_FOLDER: Preparing the cell object in {cell_path}")
            new_cell.get_FeN6_molecules()
            new_cell.get_spin_and_phase_data()
            ## Link Cell and Cif objects
            new_cell.associate_cif(new_cif)
            new_cif.associate_cell(new_cell)
            ## Stores Path
            new_cell.set_path(cell_path)                               
            ## Create SCO System, and incorporate the cell
            if hasattr(new_cell,"name"): source_name = new_cell.name
            else:                        source_name = folder.split('/')[-2]
            self.add_source(source_name, new_cell, overwrite=overwrite, debug=debug)
        return self

    ###### Same for multiple folders
    def load_multiple_cell2mol_folders(self, path: str, overwrite: bool=False, debug: int=0):
        if path[-1] != '/': path += '/'
        if not os.path.isdir(path): return None
        for folder in sorted(os.listdir(path)):
            if os.path.isdir(path+folder):
                if debug > 0: print(f"Loading cell2mol folder {folder}")
                self.load_single_cell2mol_folder(path+folder, overwrite=overwrite, debug=debug)
        return self

    ######    
    def set_reference_molecs(self, overwrite: bool=False, debug: int=0):
        found_hs, hs = self.find_source("ref_hs_mol")
        found_ls, ls = self.find_source("ref_ls_mol")
        if not overwrite and found_hs and found_ls:
            if debug > 0: print(f"SET_REF_MOLEC: reference molecules already present. Use overwrite=True to reset them")
            return True

        ## Collect all FeN6 molecules in the cells of the system 
        pool = []
        for source in self.sources:
            if source.type == "cell": pool.extend(source.get_FeN6_molecules(overwrite=overwrite, debug=0))
        if debug > 0:      print(f"SET_REF_MOLEC: {len(pool)} FeN6 molecules in pool")
        if len(pool) == 0: print(f"SET_REF_MOLEC: Empty pool of reference molecules"); return False           

        ## Evaluates
        for mol in pool:
            mol.scope_FeNdist, mol.scope_FeNangle, mol.scope_epsylon = geom_sco_from_xyz(mol.labels, mol.coord, debug=0)
            mol.scope_guess_spin = guess_spin_state(int(mol.metals[0].charge), mol.scope_FeNdist[0], debug=0)

        ## Evaluates in pairs to make sure we get the 'same' molecule (same charge and chemical composition) for both spin states
        found_refs = False
        for idx, mol1 in enumerate(pool):
            for jdx, mol2 in enumerate(pool):
                if jdx != idx and not found_refs:
                    if mol1.scope_guess_spin == 'HS' and mol2.scope_guess_spin == 'LS' and mol1.charge == mol2.charge and mol1 == mol2:
                        found_refs = True
                        hs = deepcopy(mol1) #I'm making a deepcopy to avoid issues with references
                        ls = deepcopy(mol2)
        if not found_refs:
            hs = deepcopy(pool[0])
            ls = deepcopy(pool[0])

        ## Prepares the species:
        hs.name = "ref_hs_mol"
        hs.set_bonds()              
        hs.set_spin_config(4, typ='metals') 
        hs.fix_ligands_rdkit_obj()
        self.add_source(hs.name, hs)

        ls.name = "ref_ls_mol"
        ls.set_bonds()              
        ls.set_spin_config(0, typ='metals') # Not strictly necessary, but for clarity 
        ls.fix_ligands_rdkit_obj()
        self.add_source(ls.name, ls)

        # Creates Initial States:
        hs_ini_state = hs.add_state("initial")
        hs_ini_state.set_geometry(hs.labels, hs.coord)

        ls_ini_state = ls.add_state("initial")
        ls_ini_state.set_geometry(ls.labels, ls.coord)

        return True
        
    ######
    def set_reference_cells(self, overwrite: bool=False, debug: int=0):
        found_hs, hs = self.find_source("ref_hs_cell")
        found_ls, ls = self.find_source("ref_ls_cell")
        if not overwrite and found_hs and found_ls:
            if debug > 0: print(f"SET_REF_CELLS: reference HS or LS cells already present. Use overwrite=True to reset them")
            return True

        hs_ref_temp = 1000  ## Initial high value to find the lowest temperature
        ls_ref_temp = 1000  ## Initial high value to find the lowest temperature
        found_hs = False
        found_ls = False
        for source in self.sources:
            if source.type == "cell":
                cell = source
                if hasattr(cell,"warning_list"):        ## In cell2mol cells, this variable contains any warning found during the interpretation
                    if not any(cell.warning_list):      ## For a cell to be valid, no error is allowed
                        cell.get_spin_and_phase_data(debug=debug)

                    ## This try/except block is to skip cases in which the diffraction temperature is not known
                    if hasattr(cell,"cif"):
                        if hasattr(cell.cif,"diff_temp"):
                            try: cell.cif.diff_temp = float(cell.cif.diff_temp)
                            except: continue
                        else: continue
                    else: continue

                    ## Locate reference HS and LS crystals based on their diffraction temperature
                    if   cell.phase == "HS" and cell.cif.diff_temp < hs_ref_temp:
                        hs_ref_temp = cell.cif.diff_temp 
                        hs       = deepcopy(cell)
                        found_hs = True
                    elif cell.phase == "LS" and cell.cif.diff_temp < ls_ref_temp:
                        ls       = deepcopy(cell)
                        found_ls = True
                    elif cell.phase == "IS":
                        if debug > 0: print(f"SET_REF_CELLS: IS phase found but not implemented")

        if not found_hs:
            for source in self.sources:
                if source.type == "cell":
                    cell = source
                    if hasattr(cell,"warning_list"):
                        if not any(cell.warning_list):
                            hs = deepcopy(cell)
        if not found_ls:
            for source in self.sources:
                if source.type == "cell":
                    cell = source
                    if hasattr(cell,"warning_list"):
                        if not any(cell.warning_list):
                            ls = deepcopy(cell)

        if ls is None or hs is None: # Then something went wrong
            return False

        ## Prepares the cells:
        hs.name = "ref_hs_cell"
        hs.set_spin_config(4, typ='metals')
        for mol in hs.moleclist:
            mol.set_bonds()
            if mol.iscomplex: mol.fix_ligands_rdkit_obj()
        self.add_source(hs.name, hs)

        ls.name = "ref_ls_cell"
        ls.set_spin_config(0, typ='metals') # Not strictly necessary, but for clarity 
        for mol in ls.moleclist:
            mol.set_bonds()
            if mol.iscomplex: mol.fix_ligands_rdkit_obj()
        self.add_source(ls.name, ls)

        if hs.natoms != ls.natoms: print(f"Warning: different number of atoms in crystal; HS: {hs.natoms} vs. LS: {ls.natoms}")

        # Creates "initial" states:
        hs_ini_state = hs.add_state("initial")
        hs_ini_state.set_geometry(hs.labels, hs.coord)
        hs_ini_state.set_cell(hs.cell_vector, hs.cell_param)
        hs_ini_state.get_moleclist()
        
        ls_ini_state = ls.add_state("initial")
        ls_ini_state.set_geometry(ls.labels, ls.coord)
        ls_ini_state.set_cell(ls.cell_vector, ls.cell_param)
        ls_ini_state.get_moleclist()

        return True

######################################
##### CELL Object Adapted to SCO #####
######################################
class sco_cell(cell):
    def __init__(self, name: str, labels: list, pos: list, cell_vector: list=None, cell_param: list=None) -> None:
        cell.__init__(self, name, labels, pos, cell_vector, cell_param)
        self.subtype              = "sco_cell"

    def __repr__(self) -> None:
        to_print  = ''
        to_print += '-----------------------------------\n'
        to_print += '-- >>> SCOPE SCO-CELL Object >>> --\n'
        to_print += '-----------------------------------\n'
        to_print += cell.__repr__(self, indirect=True)
        if hasattr(self,"phase"):              to_print += f' Phase                 = {self.phase}\n'
        if hasattr(self,"hs_molar_fraction"):  to_print += f' HS Molar Fraction     = {self.hs_molar_fraction}\n'
        to_print += '\n'
        return to_print

    ######
    def get_FeN6_molecules(self, overwrite: bool=False, debug: int=0):
        if not hasattr(self,"moleclist"): self.get_moleclist()

        ## If already computed, returns them
        if not overwrite and hasattr(self,"FeN6s"):
            if len(self.FeN6s) > 0: return self.FeN6s

        self.FeN6s = []
        for mol in self.moleclist:
            if mol.iscomplex:
                keepit = False
                for met in mol.metals:
                    if hasattr(met,"charge"):
                        if debug > 0: print(f"    GET_FeN6s: Evaluating metal: {met}")
                        coord_sphere = met.get_coord_sphere_formula(debug=debug)
                        if debug > 0: print(f"    GET_FeN6s: Received {coord_sphere=}")
                        if coord_sphere == "N6": keepit = True
                    else: print(f"    GET_FeN6s: Metal does not have charge variable")
                if keepit:
                    if debug > 1: print(f"    GET_FeN6s: Entering geom_sco with:")
                    if debug > 1: print_xyz(mol.labels, mol.coord)
                    mol.scope_FeNdist, mol.scope_FeNangle, mol.scope_epsylon = geom_sco_from_xyz(mol.labels, mol.coord, debug=0)
                    mol.scope_guess_spin = guess_spin_state(int(mol.metals[0].charge), mol.scope_FeNdist[0], debug=debug)
                    self.FeN6s.append(mol)
                    if debug > 0: print(f"    GET_FeN6s: found {mol.scope_guess_spin} molecule with oxidation state: {mol.metals[0].charge=}") 
        return self.FeN6s 

    ######
    def get_spin_and_phase_data(self, debug: int=0):
        if not hasattr(self,"FeN6s"): self.get_FeN6_molecules()
        guess_spins = []
        self.phase = ''
        self.hs_molar_fraction = 0.0
        for mol in self.FeN6s:
            guess_spins.append(mol.scope_guess_spin)

        if all(sp == "HS" for sp in guess_spins):
            self.phase = "HS"
            self.hs_molar_fraction = 1.0
        elif all(sp == "LS" for sp in guess_spins):
            self.phase = "LS"
            self.hs_molar_fraction = 0.0
        elif "HS" in guess_spins and "LS" in guess_spins:
            self.phase = "IS"
            self.hs_molar_fraction = float(guess_spins.count("HS")/len(guess_spins))
        if debug > 0: print(f"Get_Spin&Phase: phase:{self.phase} and molar_frac: {self.hs_molar_fraction}")

#############################################
##### CONVERSION TO SCO-Adapted Classes #####
#############################################
def convert_to_sco_cell(generic_cell):
    if not isinstance(generic_cell, cell):
        print(f"CONVERT_TO_SCO_CELL: Input is not a 'cell' object")
        return None
    if isinstance(generic_cell, sco_cell):
        if debug > 0: print(f"CONVERT_TO_SCO_CELL: Input is already a 'sco_cell' object")
        return generic_cell
    new_cell = sco_cell(generic_cell.name, generic_cell.labels, generic_cell.coord, generic_cell.cell_vector, generic_cell.cell_param)
    for attr in generic_cell.__dict__.keys():
        if not hasattr(new_cell,attr):
            setattr(new_cell, attr, getattr(generic_cell,attr))
    return new_cell
