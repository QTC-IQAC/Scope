import pickle
import sys
import copy
import os
import shutil
import csv
from copy import deepcopy

from Scope.Classes_Job import recipe, job
from Scope.Parse_Cif import get_cif_diffraction_data, get_cif_authors, get_cif_journal 
from Scope.Parse_General import search_string, read_lines_file 
from Scope.Geom_SCO_V1 import geom_sco_from_xyz, guess_spin_state
from Scope.Read_Write import save_binary, load_binary

from cell2mol.tmcharge_common import labels2formula
from cell2mol.elementdata import ElementData
elemdatabase = ElementData()

############################
##### SYSTEM & CRYSTAL #####
############################
class sco_system(object):
    def __init__(self, refcode: str, cell2mol_path: str, scope_path: str) -> None:
        self.type = "sco_system" 
        self.version = "0.2" 
        self.refcode = refcode
        self.refcode_wo_digits = ''.join([i for i in refcode if not i.isdigit()])
        self.cell2mol_path = cell2mol_path
        self.scope_path = scope_path
        self.list_of_crystals = []
        self.list_of_recipes = []
    
    ##########
    def set_reference_molecs(self, debug: int=0):
        self.hasLS = False
        self.hasHS = False

        pool = []
        ##### All these below could go somewhere else #####
        ### Colects reference molecules from all crystals
        for crys in self.list_of_crystals:
            if len(crys.list_of_molecules) == 0: crys.get_FeN6_molecules(debug=debug)
            for mol in crys.list_of_molecules:
                pool.append(mol) 
            
        ### Selects reference molecules
        if debug > 0: print(f"{len(pool)} tmcs in pool")
        if len(pool) > 0:
            for idx, mol in enumerate(pool):
                if mol.scope_guess_spin == 'HS' and not self.hasHS:
                    self.hasHS = True
                    self.HSid = idx
                    self.HSref = deepcopy(mol)
                    if debug > 0: print(f"HS reference molecule found")
                elif mol.scope_guess_spin == 'LS' and not self.hasLS:
                    self.hasLS = True
                    self.LSid = idx
                    self.LSref = deepcopy(mol)
                    if debug > 0: print(f"LS reference molecule found")
            ### If it hasn't found any molecule that can be classified as HS and LS... then takes anything
            if not self.hasHS:
                self.HSid = 0 
                self.HSref = deepcopy(pool[0])
                self.HSref.scope_guess_spin == 'HS'
                if debug > 0: print(f"HS reference molecule assumed")
            if not self.hasLS:
                self.LSid = 0 
                self.LSref = deepcopy(pool[0])
                self.LSref.scope_guess_spin == 'LS'
                if debug > 0: print(f"LS reference molecule assumed")
        else: print("Empty pool of reference molecules")

    ##########
    def set_reference_crystals(self, debug: int=0):
        self.has_HS_ref_crys = False
        self.has_LS_ref_crys = False

        #### Selects reference crystals:
        HS_ref_crys_temp = 1000
        LS_ref_crys_temp = 1000
        for idx, crys in enumerate(self.list_of_crystals):
  
            if not hasattr(crys,"phase"): crys.get_spin_and_phase_data(debug=debug)

            ### This try/except block is to skip cases in which the diffraction temperature is not known
            try: crys.diff_temp = float(crys.diff_temp)
            except: continue

            if crys.phase == "HS" and crys.diff_temp < HS_ref_crys_temp: 
                self.HS_ref_crys = deepcopy(crys)
                self.HS_ref_crys_id = idx
                self.HS_ref_crys_temp = crys.diff_temp                
                self.has_HS_ref_crys = True
                if debug > 0: print(f"HS reference crystal found")
            elif crys.phase == "LS" and crys.diff_temp < LS_ref_crys_temp: 
                self.LS_ref_crys = deepcopy(crys)
                self.LS_ref_crys_id = idx
                self.LS_ref_crys_temp = crys.diff_temp                
                self.has_LS_ref_crys = True
                if debug > 0: print(f"LS reference crystal found")
            elif crys.phase == "IS":
                if debug > 0: print(f"    Set_REF_Crystals: IS phase found but not implemented")

        if not self.has_HS_ref_crys:
            self.HS_ref_crys = deepcopy(self.list_of_crystals[0])
            self.HS_ref_crys_id = 0
            self.HS_ref_crys.phase == 'HS'
            if debug > 0: print(f"HS reference crystal assumed")
        if not self.has_LS_ref_crys:
            self.LS_ref_crys = deepcopy(self.list_of_crystals[0])
            self.LS_ref_crys_id = 0
            self.LS_ref_crys.phase == 'LS'
            if debug > 0: print(f"LS reference crystal assumed")

    ##########
    def get_recipe(self, obj, keyword, debug: int=0):
        ## obj is the actual object that will be computed. For instance, a molecule, or the cell.
        recipe_exists = False
        for idx, rec in enumerate(self.list_of_recipes):
            if rec.keyword == keyword and rec.object == obj and os.path.isdir(rec.path):
                recipe_exists = True
                current_recipe = rec

        ## If such entry doesn't exist. Then we create it
        if not recipe_exists:
            current_recipe = recipe(obj, self.scope_path+'/'+keyword, keyword)
            self.list_of_recipes.append(current_recipe)
            if not os.path.isdir(self.scope_path+'/'+keyword): os.makedirs(self.scope_path+'/'+keyword)
        return recipe_exists, current_recipe

class crystal(object):
    def __init__(self, refcode: str, name: str, cell2mol_path: str, cell: object) -> None:
        self.type = "crystal" 
        self.version = "0.2" 
        self.refcode = refcode
        self.name = name
        self.cell = cell 
        self.cell2mol_path = cell2mol_path
        self.list_of_molecules = []

    def read_cif_data(self, cifpath: str) -> None:
        self.diff_temp = get_cif_diffraction_data(cifpath) 
        self.authors = get_cif_authors(cifpath)
        self.journal_year =   get_cif_journal(cifpath)[0]
        self.journal_name =   get_cif_journal(cifpath)[1] 
        self.journal_volume = get_cif_journal(cifpath)[2]
        self.journal_page =   get_cif_journal(cifpath)[3]

    def add_iso_calcs_path(self, iso_calcs_path: str='Unk'):
        self.iso_calcs_path = iso_calcs_path

    def get_spin_and_phase_data(self, debug: int=0):
        self.guess_spins = [] 
        self.phase = ''
        self.HSmolarfraction = 0.0
        for idx, mol in enumerate(self.list_of_molecules):
            self.guess_spins.append(mol.scope_guess_spin)

        if all(sp == "HS" for sp in self.guess_spins): 
            self.phase = "HS"
            self.HSmolarfraction = 1.0
        elif all(sp == "LS" for sp in self.guess_spins): 
            self.phase = "LS"
            self.HSmolarfraction = 0.0
        elif "HS" in self.guess_spins and "LS" in self.guess_spins: 
            self.phase = "IS"
            self.HSmolarfraction = float(self.guess_spins.count("HS")/len(self.guess_spins))
        if debug > 0: print(f"    Get_Spin&Phase: phase:{self.phase} and molar_frac: {self.HSmolarfraction}")

    def get_FeN6_molecules(self, debug: int=0):
        for idx, mol in enumerate(self.cell.moleclist):
            if mol.type == "Complex":
                keepit = False
                for met in mol.metalist:
                    if hasattr(met, "coord_sphere"):
                        formula = labels2formula(met.coord_sphere)
                        if formula == "N6": keepit = True
                    else: print("No coord_sphere variable in metal object")
                if keepit:
                    ox_state = mol.metalist[0].totcharge
                    mol.scope_FeNdist, mol.scope_FeNangle = geom_sco_from_xyz(mol.labels, mol.coord)
                    mol.scope_guess_spin = guess_spin_state(int(ox_state), mol.scope_FeNdist[0])
                    self.list_of_molecules.append(mol)
                    if debug > 0: print(f"    GET_FeN6: found {mol.scope_guess_spin} molecule of OS: {ox_state}") 
       
################################
##### ASSOCIATED FUNCTIONS #####
################################
def find_crystals(name: str, corepath: str, debug: int=0):

    ### Reads all directories in a folder, searches for those containing both a Cell ".gmol" object and a ".cif" file.
    name_wo = ''.join([i for i in name if not i.isdigit()])

    crystal_filenames = []
    crystal_objects = []
    paths_list = []
    names_list = []
    cif_paths = []
    if debug > 0: print("FIND: Corepath:", corepath)
    for crystal in sorted(os.listdir(corepath)):
        if os.path.isdir(corepath+"/"+crystal):
            crystal_path = corepath+"/"+crystal
            cell_found = False
            cell_loaded = False
            cif_found = False

            # A first loop over the list of files, to ensure that the data is there
            for fil in sorted(os.listdir(crystal_path)):
                if fil.endswith(".gmol") and fil.startswith("Cell") and name_wo in fil:
                    cell_found = True
                    try: 
                        if debug > 0: print("Trying to load cell from", crystal_path+"/"+fil)
                        cell = load_binary(crystal_path+"/"+fil)
                        cell_loaded = True
                    except: 
                        if debug > 0: print("Cell could not be loaded from", crystal_path+"/"+fil)
                        else: pass
 
                if fil.endswith(".cif") and name_wo in fil:
                    cif_found = True
            if debug > 0: print("FIND: crystal_path", crystal_path, "results:", cell_found, cell_loaded, cif_found)

            ## If the information could be properly retrieved, then it stores it. 
            if cell_found and cell_loaded and cif_found: 
                for fil in sorted(os.listdir(crystal_path)):
                    if fil.endswith(".gmol") and fil.startswith("Cell") and name_wo in fil:
                        crystal_filenames.append(fil)
                        crystal_objects.append(load_binary(crystal_path+"/"+fil))
                        paths_list.append(crystal_path)
                        if crystal == "Original": crystal = name
                        names_list.append(crystal)
                    if fil.endswith(".cif") and name_wo in fil:
                        cif_paths.append(crystal_path+"/"+fil)

    return crystal_filenames, crystal_objects, paths_list, names_list, cif_paths

##############################
def create_sco_system(name, cell2mol_path: str="default", scope_path: str="default", debug: int=0) -> object:

    ### Creates a system. To do so, it...
    ### Reads all directories in a folder, searches for those containing both a Cell ".gmol" object and a ".cif" file.
    ### Each pair of ".gmol" and ".cif" will be considered a "crystal" and associates them to the "system" 
 
    # creates the system with basic info
    if cell2mol_path == "default": cell2mol_path = os.getcwd()
    if scope_path == "default": scope_path = os.getcwd()

    # finds the crystals with Cell and cif file
    newsystem = sco_system(name, cell2mol_path, scope_path)
    if debug > 0: print("Searching crystals in", cell2mol_path)
    crystal_files, crystal_objects, crystal_paths, crystal_names, cif_paths  = find_crystals(name, cell2mol_path, debug=debug)
    if debug > 0: print("Found", len(crystal_files),"crystals in",cell2mol_path,":")
    # for each crystal:
    for idx, crys in enumerate(crystal_files):
        if debug > 0: print(idx, name, crystal_names[idx], crystal_paths[idx])
        newcrystal = crystal(name, crystal_names[idx], crystal_paths[idx], crystal_objects[idx])
        newcrystal.read_cif_data(cif_paths[idx])
        newcrystal.get_FeN6_molecules(debug=debug)
        newcrystal.get_spin_and_phase_data(debug=debug)
                            
        # At the end, it stores the data in the system object
        newsystem.list_of_crystals.append(newcrystal)

    # Determines Reference HS and LS molecules and crystals, if there are
    if debug > 0: print(len(newsystem.list_of_crystals), "crystals in system. Setting reference molecules and crystals")
    newsystem.set_reference_molecs(debug=debug)
    newsystem.set_reference_crystals(debug=debug)
    return newsystem
##############################

