import pickle
import sys
import os
import shutil
import csv

from Scope.Parse_General import search_string, read_lines_file 
from Scope.Read_Write import save_binary, load_binary
from Scope.Scope_Classes import crystal, sco_system

###########################
##### Logistics Tools #####
###########################

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
    if debug > 0: print("Found", len(crystal_files),"in",cell2mol_path,":")
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
    if debug > 0: print(len(newsystem.list_of_crystals), "crystals in system. Setting reference molecules.")
    newsystem.set_reference_molecs(debug=debug)
    newsystem.set_reference_crystals(debug=debug)
    return newsystem
##############################

