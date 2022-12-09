import pickle
import sys
import os
import shutil
import csv

#scopepath = '/home/g4vela/SCOPE/Database_SCO/Scripts/Scope'
#sys.path.append(scopepath)

from Parse_General import search_string, read_lines_file 
from Read_Write import save_binary, load_binary
from Scope_Classes import crystal, sco_system

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
    for crystal in sorted(os.listdir(corepath)):
        if os.path.isdir(crystal):
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
def create_sco_system(name, cell2mol_path: str="default", debug: int=0) -> object:

    ### Creates a system. To do so, it...
    ### Reads all directories in a folder, searches for those containing both a Cell ".gmol" object and a ".cif" file.
    ### Each pair of ".gmol" and ".cif" will be considered a "crystal" and associates them to the "system" 
 
    # creates the system with basic info

    if cell2mol_path == "default": cell2mol_path = os.getcwd()
    newsystem = sco_system(name, cell2mol_path)
    # finds the crystals with Cell and cif file
    crystal_files, crystal_objects, crystal_paths, crystal_names, cif_paths  = find_crystals(name, cell2mol_path, debug=debug)
    # for each crystal:
    for idx, crys in enumerate(crystal_files):
        newcrystal = crystal(name, crystal_names[idx], crystal_paths[idx], crystal_objects[idx])
        newcrystal.read_cif_data(cif_paths[idx])
                            
        # At the end, it stores the data in the system object
        newsystem.list_of_crystals.append(newcrystal)

    # Determines Reference HS and LS molecules, if there are
    if debug > 0: print(len(newsystem.list_of_crystals), "crystals in system. Setting reference molecules.")
    newsystem.set_reference_molecs(debug=debug)
    return newsystem
##############################

