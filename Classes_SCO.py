import pickle
import sys
import copy
import os
import shutil
import csv
from copy import deepcopy

from Scope.Parse_Cif import get_cif_diffraction_data, get_cif_authors, get_cif_journal 
from Scope.Parse_General import search_string, read_lines_file 
from Scope.Structure_SCO import geom_sco_from_xyz, guess_spin_state
from Scope.Read_Write import save_binary, load_binary, print_xyz
from Scope.Classes_State import state
from Scope.Classes_Molecule import *
#from Scope.Bibliography import *

from Scope.Workflow import Branch
from Scope.Workflow.Branch import *

from Scope.Adapted_from_cell2mol import labels2formula
from Scope.Elementdata import ElementData
elemdatabase = ElementData()

###########################
def cell_postocoord(cell: object) -> None:
    if not hasattr(cell,"coord") and hasattr(cell,"pos"): setattr(cell,"coord",cell.pos)

############################
##### SYSTEM & CRYSTAL #####
############################
class sco_system(object):
    def __init__(self, refcode: str, cell2mol_path: str, calcs_path: str, sys_path: str) -> None:
        self.type                 = "sco_system" 
        self.version              = "V3" 
        self.refcode              = refcode
        self.refcode_wo_digits    = ''.join([i for i in refcode if not i.isdigit()])
        self.cell2mol_path        = cell2mol_path
        self.calcs_path           = calcs_path
        self.sys_path             = sys_path
        self.crystals             = []
        self.branches             = []
        self.results              = dict()
        self.status               = "active"

    def save(self, filepath: str=None):
        if filepath is None: filepath = self.sys_path+self.refcode+'.sys'
        from Scope.Read_Write import save_binary
        save_binary(self, filepath)

###################################
### Functions to restart system ###
###################################
    # Eventually, reset_paths could be "change_cluster" function
    def reset_paths(self, environment: object, debug: int=0) -> None: 
        reset = False

        ## Environment Sends the Global Paths. The Paths Specifics to the System are obtained here: 
        target_sys_path             = environment.sys_path+self.refcode+'/'
        target_calcs_path           = environment.calcs_path+self.refcode+'/'
        target_cell2mol_path        = environment.cell2mol_path+self.refcode+'/'

        ## We make sure that those exist:
        if os.path.isdir(target_sys_path) and os.path.isdir(target_calcs_path) and os.path.isdir(target_cell2mol_path):      
            reset = True
            self.sys_path             = target_sys_path
            self.calcs_path           = target_calcs_path 
            self.cell2mol_path        = target_cell2mol_path 
            if debug > 0: print(f"RESET_PATHS: new system paths: {self.cell2mol_path}, {self.calcs_path}, {self.sys_path}")

            for crys in self.crystals:
                if debug > 0: print(f"RESET_PATHS: new cell2mol path for crystal {crys.name}")
                crys.reset_paths(target_cell2mol_path)

            ## We go down the hierarchy to change branch, recipe, jobs, and calculation paths:
            for br in self.branches:
                tmp = target_calcs_path+br.keyword+'/'
                if os.path.isdir(tmp): br.path = tmp 
                else: print(f"RESET_PATHS: {tmp} path does not exist")
                for rec in br.recipes:
                    rec.path = br.path
                    for job in rec.jobs:
                        job.path = rec.path
                        for comp in job.computations:
                            if job.setup == "findiff" and os.path.isdir(job.path+"findiff"): comp.path = job.path+"findiff_test2/"
                            else:                                                            comp.path = job.path
                            comp.inp_path = comp.path+comp.inp_name
                            comp.out_path = comp.path+comp.out_name
                            comp.sub_path = comp.path+comp.sub_name
                            comp.check_files()
                            if debug > 0: print(f"RESET_PATHS: new computation path: {comp.inp_path}")
        return reset

    ################################
    def remove_all_branches(self) -> None:
        if hasattr(self,"branches"): delattr(self,"branches"); setattr(self,"branches",[])
        return self.branches
       
    ## to make it work, target must be an attribute of branch, and recipes created at branch level, not at system level as in def add_branch below
    ###def reset_branch(self, keyword: str) -> None:
    ###    exists, old_branch = self.find_branch(keyword)
    ###    if exists:
    ###        target = old_branch.target
    ###    self.remove_branch(self, br_keyword=keyword):
    ###    new_branch = self.add_branch(keyword, target: str, debug: int=0):
    ###    return new_branch

    ################################
    def reset_crystals(self) -> None:
        if hasattr(self,"crystals"): delattr(self,"crystals"); setattr(self,"crystals",[])

    ################################
    def find_branch(self, keyword: str, debug: int=0):
        if debug > 1: print("finding branch with keyword:", keyword)
        if debug > 1: print("there are", len(self.branches), "branches in system")
        if len(self.branches) == 0: return False, None
        for idx, br in enumerate(self.branches):
            if debug > 1: print("evaluating branch with keyword:", br.keyword, "and path:", br.path)
            if br.keyword.lower() == keyword.lower():
                if not os.path.isdir(br.path): 
                    print(f"WARNING: The path associated with the computations of this branch (below) does not exist. Loading the branch anyway")
                    print(f"WARNING: {br.path=}")
                return True, br
        return False, None

    ################################
    def find_computation(self, branch_keyword: str, recipe_keyword: str, job_keyword: str, comp_keyword: str='', comp_step: int=1, comp_run_number: int=1, debug: int=0):
        if len(self.branches) == 0: return False, None
        for br in self.branches:
            if br.keyword.lower() == branch_keyword.lower():
                if len(br.recipes) == 0: return False, None
                for rec in br.recipes:
                    if rec.subject.spin.lower() == recipe_keyword.lower():
                        if len(rec.jobs) == 0: return False, None
                        for job in rec.jobs:
                            if job.keyword.lower() == job_keyword.lower():
                                if len(job.computations) == 0: return False, None
                                for idx, comp in enumerate(job.computations):
                                    if comp.run_number == comp_run_number and comp.step == comp_step and comp.keyword.lower() == comp_keyword.lower(): return True, comp
        return False, None

    ################################
    def add_branch(self, keyword: str, debug: int=0):
        new_branch = branch(self.calcs_path+keyword, keyword, self, debug=debug)
        if not os.path.isdir(self.calcs_path+keyword): 
            try: os.makedirs(self.calcs_path+keyword)
            except Exception as exc:
                 print(f"Error creating branch folder in {self.calcs_path+keyword}")
                 print(exc)

#        if   target.lower()  == "ref_mol"  and hasattr(self,"HS_ref_mol")  and hasattr(self,"LS_ref_mol"):   object_list = list([self.HS_ref_mol,self.LS_ref_mol])
#        elif target.lower()  == "ref_crys" and hasattr(self,"HS_ref_crys") and hasattr(self,"LS_ref_crys"):  object_list = list([self.HS_ref_crys.cell,self.LS_ref_crys.cell])
#        else: print("Get_Branch: target could not be identified or not available. Recipes were not created")
#
#        ## Creates recipes for the branch. One for each object. 
#        for idx, gmol in enumerate(object_list):
#            new_recipe = new_branch.add_recipe(gmol)
        self.branches.append(new_branch)
        return new_branch

    ################################
    def remove_branch(self, br_keyword=None):
        found = False
        for idx, br in enumerate(self.branches):
            if br.keyword == str(br_keyword): found = True; found_idx = idx
        if found:
            to_delete = self.branches[found_idx]
            del self.branches[found_idx]

    ################################
    def set_reference_molecs(self, debug: int=0):
 
        if debug > 0: print("Setting Reference Molecules")
        self.hasLS = False
        self.hasHS = False

        pool = []
        ##### All these below could go somewhere else #####
        ### Colects reference molecules from all crystals
        print(len(self.crystals))
        for crys in self.crystals:
            if len(crys.list_of_molecules) == 0: 
                if debug > 0: print("Setting Reference Molecules; Collecting Molecules with FeN6 coordination Sphere") 
                crys.get_FeN6_molecules(debug=debug)
            for mol in crys.list_of_molecules:
                if not mol.check_fragmentation(): pool.append(mol)

        ### Selects reference molecules
        if debug > 0: print(f"{len(pool)} tmcs in pool")
        if len(pool) > 0:
            for idx, mol in enumerate(pool):
                if mol.scope_guess_spin == 'HS' and not self.hasHS:
                    self.hasHS = True
                    self.HS_ref_mol_id = idx
                    self.HS_ref_mol = deepcopy(mol)
                    #self.HS_ref_mol = import_molecule(mol, parent=self)
                    self.HS_ref_mol.spin = 'HS'
                    self.HS_ref_mol._sys = self  #Duplicate with parent to avoid crashes
                    if debug > 0: print(f"HS reference molecule found")
                elif mol.scope_guess_spin == 'LS' and not self.hasLS:
                    self.hasLS = True
                    self.LS_ref_mol_id = idx
                    self.LS_ref_mol = deepcopy(mol)
                    #self.LS_ref_mol = import_molecule(mol, parent=self)
                    self.LS_ref_mol.spin = 'LS'
                    self.LS_ref_mol._sys = self   #Duplicate with parent to avoid crashes
                    if debug > 0: print(f"LS reference molecule found")
            ### If it hasn't found any molecule that can be classified as HS and LS... then takes anything
            if not self.hasHS:
                self.HS_ref_mol_id = 0
                #self.HS_ref_mol = import_molecule(pool[0], parent=self)
                self.HS_ref_mol = deepcopy(pool[0])
                self.HS_ref_mol.spin = 'HS'
                self.HS_ref_mol._sys = self  #Duplicate to avoid crashes
                if debug > 0: print(f"HS reference molecule assumed")
            if not self.hasLS:
                self.LS_ref_mol_id = 0
                self.LS_ref_mol = deepcopy(pool[0])
                #self.LS_ref_mol = import_molecule(pool[0], parent=self)
                self.LS_ref_mol.spin = 'LS'
                self.LS_ref_mol._sys = self   #Duplicate to avoid crashes
                if debug > 0: print(f"LS reference molecule assumed")

        else: print("Empty pool of reference molecules"); return False

        assert self.HS_ref_mol.natoms != self.LS_ref_mol.natoms; f"Warning: different number of atoms in molecule; HS: {self.HS_ref_mol.natoms} vs. LS: {self.LS_ref_mol.natoms}"

        ## Fixes RDKIT objects from cell2mol
        self.HS_ref_mol.set_bonds()
        self.HS_ref_mol.fix_ligands_rdkit_obj()
        self.LS_ref_mol.set_bonds()
        self.LS_ref_mol.fix_ligands_rdkit_obj()

        ## Creates "initial" states:
        if hasattr(self,"HS_ref_mol"):
            HS_ini_state = state(self.HS_ref_mol, "initial")
            HS_ini_state.set_geometry(self.HS_ref_mol.labels, self.HS_ref_mol.coord)
        if hasattr(self,"LS_ref_mol"):
            LS_ini_state = state(self.LS_ref_mol, "initial")
            LS_ini_state.set_geometry(self.LS_ref_mol.labels, self.LS_ref_mol.coord)
        return True
        
    ################################
    def set_reference_crystals(self, debug: int=0):
        self.has_HS_ref_crys = False
        self.has_LS_ref_crys = False

        #### Selects reference crystals:
        HS_ref_crys_temp = 1000
        LS_ref_crys_temp = 1000
        for idx, crys in enumerate(self.crystals):
  
            #if not hasattr(crys.cell,"warning_list"): print("NO WARNING LIST")
            if hasattr(crys.cell,"warning_list"):
                if not any(crys.cell.warning_list):
                    if not hasattr(crys,"phase"): crys.get_spin_and_phase_data(debug=debug)

                    ### This try/except block is to skip cases in which the diffraction temperature is not known
                    try: crys.diff_temp = float(crys.diff_temp)
                    except: continue

                    if crys.phase == "HS" and crys.diff_temp < HS_ref_crys_temp: 
                        self.HS_ref_crys = deepcopy(crys)
                        self.HS_ref_crys_id = idx
                        self.HS_ref_crys_temp = crys.diff_temp                
                        HS_ref_crys_temp = crys.diff_temp 
                        self.HS_ref_crys._sys  == self 
                        self.HS_ref_crys.type  == "crystal"
                        setattr(self.HS_ref_crys.cell,"_sys",self)
                        setattr(self.HS_ref_crys.cell,"type","cell")
                        setattr(self.HS_ref_crys.cell,"spin","HS")
                        self.has_HS_ref_crys = True
                        if debug > 0: print(f"HS reference crystal found, new temperature is:", crys.diff_temp)
                    elif crys.phase == "LS" and crys.diff_temp < LS_ref_crys_temp: 
                        self.LS_ref_crys = deepcopy(crys)
                        self.LS_ref_crys_id = idx
                        self.LS_ref_crys_temp = crys.diff_temp                
                        LS_ref_crys_temp = crys.diff_temp                
                        self.LS_ref_crys._sys  == self 
                        self.LS_ref_crys.type  == "crystal"
                        setattr(self.LS_ref_crys.cell,"_sys",self)
                        setattr(self.LS_ref_crys.cell,"spin","LS")
                        setattr(self.LS_ref_crys.cell,"type","cell")
                        self.has_LS_ref_crys = True
                        if debug > 0: print(f"LS reference crystal found, new temperature is:", crys.diff_temp)
                    elif crys.phase == "IS":
                        if debug > 0: print(f"    Set_REF_Crystals: IS phase found but not implemented")

        if not self.has_HS_ref_crys:
            for idx, crys in enumerate(self.crystals):
                if hasattr(crys.cell,"warning_list"):
                    if not any(crys.cell.warning_list):
                        self.HS_ref_crys = deepcopy(crys)
                        self.HS_ref_crys_id = 0
                        self.HS_ref_crys.phase == 'HS'
                        self.HS_ref_crys._sys  == self 
                        setattr(self.HS_ref_crys.cell,"type","cell")
                        setattr(self.HS_ref_crys.cell,"_sys",self)
                        setattr(self.HS_ref_crys.cell,"spin","HS")
                        if debug > 0: print(f"HS reference crystal assumed")

        if not self.has_LS_ref_crys:
            for idx, crys in enumerate(self.crystals):
                if hasattr(crys.cell,"warning_list"):
                    if not any(crys.cell.warning_list):
                        self.LS_ref_crys = deepcopy(crys)
                        #self.LS_ref_crys = deepcopy(self.crystals[0])
                        self.LS_ref_crys_id = 0
                        self.LS_ref_crys.phase == 'LS'
                        self.LS_ref_crys._sys  == self 
                        setattr(self.LS_ref_crys.cell,"type","cell")
                        setattr(self.LS_ref_crys.cell,"_sys",self)
                        setattr(self.LS_ref_crys.cell,"spin","LS")
                        if debug > 0: print(f"LS reference crystal assumed")

        # Creates "initial" states:
        if hasattr(self,"HS_ref_crys") and hasattr(self,"LS_ref_crys"): 
            if self.HS_ref_crys.cell.natoms != self.LS_ref_crys.cell.natoms: 
                print(f"Warning: different number of atoms in crystal; HS: {self.HS_ref_crys.cell.natoms} vs. LS: {self.LS_ref_crys.cell.natoms}")

        if hasattr(self,"HS_ref_crys"):
            if debug > 0: print(f"SCO.SET_REF_CRYS: creating initial_state for HS crystal")
            HS_ini_state = state(self.HS_ref_crys.cell, "initial")
            HS_ini_state.set_geometry(self.HS_ref_crys.cell.labels, self.HS_ref_crys.cell.coord)
            HS_ini_state.set_cell(self.HS_ref_crys.cell.cell_vector, self.HS_ref_crys.cell.cell_param)
            HS_ini_state.get_moleclist()
        
        if hasattr(self,"LS_ref_crys"):
            if debug > 0: print(f"SCO.SET_REF_CRYS: creating initial_state for LS crystal")
            LS_ini_state = state(self.LS_ref_crys.cell, "initial")
            LS_ini_state.set_geometry(self.LS_ref_crys.cell.labels, self.LS_ref_crys.cell.coord)
            LS_ini_state.set_cell(self.LS_ref_crys.cell.cell_vector, self.LS_ref_crys.cell.cell_param)
            LS_ini_state.get_moleclist()
        return True

    ################################
    def __repr__(self) -> None:
        to_print  = f'---------------------------------------------------\n'
        to_print += f' Formatted input interpretation of SCO System()\n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Refcode               = {self.refcode}\n'
        to_print += f' Cell Path             = {self.cell2mol_path}\n'
        to_print += f' Calculations Path     = {self.calcs_path}\n'
        to_print += f' System Path           = {self.sys_path}\n'
        to_print += f' #Crystals             = {len(self.crystals)}\n'
        to_print += f' #Branches             = {len(self.branches)}\n'
        to_print += f'---------------------------------------------------\n'
        return to_print

###########################
###### Class Crystal ######
###########################
class crystal(object):
    def __init__(self, refcode: str, name: str, cell2mol_path: str, cell_name: str, sys: object) -> None:
        self.type              = "crystal" 
        self.version           = "V3" 
        self.refcode           = refcode
        self.name              = name
        self.cell2mol_path     = cell2mol_path
        self.cell_name         = cell_name
        self.cell_path         = cell2mol_path+cell_name
        self.cell              = import_cell(load_binary(self.cell_path))
        self.list_of_molecules = []
        self._sys              = sys

        self.fix_cell_coord()

    ################################
    def reset_paths(self, new_path) -> None:
        if os.path.isdir(new_path): self.cell2mol_path = new_path
        if os.path.isfile(new_path+self.cell_name): self.cell_path = new_path+self.cell_name

    ################################
    def reload_cell(self) -> None:
        self.cell              = import_cell(load_binary(self.cell_path))

    ################################
    ## In cell2mol, the cell object does not have the coordinates of the reconstructed cell. 
    ## However, the molecule and atom objects are updated (i.e. reconstructed). We use this info to update the cell
    def fix_cell_coord(self, debug: int=0) -> None:
        self.cell.labels = []
        self.cell.pos    = []
        self.cell.coord  = []
        indices          = []
        for mol in self.cell.moleclist:
            #if not hasattr(mol,"atoms"): 
            #    mol.set_atoms(create_adjacencies=True, debug=debug) 
            #else:
            #    atoms = []
            #    for at in mol.atoms:
            #        new_atom = import_atom(at, parent=mol)
            #        atoms.append(new_atom)
            #    mol.set_atoms(atomlist=atoms)
            for idx, a in enumerate(mol.atoms):
                self.cell.labels.append(a.label)
                self.cell.pos.append(a.coord)
                self.cell.coord.append(a.coord)
                if hasattr(a,"parent_index"): indices.append(a.parent_index)
                elif hasattr(a,"index"):      indices.append(a.index)
                else:                         indices.append(idx)

        ## Below is to order the atoms as in the original cell, using the indices stored in the molecule object
        self.cell.labels = [x for _, x in sorted(zip(indices, self.cell.labels), key=lambda pair: pair[0])]
        self.cell.pos    = [x for _, x in sorted(zip(indices, self.cell.pos), key=lambda pair: pair[0])]
        self.cell.coord  = [x for _, x in sorted(zip(indices, self.cell.coord), key=lambda pair: pair[0])]
        assert len(self.cell.labels) == len(self.cell.pos)

    ################################
    def read_cif_data(self, cifpath: str) -> None:
        self.diff_temp         = get_cif_diffraction_data(cifpath) 
        self.authors           = get_cif_authors(cifpath)
        self.journal_year      = get_cif_journal(cifpath)[0]
        self.journal_name      = get_cif_journal(cifpath)[1] 
        self.journal_volume    = get_cif_journal(cifpath)[2]
        self.journal_page      = get_cif_journal(cifpath)[3]

    ################################
    def get_FeN6_molecules(self, debug: int=0):
        for idx, mol in enumerate(self.cell.moleclist):
            if mol.iscomplex:
                keepit = False
                for met in mol.metals:
                    if hasattr(met,"charge"):
                        if debug > 0: print(f"    GET_FeN6: Evaluating metal: {met}")
                        coord_sphere = met.get_coord_sphere_formula(debug=debug)
                        if debug > 0: print(f"    GET_FeN6: Received {coord_sphere=}")
                        if coord_sphere == "N6": keepit = True
                    else:
                        print(f"    GET_FeN6: Metal does not have charge variable")
                if keepit:
                    ox_state = mol.metals[0].charge
                    if debug > 1: print(f"    GET_FeN6: Entering geom_sco with:")
                    if debug > 1: print_xyz(mol.labels, mol.coord)
                    mol.scope_FeNdist, mol.scope_FeNangle, mol.scope_epsylon = geom_sco_from_xyz(mol.labels, mol.coord, debug=0)
                    mol.scope_guess_spin = guess_spin_state(int(ox_state), mol.scope_FeNdist[0], debug=debug)
                    self.list_of_molecules.append(mol)
                    if debug > 0: print(f"    GET_FeN6: found {mol.scope_guess_spin} molecule with {ox_state=}") 
        return self.list_of_molecules

    ################################
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

    ####### Bibliography options:
    def get_search(self, verbose: bool=False, download: bool=True, debug: int=0):
        self.search = get_search(self, verbose=verbose, download=download, debug=debug) 
        return self.search

    def get_abstract(self, debug: int=0):
        if not hasattr(self,"search"): self.get_search(debug=debug)
        self.abstract = get_abstract(self.search, debug=debug)
        return self.abstract

    def get_title(self, debug: int=0):
        if not hasattr(self,"search"): self.get_search(debug=debug)
        self.title = get_title(self.search, debug=debug) 
        return self.title

    def get_doi(self, debug: int=0):
        if not hasattr(self,"search"): self.get_search(debug=debug)
        self.doi = get_doi(self.search, debug=debug) 
        return self.doi
    ####### 

    def __repr__(self) -> None:
        to_print  = f'---------------------------------------------------\n'
        to_print += f' Formatted input interpretation of Crystal()\n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Refcode               = {self.refcode}\n'
        to_print += f' Cell2mol Path         = {self.cell2mol_path}\n'
        if hasattr(self,"diff_temp"):        to_print += f' Diffraction Temp      = {self.diff_temp}\n'         
        if hasattr(self,"authors"):          to_print += f' Authors               = {self.authors}\n'           
        if hasattr(self,"journal_year"):     to_print += f' Year of Publication   = {self.journal_year}\n'      
        if hasattr(self,"journal_name"):     to_print += f' Journal Name          = {self.journal_name}\n'      
        if hasattr(self,"journal_volume"):   to_print += f' Journal Volume        = {self.journal_volume}\n'    
        if hasattr(self,"journal_page"):     to_print += f' Journal Page          = {self.journal_page}\n'      
        if hasattr(self,"phase"):            to_print += f' Phase                 = {self.phase}\n'
        if hasattr(self,"HSmolarfraction"):  to_print += f' HS Molar Fraction     = {self.HSmolarfraction}\n'
        return to_print
       
################################
##### ASSOCIATED FUNCTIONS #####
################################
def find_crystals(name: str, corepath: str, debug: int=0):

    if corepath[-1] != '/': corepath += '/'
    ### Reads all directories in a folder, searches for those containing both a Cell ".gmol" object and a ".cif" file.
    name_wo = ''.join([i for i in name if not i.isdigit()])

    crystal_filenames = []
    crystal_objects = []
    paths_list = []
    names_list = []
    cif_paths = []
    if debug > 0: print("FIND: Corepath:", corepath)
    for crystal in sorted(os.listdir(corepath)):
        if os.path.isdir(corepath+crystal):
            crystal_path = corepath+crystal+'/'
            cell_found = False
            cell_loaded = False
            cif_found = False

            # A first loop over the list of files, to ensure that the data is there
            for fil in sorted(os.listdir(crystal_path)):
                if fil.endswith(".gmol") and fil.startswith("Cell") and name_wo in fil:
                    cell_found = True
                    try: 
                        if debug > 0: print("Trying to load cell from", crystal_path+fil)
                        cell = load_binary(crystal_path+fil)
                        cell_loaded = True
                    except: 
                        if debug > 0: print("Cell could not be loaded from", crystal_path+fil)
                        else: pass
 
                if fil.endswith(".cif") and name_wo in fil:
                    cif_found = True
            if debug > 0: print("FIND: crystal_path", crystal_path, "results:", cell_found, cell_loaded, cif_found)

            ## If the information could be properly retrieved, then it stores it. 
            if cell_found and cell_loaded and cif_found: 
                for fil in sorted(os.listdir(crystal_path)):
                    if fil.endswith(".gmol") and fil.startswith("Cell") and name_wo in fil:
                        crystal_filenames.append(fil)
                        crystal_objects.append(load_binary(crystal_path+fil))
                        paths_list.append(crystal_path)
                        if crystal == "Original": crystal = name
                        names_list.append(crystal)
                    if fil.endswith(".cif") and name_wo in fil:
                        cif_paths.append(crystal_path+fil)

    return crystal_filenames, crystal_objects, paths_list, names_list, cif_paths

##############################
def create_sco_system(name, cell2mol_path: str, calcs_path: str, sys_path: str="default", debug: int=0) -> object:
    ## Creates a system. To do so, it...
    ## Reads all directories in cell2mol_path, searches for those containing both a Cell ".gmol" object and a ".cif" file.
    ## Each pair of ".gmol" and ".cif" will be considered a "crystal" and associates them to the "system" 

    ## calcs_path is the place where the analysis and computations will be done. 
    if cell2mol_path[-1] != '/': cell2mol_path += '/'
    if calcs_path[-1] != '/': calcs_path += '/'
    if sys_path == "default": sys_path = calcs_path
    if sys_path[-1] != '/':   sys_path += '/'
 
    # finds the crystals with Cell and cif file
    newsystem = sco_system(name, cell2mol_path, calcs_path, sys_path)
    if debug > 0: print("Searching crystals in", cell2mol_path)
    crystal_files, crystal_objects, crystal_paths, crystal_names, cif_paths  = find_crystals(name, cell2mol_path, debug=debug)
    if debug > 0: print("Found", len(crystal_files),"crystals in",cell2mol_path)
    # for each crystal:
    for idx, crys in enumerate(crystal_files):
        if debug > 0: print("Crystal",idx, name, crystal_names[idx], crystal_paths[idx]+crystal_files[idx])
        newcrystal = crystal(name, crystal_names[idx], crystal_paths[idx], crystal_files[idx], sys=newsystem)

        # Adds coord as variable. Legacy at some point
        #cell_postocoord(newcrystal.cell)
        # Stores reconstructed coordinates as a new coord and label variables. Original ones are stored separately
        #cell_reconstructed_coord(newcrystal.cell)

        if debug > 0: print(f"CREATE_SCO_SYSTEM: reading .cif data")
        newcrystal.read_cif_data(cif_paths[idx])
        if debug > 0: print(f"CREATE_SCO_SYSTEM: Getting FeN6 mols")
        newcrystal.get_FeN6_molecules(debug=debug)
        if debug > 0: print(f"CREATE_SCO_SYSTEM: Getting Phase data")
        newcrystal.get_spin_and_phase_data(debug=debug)
                            
        # If it worked, at the end it stores the data in the system object
        newsystem.crystals.append(newcrystal)

    # Determines Reference HS and LS molecules and crystals, if there are
    if debug > 0: print(len(newsystem.crystals), "crystals in system. Setting reference molecules and crystals")
    worked1 = newsystem.set_reference_molecs(debug=debug)
    worked2 = newsystem.set_reference_crystals(debug=debug)
    if worked1 and worked2: return newsystem
    else:                   return None
##############################

####################################
#### Functions to restart system ###
####################################
#def reset_paths_out(sys: object, cell2mol_path: str, calcs_path: str, sys_path: str, debug: int=0) -> None:
#    if os.path.isdir(cell2mol_path): sys.cell2mol_path        = cell2mol_path
#    if os.path.isdir(calcs_path):    sys.calcs_path           = calcs_path
#    if os.path.isdir(sys_path):      sys.sys_path             = sys_path
#    if debug > 0: print(f"RESET_PATHS_OUT: new system paths: {sys.cell2mol_path}, {sys.calcs_path}, {sys.sys_path}")
#    for br in sys.branches:
#        if os.path.isdir(calcs_path+br.keyword+'/'): br.path = calcs_path+br.keyword+'/'
#        else: print(f"RESET_PATHS_OUT: {calcs_path+br.keyword+'/'} path does not exist")
#        for rec in br.recipes:
#            rec.path = br.path
#            for job in rec.jobs:
#                job.path = rec.path
#                for comp in job.computations:
#                    if job.setup == "findiff" and os.path.isdir(job.path+"findiff"): comp.path = job.path+"findiff/"
#                    else:                                                            comp.path = job.path
#                    comp.inp_path = comp.path+comp.inp_name
#                    comp.out_path = comp.path+comp.out_name
#                    comp.sub_path = comp.path+comp.sub_name
#                    comp.check_files()
#                    if debug > 0: print(f"RESET_PATHS_OUT: new computation path: {comp.inp_path}")  
#
