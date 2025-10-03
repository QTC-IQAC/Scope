#####################################
##### Contains the SYSTEM Class #####
#####################################

import os
import numpy as np
from Scope.Read_Write import save_binary
from Scope.Workflow.Branch import branch

class system(object):
    def __init__(self, name: str, environment: object) -> None:
        self.version              = "1.0"
        self.type                 = "system"
        self.subtype              = "system"
        self.origin               = "created"
        self.name                 = name
        self.sources              = []

        ## Connection with Computational Workflokw
        self.results              = dict()
        self.branches             = []
        self.states               = []

        ## Paths
        self.read_paths_from_environment(environment)

    ######
    def __repr__(self, indirect: bool=False):
        to_print  = ''
        if not indirect: to_print += '---------------------------------\n'
        if not indirect: to_print += '   >>> SCOPE System Object >>>   \n'
        if not indirect: to_print += '---------------------------------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        to_print += f' Subtype               = {self.subtype}\n'
        to_print += f' Name                  = {self.name}\n'
        to_print += f' Source Path           = {self.sources_path}\n'  ## Path where files with molecular or cell structures are stored
        to_print += f' Calculations Path     = {self.calcs_path}\n'    ## Path where folders with calculations will be stored
        to_print += f' System File Path      = {self.sys_path}\n'      ## Path where the system object is stored
        if len(self.sources) > 0:
            to_print += '\n'
            to_print += f' # of Sources          = {len(self.sources)}\n'
            to_print += f'     idx, type, name, formula               \n'
            for idx, spec in enumerate(self.sources):
                to_print += f'     {idx}: {spec.type} {spec.name} {spec.formula} \n'
        if not indirect: to_print += '\n'
        return to_print

    ######
    def create_folders(self):
        """Create system and calculations folders if they do not exist."""
        for path in [self.sys_path, self.calcs_path]:
            try:
                os.makedirs(path, exist_ok=True)
            except Exception as exc:
                print(f"Error creating folder: {path}")
                print(exc)

    ######
    def save(self, filepath: str=None):
        if filepath is None: filepath = self.sys_file
        save_binary(self, filepath)

    ######
    def find_source(self, name: str, debug: int=0):
        if debug > 0: 
            print(f"FIND_SOURCE. Searching for source with {name} in system with {len(self.sources)} sources:")
            for sour in self.sources:
                print(f"    {sour.name} {sour.type}")
        if len(self.sources) == 0: return False, None
        for sour in self.sources:
            if sour.name.lower() == name.lower(): return True, sour 
        return False, None

    ######
    def add_source(self, name: str, new_source: object, overwrite: bool=False, debug: int=0):
        ## Links the system to the source
        new_source._sys = self
        ## Search if source with the same name already exists
        found, source = self.find_source(name, debug=debug)
        ## If not, it is added
        if not found: 
            self.sources.append(new_source)
        ## If it exists, it is overwritten if specified 
        elif overwrite: 
            self.sources = [s for s in self.sources if s.name.lower() != name.lower()]
            self.sources.append(new_source)
        else: 
            print(f"ADD_SOURCE: Source with name {new_source.name} already exists in system") 
            print(f"If you would like to Overwrite, specify overwrite=True")
        return self.sources

    #############
    ### Paths ###
    #############
    def check_paths(self, debug: int=0) -> bool:
        if not os.path.isfile(self.sys_file) or not os.path.isdir(self.sys_path) or not os.path.isdir(self.calcs_path) or not os.path.isdir(self.sources_path):  
            if debug > 0: print("SYSTEM.CHECK_PATHS: WARNING: folders do not exist")
            return False
        return True

    ######
    def read_paths_from_environment(self, environment: object, debug: int=0) -> None: 
        reset = False

        ## Fix for older versions
        if not hasattr(environment,"sources_path"):
            if hasattr(environment,"cell2mol_path"):  environment.sources_path = environment.cell2mol_path
            else: print("SYSTEM.READ_PATHS_FROM_ENV: WARNING: no sources_path in environment"); return False

        ## Environment Sends the Global Paths
        target_sources_path         = f"{environment.sources_path}{self.name}/"
        target_calcs_path           = f"{environment.calcs_path}{self.name}/"
        target_sys_path             = f"{environment.sys_path}{self.name}/"
        target_sys_file             = f"{environment.sys_path}{self.name}/{self.name}.npy"

        ## We make sure that those exist:
        #if not os.path.isdir(target_sys_path) or os.path.isdir(target_calcs_path) or os.path.isdir(target_sources_path):  
        #    print("SYSTEM.READ_PATHS_FROM_ENV: WARNING: folders do not exist")
        #    return False

        reset = True
        self.sources_path         = target_sources_path 
        self.calcs_path           = target_calcs_path 
        self.sys_path             = target_sys_path
        self.sys_file             = target_sys_file
        if debug > 0: 
            print(f"RESET_PATHS: new paths:")
            print(f"Source path: {self.sources_path}")
            print(f"Calcs path:  {self.calcs_path}")
            print(f"System path: {self.sys_path}")
            print(f"System file: {self.sys_file}")

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

    #########################################################
    ### Functions to Interact with Computational Workflow ###
    #########################################################
    def add_branch(self, name: str, debug: int=0):
        new_branch = branch(self.calcs_path+name, name, self, debug=debug)
        if not os.path.isdir(self.calcs_path+name): 
            try: os.makedirs(self.calcs_path+name)
            except Exception as exc:
                 print(f"Error creating branch folder in {self.calcs_path+name}")
                 print(exc)
        self.branches.append(new_branch)
        return new_branch

    ######
    def remove_branch(self, br_name=None):
        found = False
        for idx, br in enumerate(self.branches):
            if br.name == str(br_name): found = True; found_idx = idx
        if found:
            to_delete = self.branches[found_idx]
            del self.branches[found_idx]

    ######
    def remove_all_branches(self) -> None:
        if hasattr(self,"branches"): delattr(self,"branches"); setattr(self,"branches",[])
        return self.branches

    ######
    def find_branch(self, name: str, debug: int=0):
        if debug > 1: print(f"FIND_BRANCH. Finding branch with name:", name)
        if debug > 1: print(f"FIND_BRANCH. There are {len(self.branches)} branches in system")
        if len(self.branches) == 0: return False, None
        for idx, br in enumerate(self.branches):
            if debug > 1: print(f"FIND_BRANCH. Evaluating branch {idx} with name: {br.name} and path: {br.path}")
            if br.name.lower() == name.lower():
                if debug > 1: print(f"FIND_BRANCH. Branch was found. Checking path...")
                if not os.path.isdir(br.path) and debug > 0: 
                    print(f"WARNING: The path associated with this branch (below) does not exist. Loading the branch anyway")
                    print(f"WARNING: {br.path=}")
                return True, br
        return False, None

    ######
    def find_computation(self, branch_keyword: str, recipe_keyword: str, job_keyword: str, comp_keyword: str='', comp_step: int=1, comp_run_number: int=1, debug: int=0):
        if len(self.branches) == 0: return False, None
        for br in self.branches:
            if br.keyword.lower() == branch_keyword.lower():
                if len(br.recipes) == 0: return False, None
                for rec in br.recipes:
                    if rec.source.spin.lower() == recipe_keyword.lower():
                        if len(rec.jobs) == 0: return False, None
                        for job in rec.jobs:
                            if job.keyword.lower() == job_keyword.lower():
                                if len(job.computations) == 0: return False, None
                                for idx, comp in enumerate(job.computations):
                                    if comp.run_number == comp_run_number and comp.step == comp_step and comp.keyword.lower() == comp_keyword.lower(): return True, comp
        return False, None

    #########################################
    ### Functions to Interact with States ###
    #########################################
    def add_state(self, state: object, debug: int=0):
        ## Verifies that a state with the same name does not exist already. If not, appends it
        for my_state in self.states:
            if my_state.name == state.name: 
                found = True
                print("SYSTEM.ADD_STATE: you're trying to create a state that already exists in _source.")
                print("SYSTEM.ADD_STATE: use function called 'find_state' instead to retrieve existing state")
        if not found: self.states.append(state)
