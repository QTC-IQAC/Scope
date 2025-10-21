import os
from datetime import datetime
from Scope.Classes_Environment  import set_user
from Scope.Workflow.Recipe      import *

##########################
###### BRANCH CLASS ######
##########################
class branch(object):
    def __init__(self, path: str, name: str, _system: object, debug: int=0) -> None:
        self.type             = "branch"
        self.creation_time    = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        self.creation_user    = set_user()
        self.path             = path
        self.name             = name
        self.recipes          = []
        self._system          = _system
        self.isregistered     = False
        self.isgood           = False
        self.isfinished       = False
        self.results          = dict()
        self.status           = "active"

        ## Corrects self.path in case the user forgets to add '/' 
        if self.path[-1] != '/': self.path += '/'

    ########################################
    ### Add // Remove // Restart Recipes ###
    ########################################
    def reset_recipes(self):
        delattr(self,"recipes"); setattr(self,"recipes", [])
        return self.recipes

    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  self.results[result.key] = result

    def add_recipe(self, name: str):
        exists, new_recipe = self.find_recipe(name)
        if not exists: 
            exists, source = self._system.find_source(name)
            if exists: 
                new_recipe = recipe(name, source, _branch=self)
                self.recipes.append(new_recipe)
                return new_recipe
            else:
                print(f"BRANCH. Source with name {name} does not exist. Recipe could not be created")
                print(f"BRANCH. Current sources: {name} does not exist. Recipe could not be created")
    
    def remove_output_lines(self):
        for rec in self.recipes:
            for job in rec.jobs:
                for comp in job.computations:
                    comp.delete_lines()
    
    def find_recipe(self, name: str, debug: int=0):
        name = name.lower()
        if debug > 1: print(f"BRANCH.FIND_RECIPE: Searching Recipe with {name=}:") 
        for rec in self.recipes:
            if debug > 1: print(f"BRANCH.FIND_RECIPE: Comparing with",rec.name)
            if rec.name == name: return True, rec 
        return False, None

    ##############
    ### Status ###
    ##############
    def set_status(self, status: str):
        if not hasattr(self._system,"sys_path"):
            raise ValueError("BRANCH.SET_STATUS: system doesnt have sys_path")
        if status not in ['active','terminated','finished']:
            raise ValueError("BRANCH.SET_STATUS: status should be 'active','terminated' or 'finished'")
        self.status = status
        ## Create file for get_status
        filename = f"{self.name}_FINISHED"
        filepath = f"{self._system.sys_path}{filename}"
        if not os.path.exists(filepath): open(filepath, "a").close() # Should create an empty file

        return self.status
    
    def read_status(self):
        # There is an alternative to this function that does not require loading the system and branch. 
        # It is in Utils/Run_Job.py as get_status()
        sys_path = self._system.sys_path
        if   os.path.isfile(f"{sys_path}TERMINATED"):           self.status == "terminated"
        elif os.path.isfile(f"{sys_path}{self.name}_FINISHED"): self.status == "finished"
        return self.status

    def clear_status(self):
        sys_path = self._system.sys_path
        terminated_file = f"{sys_path}TERMINATED"
        finished_file   = f"{sys_path}{self.name}_FINISHED"
        if os.path.isfile(terminated_file): os.remove(terminated_file)
        if os.path.isfile(finished_file):   os.remove(finished_file)
    
    ####################
    ### Registration ###
    ####################
    def register(self, debug: int=0):
        if debug > 1: print("Registering Branch:", self.name)
        allgood     = True
        allfinished = True
        if len(self.recipes) > 0:
            for rec in self.recipes:
                if not rec.isregistered:                 rec.register(debug=debug)
                if not rec.isgood:                       allgood     = False
                if not rec.isfinished:                   allfinished = False
        else: 
            allgood  = False
            allfinished = False
        if allgood:     self.isgood     = True
        if allfinished: self.isfinished = True
        self.isregistered = True
        if debug > 1: print("Registered Branch:", self.name, "[REG, GOOD, FIN]", self.isregistered, self.isgood, self.isfinished)

    #############
    ### Other ###
    #############
    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> BRANCH                                      \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' System                = {self._system.name}\n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' self.status           = {self.status}\n'
        to_print += f' self.creation_time    = {self.creation_time}\n'
        to_print += f' self.creation_user    = {self.creation_user}\n'
        to_print += f' self.path             = {self.path}\n'
        to_print += f' self.name             = {self.name}\n'
        to_print += f' Num Recipes           = {len(self.recipes)}\n'
        to_print += '\n'
        return to_print
