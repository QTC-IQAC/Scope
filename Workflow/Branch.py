from copy import deepcopy
import os
from datetime import datetime

from ..Classes_Environment  import set_cluster, set_user
from .Recipe    import *

##########################
###### BRANCH CLASS ######
##########################
class branch(object):
    def __init__(self, path: str, keyword: str, _sys: object, debug: int=0) -> None:
        self.type             = "branch"
        self.creation_time    = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        self.creation_cluster = set_cluster()
        self.creation_user    = set_user()
        self.path             = path
        self.keyword          = keyword
        self.recipes          = []
        self._sys             = _sys
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

    def add_recipe(self, keyword: str):
        exists, new_recipe = self.find_recipe(keyword)
        if not exists: 
            exists, source = self._sys.find_source(keyword)
            if exists: 
                new_recipe = recipe(keyword, source, _branch=self)
                self.recipes.append(new_recipe)
                return new_recipe
            else:
                print(f"BRANCH. Source with name {keyword} does not exist. Recipe could not be created")
                print(f"BRANCH. Current sources: {keyword} does not exist. Recipe could not be created")
    
    ########################################
    def find_recipe(self, keyword: str, debug: int=0):
        if debug > 1: print(f"BRANCH.FIND_RECIPE: Searching Recipe with keyword:", keyword) 
        for rec in self.recipes:
            if debug > 1: print(f"BRANCH.FIND_RECIPE: Comparing with",rec.keyword)
            if rec.keyword == keyword: return True, rec 
        return False, None

    ########################################
    def set_status(self, status: str):
        status = status.lower()
        if hasattr(self,"status"):
            if self.status != status: self.status = status
        elif not hasattr(self,"status"): 
            setattr(self,"status",status)

        ## Create file for get_status
        if hasattr(self._sys,"sys_path"):  path = self._sys.sys_path
        else: print("BRANCH.SET_STATUS: system doesnt have sys_path"); return None
            
        if   self.keyword.lower() == "isolated" and status == "finished": filename="ISO_FINISHED"
        elif self.keyword.lower() == "solid" and status == "finished":    filename="SOLID_FINISHED"
    
        filepath = str(path+'/'+filename)
        if not os.path.exists(filepath): open(filepath, "a").close() # Should create an empty file

    def remove_output_lines(self):
        for rec in self.recipes:
            for job in rec.jobs:
                for comp in job.computations:
                    comp.delete_lines()
    
####################
### Registration ###
####################
    def register(self, debug: int=0):
        if debug > 1: print("Registering Branch:", self.keyword)
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
        if debug > 1: print("Registered Branch:", self.keyword, "[REG, GOOD, FIN]", self.isregistered, self.isgood, self.isfinished)

#############
### Other ###
#############
    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> BRANCH                                      \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' System                = {self._sys.refcode}\n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' self.status           = {self.status}\n'
        to_print += f' self.creation_time    = {self.creation_time}\n'
        to_print += f' self.creation_cluster = {self.creation_cluster}\n'
        to_print += f' self.creation_user    = {self.creation_user}\n'
        to_print += f' self.path             = {self.path}\n'
        to_print += f' self.keyword          = {self.keyword}\n'
        to_print += f' Num Recipes           = {len(self.recipes)}\n'
        to_print += '\n'
        return to_print

#    def delete_inactive_recipes(self, obj, keyword, debug: int=0):
#        for idx, rec in enumerate(self.recipes):
#            if debug >= 1: print("")
#            if debug >= 1: print("    DEL_INACTIVE_REC: idx:", idx)
#            if debug >= 1: print("    DEL_INACTIVE_REC: recipe:", rec)
#            if debug >= 1: print("    DEL_INACTIVE_REC: keyword:", rec.keyword)
#            if hasattr(obj,"type") and debug >= 1:              print("    DEL_INACTIVE_REC: object type:", obj.type)
#            elif not hasattr(obj,"type") and debug >= 1:              print("    DEL_INACTIVE_REC: object type UNKNOWN")
#            if hasattr(rec.object,"type") and debug >= 1:       print("    DEL_INACTIVE_REC: rec.object type:", rec.object.type)
#            elif not hasattr(rec.object,"type") and debug >= 1: print("    DEL_INACTIVE_REC: rec.object type: cell")
#
#            if   rec.keyword == 'Solid' and not hasattr(rec.object,"type"): print("would remove:", rec.keyword, rec)
#            elif rec.keyword == keyword and rec.object != obj:
#                if rec.keyword == 'Solid':
#                    same_phase = rec.object.list_of_molecules[0].scope_guess_spin == obj.list_of_molecules[0].scope_guess_spin
#                    if debug >= 1: print("    DEL_INACTIVE_REC: number of molecules:", len(rec.object.list_of_molecules))
#                    if same_phase: print("would remove:", rec.keyword, rec)
##                self.recipes.remove(rec)
#                elif rec.keyword == 'Isolated':
#                    if rec.object.coord[0] == obj.coord[0]: print("would remove:", rec.keyword, rec)
