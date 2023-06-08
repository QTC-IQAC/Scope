import sys
import copy
from copy import deepcopy
import os
import numpy as np
from datetime import datetime

from Scope.Classes_Input import interpret_software
from Scope.Workflow import Computation
from Scope.Workflow.Computation import *
from Scope.Other import where_in_array

###################
######  JOB  ######
###################
class job(object):
    def __init__(self, job_data: object, _recipe: object):        
        self.type             = "job"
        self._recipe          = _recipe
        self.path             = _recipe.path
        self.keyword          = job_data.keyword.lower()
        self.job_data         = job_data    ## I hate to do this
        self.hierarchy        = int(job_data.hierarchy)
        self.software         = interpret_software(job_data.software)
        self.suffix           = job_data.suffix
        #self.environment      = job_data.environment
        self.requisites       = job_data.requisites
        self.constrains       = job_data.constrains
        self.setup            = job_data.setup.lower()
        self.must_be_good     = job_data.must_be_good                
        self.computations     = []
        self.isregistered     = False
        self.isgood           = False
        self.isfinished       = False
        #self.run_number       = 0   ## WARNING, I'm not sure about that

        ## Corrects self.path in case the user forgets to add '/' 
        if self.path[-1] != '/': self.path += '/'
        ## Corrects self.setup in case the user forgets to change
        if self.keyword == 'findiff' or self.keyword == 'findif': self.setup == 'findiff'
        
    #def set_run_number(self) -> int:
    #    self.run_number = 0
    #    for idx, jb in enumerate(self._recipe.jobs):  ## Searches in the recipe it is contained
    #        if jb.keyword == self.job_data.keyword and hasattr(jb,"isfinished") and hasattr(jb,"run_number"):
    #            if jb.isfinished and jb.run_number > self.run_number: self.run_number = jb.run_number
    #    self.run_number += 1
    #    if not self.job_data.must_be_good: self.run_number = 1
    #    return self.run_number

    def find_computation(self, keyword: str='', index=None):
        found = False
        for idx, comp in enumerate(self.computations):
            if index is None: 
                if comp.keyword == keyword and not found: this_comp = comp; found = True
            else:
                if comp.keyword == keyword and comp.index == int(index) and not found: this_comp = comp; found = True
        if found:    return found, this_comp
        else:        return found, None

    def add_computation(self, index: int, qc_data: object, path: str='', comp_keyword: str='', is_update: bool=False, debug: int=0):
        if path == '': path == self.path
        new_computation       = computation(index, comp_keyword, qc_data, path, _job=self, is_update=is_update, debug=debug)
        self.computations.append(new_computation)
        return new_computation 
    
    def remove_computation(self, comp_keyword: str):
        found = False
        for idx, comp in enumerate(self.computations):
            if comp.keyword == comp_keyword:  found = True; found_idx = idx
        if found: del self.computations[found_idx]
            
    def check_requisites(self, debug: int=0) -> None:
        self.requisites_fulfilled = False
        self.constrains_fulfilled = False
        requisites_fulfilled = np.zeros((len(self.requisites)))  ## To be correct, all must be 1
        constrains_fulfilled = np.zeros((len(self.constrains)))  ## To be correct, all must be 0
        if debug > 1: print("Checking Requisites", self.requisites)
        if debug > 1: print("Checking Constrains", self.constrains)
        for idx, job in enumerate(self._recipe.jobs):

            ## If necessary, it registers any related job
            if debug > 1: print("Evaluating Job with keyword:", job.keyword)
            if debug > 1: print("Evaluating Job, isregistered:", job.isregistered)
            if debug > 1: print("Evaluating Job, isgood:", job.isgood)
            if (job.keyword in self.requisites or job.keyword in self.constrains) and not job.isregistered: 
                if debug > 1: print("Registering Previous Unregistered Job", job.keyword)
                job.register(debug=debug)
                if debug > 1: print("Registered Job while checking requisites", job.keyword)
                if debug > 1: print(job.keyword, job.isregistered, job.isgood, job.isfinished)

            ## Evaluates Requisites and Constrains
            if job.keyword in self.requisites and job.isfinished:
                if job.must_be_good and job.isgood: 
                    requisites_fulfilled[where_in_array(self.requisites,job.keyword)[0]] = 1
                    if debug > 1: print("Requisite: ", job.keyword, "fulfilled 1")
                if not job.must_be_good:
                    requisites_fulfilled[where_in_array(self.requisites,job.keyword)[0]] = 1
                    if debug > 1: print("Requisite: ", job.keyword, "fulfilled 2")
            elif job.keyword in self.constrains and job.isfinished and job.isgood: 
                constrains_fulfilled[where_in_array(self.constrains,job.keyword)[0]] = 1
                if debug > 1: print("Constrain: ", job.keyword, "not fulfilled")
            elif job.keyword == self.keyword and 'self' in self.constrains and job.isfinished and job.isgood: 
                constrains_fulfilled[where_in_array(self.constrains,'self')[0]] = 1
                if debug > 1: print("Constrain: ", job.keyword, "in self not fulfilled")
            else:
                if debug > 1: print("Unrelated or Unregisterd Job with keyword:", job.keyword)

        # Takes Decision
        if all(a == 1 for a in requisites_fulfilled) or len(self.requisites) == 0: self.requisites_fulfilled = True
        if all(b == 0 for b in constrains_fulfilled) or len(self.constrains) == 0: self.constrains_fulfilled = True
        if self.requisites_fulfilled and self.constrains_fulfilled: print("Requisites fulfilled"); return True
        else: return False

###############################
### Create Computations Set ###
###############################
    def set_computations_from_setup(self, qc_data: object, debug: int=0): 
        ## Requires qc_data as external input

        ## Setup for regular computations: "1 job=1 computation"
        if self.setup == "regular" or self.setup == "reg":
            exists, comp = self.find_computation()
            if not exists: new_computation = self.add_computation(int(1), qc_data, self.path, comp_keyword="", is_update=False, debug=debug)

        ## Setup for finite Differences
        elif self.setup == "displacement" or self.setup == "disp":

            #### Only applies to geometries that are not minimum
            gmol = self._recipe.subject
            if hasattr(gmol,"isminimum"):
                if not gmol.isminimum: 
                    from Scope.Gmol_ops import displace_neg_freqs 
                    exists, comp = self.find_computation()
                    if not exists: new_computation = self.add_computation(int(1), qc_data, self.path, comp_keyword="", is_update=False, debug=debug)

                    ## Displaces Coordinates Following Negative Freqs ###
                    disp_coord = displace_neg_freqs(gmol,ini_coord_tag=qc_data.coord_tag,debug=debug)
                    newtag = "disp_coord"
                    if hasattr(gmol,newtag): gmol_update_geom(gmol, disp_coord, tag=newtag, debug=debug)
                    else:                    gmol_create_geom(gmol, disp_coord, tag=newtag, debug=debug)
                else: _recipe.remove_job(keyword=self.keyword)    # not sure if this is possible
            else:     _recipe.remove_job(keyword=self.keyword)    # I'm trying to delete the job when it is not necessary
        
        ## Setup for finite Differences
        elif self.setup == "findiff": 
            from Scope.Findiff import findiff_displacements
            backup_coord_tag = qc_data.coord_tag
            findiff_path = self.path+"findiff"
            geoms, names = findiff_displacements(qc_data.coord_tag, debug=debug)
            for idx, geo in enumerate(geoms):
                if hasattr(self._job._recipe._obj,names[idx]):  gmol_update_geom(self._job._recipe._obj, geo, tag=names[idx], debug=debug)
                else:                                           gmol_create_geom(self._job._recipe._obj, geo, tag=names[idx], debug=debug)
                qc_data.coord_tag = names[idx]

                # Computations are only added if they do not exist
                exists, comp = self.find_computation(keyword=names[idx])
                if not exists: new_computation = self.add_computation(idx, qc_data, findiff_path, comp_keyword=names[idx], is_update=False, debug=debug)

        else: pass
                 
####################
### Registration ###
####################
    def register(self, debug: int=0):
        if debug > 1: print("Registering Job:", self.keyword)

        ##########################################################################
        #### Irrespectively of the setup, we try to register all computations ####
        ##########################################################################
        allgood     = True
        allfinished = True
        if len(self.computations) > 0:
            for idx, comp in enumerate(self.computations):
                if not comp.has_update:    # if has_update means that there is another computation of the same type with a different run_number
                    if debug > 1: print("Registering Job: Evaluating computation with run_number:", comp.run_number)
                    comp.check_files()
                    if comp.output_exists:
                        #if debug > 1: print("Registering Job:", self.keyword, "output of:", comp.keyword, "exists")
                        if not comp.isregistered:                     worked = comp.register(debug=debug)
                        else:                                         worked = True
                        if not worked or not comp.isgood:             allgood     = False
                        if not worked or not comp.isfinished:         allfinished = False
                    else:                            
                        allgood     = False 
                        allfinished = False
        else:                            
            allgood     = False 
            allfinished = False

#        ############################################################
#        ## Findiff Setup: we extract frequencies at the job level ##
#        ############################################################
#        if allgood and self.setup == 'findiff': 
#            if debug > 1: print("------------------------------------------")
#            if debug > 1: print("Registering Finite Differences of this job")
#            if debug > 1: print("------------------------------------------")
#            from Scope.Register_Data import reg_findiff 
#            worked = reg_findiff(self)
#            if not worked: allgood = False

        if allgood:                                           self.isgood       = True
        if allfinished:                                       self.isfinished   = True
        self.isregistered = True
        if debug > 1: print("Registered Job:", self.keyword, "[REG, GOOD, FIN]", self.isregistered, self.isgood, self.isfinished)

#############
### Other ###
#############
    def get_info(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> >>> >>> JOB                                 \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Crystal               = {self._recipe.subject.refcode}\n'
        to_print += f' Type of Object        = {self._recipe.subject.type}\n'
        to_print += f' Spin                  = {self._recipe.subject.spin}\n'
        to_print += f' Recipe                = {self._recipe.keyword}\n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' self.path             = {self.path}\n'
        to_print += f' self.code             = {self.keyword}\n'
        to_print += f' self.hierarchy        = {self.hierarchy}\n'
        #to_print += f' self.run_number       = {self.run_number}\n'
        to_print += f' self.requisites       = {self.requisites}\n'
        to_print += f' self.constrains       = {self.constrains}\n'
        to_print += f' self.software         = {self.software}\n'
        #to_print += f' self.suffix           = {self.suffix}\n'
        #to_print += f' self.must_be_good     = {self.must_be_good}\n'
        to_print += f' self.setup            = {self.setup}\n'
        to_print += f' Num Computations      = {len(self.computations)}\n'
        #if len(self.computations) > 1: 
        #    self.computations.sort(key=lambda x: (x.index))
        #    to_print += f'     Last Computation Index = {self.jobs[-1].keyword}\n'
        #    if self.jobs[-1].run_number > 1: to_print += f' Last Job Run Number   = {self.jobs[-1].run_number}\n'
        to_print += '----------------------------------------------------\n'
        to_print += f' self.isregistered (Temp) = {self.isregistered}\n'
        to_print += f' self.isgood       (Temp) = {self.isgood}\n' 
        to_print += f' self.isfinished   (Temp) = {self.isfinished}\n' 
        to_print += '----------------------------------------------------\n'
        print(to_print)
