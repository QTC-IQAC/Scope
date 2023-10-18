import sys
import copy
from copy import deepcopy
import os
import numpy as np
from datetime import datetime

from Scope.Classes_State import state, find_state
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
        self.job_data         = job_data    
        self.computations     = []
        self.isregistered     = False
        self.isgood           = False
        self.isfinished       = False
        
        ## I hate to do this; repeat variables from job_data
        self.keyword          = job_data.keyword
        self.istate           = job_data.istate
        self.fstate           = job_data.fstate
        self.hierarchy        = int(job_data.hierarchy)
        self.suffix           = job_data.suffix
        self.requisites       = job_data.requisites
        self.constrains       = job_data.constrains
        self.setup            = job_data.setup.lower()
        self.must_be_good     = job_data.must_be_good                

        ## Corrects self.path in case the user forgets to add '/' 
        if self.path[-1] != '/': self.path += '/'
        ## Corrects self.setup in case the user forgets to change
        if self.keyword == 'findiff' or self.keyword == 'findif': self.setup == 'findiff'
        
    def check_input(self, job_path: str, debug: int=0):
        from Scope.Classes_Input import set_job_data, set_qc_data
        ## job_data is for JOB
        new_job_data    = set_job_data(job_path, section="&job_data" , debug=0)
        old_job_data    = self.job_data 
        if new_job_data != old_job_data: 
            print(f"CHECK_INPUT: identified changes in job_data for job.keyword={self.keyword}")
            self.update_job_data(new_job_data) 
            new_qc_data    = set_qc_data(job_path, section="&qc_data" , debug=0)
            for comp in self.computations:
                comp.check_input(job_path=job_path, debug=debug)

    def update_job_data(self, new_job_data, debug: int=0):
        self.job_data         = new_job_data
        self.keyword          = new_job_data.keyword.lower()
        self.suffix           = new_job_data.suffix
        self.requisites       = new_job_data.requisites
        self.constrains       = new_job_data.constrains
        self.setup            = new_job_data.setup.lower()
        self.must_be_good     = new_job_data.must_be_good                

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
    
    def remove_computation(self, comp_keyword=None, comp_index=None):
        found = False
        for idx, comp in enumerate(self.computations):
            if comp_index is None and comp_keyword is not None: 
                if comp.keyword == str(comp_keyword): found = True; found_idx = idx
            elif comp_index is not None and comp_keyword is None: 
                if comp.index == int(comp_index): found = True; found_idx = idx
        if found: 
            to_delete = self.computations[found_idx]
            to_delete.check_files() 
            if to_delete.input_exists:   os.remove(to_delete.inp_path) 
            if to_delete.output_exists:  os.remove(to_delete.out_path) 
            if to_delete.subfile_exists: os.remove(to_delete.sub_path) 
            del self.computations[found_idx]
            
    def check_requisites(self, debug: int=0) -> None:
        self.requisites_fulfilled = False
        self.constrains_fulfilled = False
        requisites_fulfilled = np.zeros((len(self.requisites)))  ## To be correct, all must be 1
        constrains_fulfilled = np.zeros((len(self.constrains)))  ## To be correct, all must be 0
        if debug > 1: print("Checking Requisites", self.requisites)
        if debug > 1: print("Checking Constrains", self.constrains)
        for idx, job in enumerate(self._recipe.jobs):

            if self != job:
                ## If necessary, it registers any related job
                if debug > 1: print("Evaluating Job with keyword:", job.keyword)
                if debug > 1: print("Evaluating Job, isregistered:", job.isregistered)
                if debug > 1: print("Evaluating Job, isgood:", job.isgood)
                if debug > 1: print("Evaluating Job, isfinished:", job.isfinished)
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
                elif job.keyword not in self.requisites and job.keyword not in self.constrains:
                    if debug > 1: print("Unrelated or Unregisterd Job with keyword:", job.keyword)
                elif job.keyword in self.requisites and not job.isfinished:
                    if debug > 1: print("Job in Requisites has not finished:", job.keyword)

        # Takes Decision
        if all(a == 1 for a in requisites_fulfilled) or len(self.requisites) == 0: self.requisites_fulfilled = True
        if all(b == 0 for b in constrains_fulfilled) or len(self.constrains) == 0: self.constrains_fulfilled = True
        if self.requisites_fulfilled and self.constrains_fulfilled: 
            if debug > 1: print("Requisites fulfilled")
            return True
        else: 
            if not   self.requisites_fulfilled: 
                if debug > 1: print("Requisites NOT fulfilled")
            elif not self.constrains_fulfilled: 
                if debug > 1: print("Constrains NOT fulfilled")
            return False

    def remove_output_lines(self):
        for comp in self.computations:
            comp.delete_lines()

###############################
### Create Computations Set ###
###############################
    def set_computations_from_setup(self, qc_data: object, debug: int=0): 

        ## For Simplicity
        gmol = self._recipe.subject

        #####################
        ## 1- Setup for regular computations: "1 job => 1 computation"
        #####################
        if self.setup == "regular" or self.setup == "reg":
            exists, comp = self.find_computation()
            if not exists: comp = self.add_computation(int(1), qc_data, self.path, comp_keyword="", is_update=False, debug=debug)
            comp.qc_data._add_attr("istate",self.istate)         
            comp.qc_data._add_attr("fstate",self.fstate)         

        #####################
        ## 2- Setup to create displaced geometries using VNM to escape from local minima
        #####################
        elif self.setup == "displacement" or self.setup == "disp":
            
            if debug > 0: print(f"SET COMPUTATIONS FROM SETUP: setting displacement starting from {self.istate}")
            exists, initial_state = find_state(gmol, self.istate)
            if not exists: print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' was not found")

            #### Only applies to geometries that are not minimum
            if hasattr(initial_state,"coord") and hasattr(initial_state,"VNMs") and hasattr(initial_state,"isminimum"):
                if initial_state.isminimum: 
                    if debug > 0: print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' is already a minimum")
                    self._recipe.remove_job(keyword=self.keyword)    # not sure if this is possible
                else: 


                    ## 0-Checks that the VNM exists, and that they have eigenvalues. If not, registers those eigenvalues.
                    if hasattr(initial_state,"VNMs"):
                        if not hasattr(initial_state.VNMs[0],"atomidxs"):  #Actually only checking for the first one, but its ok
                            print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' does not have VNMs with eigenvectors.") 
                            found = False
                            for idx, comp in enumerate(initial_state.computations):
                                if "freq" in comp._job.keyword and comp.isgood and not found:
                                    print(f"SET COMPUTATIONS FROM SETUP: I will try to read the eigenvectors from", idx, comp.out_path)
                                    ### 1-Parsing and storage (this block is similar to register_frequencies)
                                    gmol = comp._job._recipe.subject
                                    if not hasattr(comp,"output"): reg_general(comp)
                                    VNMs = comp.output.get_vnms(witheigen=True)
                                    if VNMs is not None: 
                                        initial_state.set_VNMs(VNMs)
                                        found = True
                                        print(f"SET COMPUTATIONS FROM SETUP: Worked!")

                    ## 1-Displaces Coordinates Following Negative Freqs ###
                    from Scope.Gmol_ops import displace_neg_freqs 
                    disp_coord = displace_neg_freqs(initial_state.coord, initial_state.VNMs,debug=debug)
                    exists, displ_state = find_state(gmol, "displaced")          # Checks if state already exists
                    if not exists: displ_state = state(gmol, "displaced")        # If not, creates it
                    displ_state.set_geometry(initial_state.labels, disp_coord)   # Creates New State with Displaced Coordinates 

                    ## 2-Initial State of the Computation must be updated, to account for the displacement of geometries
                    exists, comp = self.find_computation()
                    if not exists: comp = self.add_computation(int(1), qc_data, self.path, comp_keyword="", is_update=False, debug=debug)
                    comp.qc_data._add_attr("istate","displaced")  ## Updates the initial state of the computation, so it takes the displaced geometries
                    comp.qc_data._add_attr("fstate",self.fstate)         

            else:     
                print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' does not have the properties required to apply the displacement")
                self._recipe.remove_job(keyword=self.keyword)    # I'm trying to delete the job when it is not necessary
        
        #####################
        ## 3- Setup for finite Differences
        #####################
        elif self.setup == "findiff": 
            initial_state = find_state(gmol, self.istate)

            if hasattr(initial_state,"coord"):
                from Scope.Findiff import findiff_displacements
                findiff_path = self.path+"findiff"
                geoms, names = findiff_displacements(initial_state.coord, debug=debug)
                for idx, geo in enumerate(geoms):
                    assert len(geo) == len(initial_state.labels)
                    exists, displ_state = find_state(gmol, names[idx])          # Checks if state already exists
                    if not exists: displ_state = state(gmol, names[idx])        # If not, creates it
                    displ_state.set_geometry(initial_state.labels, geo)         # Creates New State with Displaced Coordinates 

                    # Initial State of the Computation must be updated, to account for the displacement of geometries
                    exists, comp = self.find_computation(keyword=names[idx])
                    if not exists: comp = self.add_computation(idx, qc_data, findiff_path, comp_keyword=names[idx], is_update=False, debug=debug)
                    comp.qc_data._add_attr("istate",names[idx]) ## Updates the initial state of the computation, so it takes the displaced geometries
                    comp.qc_data._add_attr("fstate",names[idx]) ## Updates the initial state of the computation, so it takes the displaced geometries

        else: pass

###############################
## Continuation Computations ##
###############################
    def set_continuation_computation(self, comp: object, typ: str, debug: int=0):
        comp.has_update = True
    
        # 0-We make sure that the new_run does not exist. If so, we return it directly:
        exists, new_comp = self.find_computation(keyword=comp.keyword, index=comp.index+1)
        if exists:
            print("Set_Continuation_Comp: Continuation Computation exists")
            return new_comp
        else:
            if debug > 1: print("Set_Continuation_Comp: Creating new computation to continue job. Type:", typ)
    
        if typ == "opt":
            new_comp = self.add_computation(comp.index+1, comp.qc_data, path=comp.path, comp_keyword=comp.keyword, is_update=True, debug=debug)
            new_comp.qc_data = deepcopy(comp.qc_data)
        elif typ == "scf":
            new_comp = self.add_computation(comp.index+1, comp.qc_data, path=comp.path, comp_keyword=comp.keyword, is_update=True, debug=debug)
            new_comp.qc_data = deepcopy(comp.qc_data)
            if new_comp.qc_data.software == "qe":
                import random
                old_value = new_comp.qc_data.mix_beta
                new_value = round(random.uniform(0.2,0.8), 2)
                new_comp.qc_data._mod_attr("mix_beta",new_value)
                print("Set_Continuation_Comp: Mixing Beta changed to:", new_comp.qc_data.mix_beta)
        else:
            print("Set_Continuation_Comp: received unknown type of continuation computation: typ=", typ)
    
        ######
        ## Irrespectively of typ, the new computation will continue from the state, which should contain the latest available geometry and properties
        #####
        if  hasattr(comp.qc_data,"fstate"):
            iscorrect = comp.verify_state(comp.qc_data.fstate, target='opt')
            if iscorrect:
                new_comp.qc_data._add_attr("istate",comp.qc_data.fstate)
                new_comp.qc_data._add_attr("fstate",comp.qc_data.fstate)
                print("Set_Continuation_Comp: istate of new computation is modified to", new_comp.qc_data.istate)
                print("Set_Continuation_Comp: fstate of new computation is modified to", new_comp.qc_data.fstate)
            else: 
                print("Set_Continuation_Comp: istate of new computation remains as:", new_comp.qc_data.istate)
                print("Set_Continuation_Comp: fstate of new computation remains as:", new_comp.qc_data.fstate)
        elif hasattr(self,"fstate"):
            iscorrect = comp.verify_state(self.fstate, target='opt')
            if iscorrect:
                new_comp.qc_data._add_attr("istate",self.fstate)
                new_comp.qc_data._add_attr("fstate",self.fstate)
                print("Set_Continuation_Comp: istate of new computation is modified to", new_comp.qc_data.istate)
                print("Set_Continuation_Comp: fstate of new computation is modified to", new_comp.qc_data.fstate)
            else:
                print("Set_Continuation_Comp: istate of new computation remains as:", new_comp.qc_data.istate)
                print("Set_Continuation_Comp: fstate of new computation remains as:", new_comp.qc_data.fstate)
    
        return new_comp

                 

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
    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> >>> >>> JOB                                 \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Crystal               = {self._recipe.subject.refcode}\n'
        to_print += f' Type of Object        = {self._recipe.subject.type}\n'
        to_print += f' Spin                  = {self._recipe.subject.spin}\n'
        to_print += f' Recipe                = {self._recipe.keyword}\n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' self.path             = {self.path}\n'
        to_print += f' self.keyword          = {self.keyword}\n'
        to_print += f' self.hierarchy        = {self.hierarchy}\n'
        to_print += f' self.requisites       = {self.requisites}\n'
        to_print += f' self.constrains       = {self.constrains}\n'
        to_print += f' self.setup            = {self.setup}\n'
        to_print += f' Num Computations      = {len(self.computations)}\n'
        to_print += '----------------------------------------------------\n'
        to_print += f' self.isregistered (Temp) = {self.isregistered}\n'
        to_print += f' self.isgood       (Temp) = {self.isgood}\n' 
        to_print += f' self.isfinished   (Temp) = {self.isfinished}\n' 
        to_print += '----------------------------------------------------\n'
        return to_print
