import os
import numpy as np
from copy import deepcopy
from .Computation       import *
from ..Classes_State    import state, find_state
from ..Other            import where_in_array

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
        
    def check_job_data(self, job_path: str, debug: int=0):
        from ..Classes_Input import set_job_data, set_qc_data
        if debug > 0: print(f"CHECK_JOB_DATA: reading job_data from path: {job_path}")
        new_job_data    = set_job_data(job_path, section="&job_data" , debug=0)
        old_job_data    = self.job_data 
        if new_job_data != old_job_data: 
            print(f"CHECK_JOB_DATA: identified changes in job_data for job.keyword={self.keyword}")
            self.update_job_data(old_job_data, new_job_data) 
            for comp in self.computations:
                comp.check_qc_data(job_path=job_path, debug=debug)
            return True
        if debug > 0: print(f"CHECK_JOB_DATA: no changes in job_data")
        return False

    def update_job_data(self, old_job_data, new_job_data, debug: int=0):
        self.job_data         = new_job_data
        self.job_data        += old_job_data
        ## This is done to mimic the __init of a job class
        self.keyword          = new_job_data.keyword.lower()
        self.suffix           = new_job_data.suffix
        self.requisites       = new_job_data.requisites
        self.constrains       = new_job_data.constrains
        self.setup            = new_job_data.setup.lower()
        self.must_be_good     = new_job_data.must_be_good                
        self.check_requisites()

    def get_max_step(self, debug: int=0):
        max_step = 0
        for idx, comp in enumerate(self.computations):
            if comp.step > max_step: 
                max_step = comp.step
        return max_step

    def check_convergence(self, debug: int=0):
        from Scope.Other import check_convergence
        self.isconverged = check_convergence(self.energies, None, self.job_data.energy_thres)
        return self.isconverged

    def find_computation(self, keyword: str='', step: int=1, run_number: int=1, debug: int=0):
        for idx, comp in enumerate(self.computations):
            if not hasattr(comp,"step"): comp.step = 1
            if comp.keyword == keyword and comp.step == step and comp.run_number == run_number: this_comp = comp; return True, this_comp
        return False, None

    def add_computation(self, qc_data: object, step: int=1, path: str='', comp_keyword: str='', is_update: bool=False, debug: int=0):
        ## Name of the computation and file paths are not created automatically. Use self.set_name and self.set_paths
        if path == '': path == self.path
        new_computation       = computation(self, qc_data, step, path, comp_keyword, is_update=is_update, debug=debug)
        self.computations.append(new_computation)
        return new_computation 
    
    def remove_computation(self, comp_keyword=None, comp_step=None, comp_index=None):
        found = False
        for idx, comp in enumerate(self.computations):
            if comp_index is None and comp_step is None and comp_keyword is not None: 
                if comp.keyword == str(comp_keyword): found = True; found_idx = idx
            elif comp_index is not None and comp_step is None and comp_keyword is None: 
                if comp.index == int(comp_index): found = True; found_idx = idx
            elif comp_index is None and comp_step is not None and comp_keyword is None: 
                if comp.step == int(comp_step): found = True; found_idx = idx
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
        if debug > 1: print("Checking Requisites", self.requisites, "for job:",self.keyword)
        if debug > 1: print("Checking Constrains", self.constrains, "for job:",self.keyword)
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
        source = self._recipe.source

        #####################
        ## 1- Setup for regular computations: "1 job => 1 computation"
        #####################
        if self.setup == "regular" or self.setup == "reg":
            exists, new_comp = self.find_computation()
            if not exists: new_comp = self.add_computation(qc_data, 1, self.path, comp_keyword="", is_update=False, debug=debug)
            new_comp.qc_data._add_attr("istate",self.istate)         
            new_comp.qc_data._add_attr("fstate",self.fstate)         
            new_comp.set_name()
            new_comp.set_paths()

        #####################
        ## 2- Setup to create displaced geometries using VNM to escape from local minima
        #####################
        elif self.setup == "displacement" or self.setup == "disp":
            
            if debug > 0: print(f"SET COMPUTATIONS FROM SETUP: setting displacement starting from {self.istate}")
            exists, initial_state = find_state(source, self.istate)
            if not exists: print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' was not found")

            if hasattr(initial_state,"VNMs") and hasattr(initial_state,"coord"):
                #### Only applies to geometries that are not minimum

                if hasattr(initial_state,"isminimum"):
                    if initial_state.isminimum: 
                        if debug > 0: print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' is already a minimum")
                        self._recipe.remove_job(keyword=self.keyword)    # not sure if this is possible
                    else: 
                        ## 0-Checks that the VNM exists, and that they have eigenvalues. If not, registers those eigenvalues.
                        if not hasattr(initial_state.VNMs[0],"atomidxs"):  #Actually only checking for the first one, but its ok
                            print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' does not have VNMs with eigenvectors.") 
                            print(f"SET COMPUTATIONS FROM SETUP: now searching a job with frequencies in the state class object")
                            found = False
                            for idx, comp in enumerate(initial_state.computations):
                                print(idx, comp._job.keyword, comp.isgood)
                                if "freq" in comp._job.keyword and comp.isgood and not found:
                                    print(f"SET COMPUTATIONS FROM SETUP: I will try to read the eigenvectors from", idx, comp.out_path)
                                    ### 1-Parsing and storage (this block is similar to register_frequencies)
                                    source = comp._job._recipe.source
                                    if not hasattr(comp,"output"): reg_general(comp)
                                    VNMs = comp.output.get_vnms(witheigen=True)
                                    if VNMs is not None: 
                                        initial_state.set_VNMs(VNMs)
                                        found = True
                                        print(f"SET COMPUTATIONS FROM SETUP: Worked!")

                        ## 1-Displaces Coordinates Following Negative Freqs ###
                        from Scope.VNM_tools import displace_neg_freqs 
                        disp_coord = displace_neg_freqs(initial_state.coord, initial_state.VNMs,debug=debug)
                        exists, displ_state = find_state(source, "displaced")          # Checks if state already exists
                        if not exists: displ_state = state(source, "displaced")        # If not, creates it
                        displ_state.set_geometry(initial_state.labels, disp_coord)   # Creates New State with Displaced Coordinates 
                        if hasattr(initial_state,"cellvec"): 
                            displ_state.set_cell(initial_state.cellvec,initial_state.cellparam) 

                        ## 2-Initial State of the Computation must be updated, to account for the displacement of geometries
                        exists, new_comp = self.find_computation()
                        if not exists: new_comp = self.add_computation(qc_data, 1, self.path, comp_keyword="", is_update=False, debug=debug)
                        new_comp.qc_data._add_attr("istate","displaced")  ## Updates the initial state of the computation, so it takes the displaced geometries
                        new_comp.qc_data._add_attr("fstate",self.fstate)         
                        new_comp.set_name()
                        new_comp.set_paths()

            else:     
                print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' does not have the properties required to apply the displacement")
                self._recipe.remove_job(keyword=self.keyword)    # I'm trying to delete the job when it is not necessary
        
        #####################
        ## 3- Setup for finite Differences
        #####################
        elif self.setup == "findiff": 
            print(f"SET COMPUTATIONS FROM SETUP: findiff selected")
            found, initial_state = find_state(source, self.istate)
            if not found: print(f"SET COMPUTATIONS FROM SETUP: initial state not found")
            else:
                if hasattr(initial_state,"coord"):
                    from Scope.Findiff import findiff_displacements
                    #self.path = self.path+"findiff_test4"
                    if not os.path.isdir(self.path): os.makedirs(self.path); print(f"SET COMPUTATIONS FROM SETUP: findiff folder created")
                    geoms, names = findiff_displacements(initial_state.coord)

                    for idx, geo in enumerate(geoms):
                        assert len(geo) == len(initial_state.labels)
                        exists, displ_state = find_state(source, names[idx])          # Checks if state already exists
                        if not exists: displ_state = state(source, names[idx])        # If not, creates it and adds it to source
                        displ_state.set_geometry(initial_state.labels, geo)         # Creates New State with Displaced Coordinates 
                        if hasattr(initial_state,"cellvec"): displ_state.set_cell(initial_state.cellvec,initial_state.cellparam) 
    
                        # Initial State of the Computation must be updated, to account for the displacement of geometries
                        exists, new_comp = self.find_computation(keyword=names[idx])
                        if not exists: new_comp = self.add_computation(qc_data, 1, self.path, comp_keyword=names[idx], is_update=False, debug=debug)
                        new_comp.qc_data = deepcopy(new_comp.qc_data)
                        new_comp.qc_data._mod_attr("istate",names[idx]) ## Updates the initial state of the computation, so it takes the displaced geometries
                        new_comp.qc_data._mod_attr("fstate",names[idx]) ## Updates the initial state of the computation, so it takes the displaced geometries
                        new_comp.qc_data._mod_attr("print_forces",True) ## Updates the initial state of the computation, so it takes the displaced geometries
                        new_comp.set_name()
                        new_comp.set_paths()

                else: 
                    print("SET COMPUTATIONS FROM SETUP: ERROR! initial_state for findiff does not have coordinates")
                    print(initial_state)

        #####################
        ## 4- Setup for repetitive optimizations (rep_opt): consecutive opt jobs until energy is stable
        ##    Here, it only needs to create the first computation. 
        ##    Following ones are created in "self.set_continuation_computation()"
        #####################
        elif self.setup == "rep_opt":
            print(f"SET COMPUTATIONS FROM SETUP: repetitive optimization (rep_opt) selected")
            exists, new_comp = self.find_computation()    # Searches for first step and first run_number computation
            if not exists: new_comp = self.add_computation(qc_data, 1, self.path, comp_keyword="", is_update=False, debug=debug)
            new_comp.qc_data._add_attr("istate",self.istate)         
            new_comp.qc_data._add_attr("fstate",self.fstate)         
            new_comp.set_name()
            new_comp.set_paths()
            if not hasattr(self,"energies"): self.energies = np.zeros((self.job_data.max_steps))

        else: print(f"SET COMPUTATIONS FROM SETUP: {self.setup=} not recognized. No computations were created")

###############################
## Continuation Computations ##
###############################
    def set_continuation_computation(self, comp: object, typ: str, debug: int=0):
        comp.has_update = True
        if debug > 1: print("Set_Continuation_Comp: Creating new computation to continue job. Type:", typ, "Path:", comp.path)
    
        ###############################
        ## Operates Depending on typ ##
        ###############################
        if typ == "opt":
            # 0-We make sure that the new_run does not exist. If so, we return it directly:
            exists, new_comp = self.find_computation(keyword=comp.keyword, step=comp.step, run_number=comp.run_number+1)
            if exists: print("Set_Continuation_Comp: Continuation Computation exists"); return new_comp
            new_comp = self.add_computation(comp.qc_data, comp.step, comp.path, comp_keyword=comp.keyword, is_update=True, debug=debug)
            new_comp.qc_data = deepcopy(comp.qc_data)
            new_comp.set_name()
            new_comp.set_paths()

        elif typ == "scf":
            exists, new_comp = self.find_computation(keyword=comp.keyword, step=comp.step, run_number=comp.run_number+1)
            if exists: print("Set_Continuation_Comp: Continuation Computation exists"); return new_comp
            new_comp = self.add_computation(comp.qc_data, comp.step, comp.path, comp_keyword=comp.keyword, is_update=True, debug=debug)
            new_comp.qc_data = deepcopy(comp.qc_data)
            new_comp.set_name()
            new_comp.set_paths()
            # If scf has failed with quantum espresso, it is worth trying a different mixing_beta value 
            if new_comp.qc_data.software == "qe":
                import random
                updated = False
                old_value = new_comp.qc_data.mix_beta
                while not updated: 
                    new_value = round(random.uniform(0.2,0.8), 2)
                    if new_value != old_value: updated = True
                new_comp.qc_data._mod_attr("mix_beta",new_value)
                print("Set_Continuation_Comp: Mixing Beta changed to:", new_comp.qc_data.mix_beta)

        elif typ == "rep_opt":
            exists, new_comp = self.find_computation(keyword=comp.keyword, step=comp.step+1, run_number=1)
            if exists: print("Set_Continuation_Comp: Continuation Computation exists"); return new_comp
            new_comp = self.add_computation(comp.qc_data, comp.step+1, comp.path, comp_keyword=comp.keyword, is_update=False, debug=debug)
            new_comp.qc_data    = deepcopy(comp.qc_data)
            #new_comp.filename   = deepcopy(comp.filename)   # Filename contains how the file must be named (e.g. refcode+suffix+step+run_number...)
            new_comp.step       = comp.step + 1
            new_comp.set_name()
            new_comp.set_paths()
            print(f"Set_Continuation_Comp: Repetitive Optimization Created with {new_comp.name=}")

        else:
            print(f"Set_Continuation_Comp: received unknown type of continuation computation: {typ=}"); return None

        print(f"Set_Continuation_Comp: set new computation with {new_comp.run_number=} and {new_comp.step=}")
    
        ######
        ## Irrespectively of typ, the new computation will continue from the fstate, which should contain the latest available geometry and properties
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

        ## In theory, this elif block could be removed
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
        to_print += f' Source Type           = {self._recipe.source.type}\n'
        to_print += f' Source Spin           = {self._recipe.source.spin}\n'
        to_print += f' Branch Name           = {self._recipe._branch.name}\n'
        to_print += f' Recipe Name           = {self._recipe.name}\n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Job path              = {self.path}\n'
        to_print += f' Job keyword           = {self.keyword}\n'
        to_print += f' Job hierarchy         = {self.hierarchy}\n'
        to_print += f' Job requisites        = {self.requisites}\n'
        to_print += f' Job constrains        = {self.constrains}\n'
        to_print += f' Job setup             = {self.setup}\n'
        to_print += f' Num Computations      = {len(self.computations)}\n'
        if self.setup == "rep_opt": to_print += f' Num Steps             = {self.get_max_step()}\n'
        to_print += '----------------------------------------------------\n'
        to_print += f' self.isregistered (Temp) = {self.isregistered}\n'
        to_print += f' self.isgood       (Temp) = {self.isgood}\n' 
        to_print += f' self.isfinished   (Temp) = {self.isfinished}\n' 
        to_print += '\n'
        return to_print


