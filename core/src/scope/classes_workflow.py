import os
from copy import deepcopy
from datetime import datetime
from scope.classes_environment        import * 
from scope.classes_state              import State, find_state
from scope.operations.dicts_and_lists import where_in_array
from scope.register_data              import reg_general, reg_optimization, reg_frequencies, reg_energy, reg_excited_states
from scope.parse_general              import read_lines_file

##########################
###### BRANCH CLASS ######
##########################
class Branch(object):
    def __init__(self, path: str, name: str, _system: object, debug: int=0) -> None:
        self.object_type      = "branch"
        self.path             = path
        self.name             = name
        self._system          = _system

        self.creation_time    = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        self.creation_user    = set_user()

        self.workflows        = []
        self.isregistered     = False
        self.isgood           = False
        self.isfinished       = False
        self.results          = dict()
        self.status           = "active"

        ## Corrects self.path in case the user forgets to add '/' 
        if self.path[-1] != '/': self.path += '/'

    ###############
    ### Results ### 
    ###############
    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  self.results[result.key] = result

    def remove_result(self, key: str):
        return self.results.pop(key, None)

    #################
    ### Workflows ###
    #################
    def reset_workflows(self):
        delattr(self,"workflows"); setattr(self,"workflows", [])
        return self.workflows

    def add_workflow(self, name: str):
        ## Workflow names are de-capitalized and spaces replaced by underscores
        name = name.lower()
        name = name.replace(" ","_")
        exists, new_workflow = self.find_workflow(name)
        if not exists: 
            exists, source = self._system.find_source(name)
            if exists: 
                new_workflow = Workflow(name, source, _branch=self)
                self.workflows.append(new_workflow)
                return new_workflow
            else:
                print(f"BRANCH.ADD_WORKFLOW: Source with name {name} does not exist. Workflow could not be created")
        else:
            print(f"BRANCH.ADD_WORKFLOW: Workflow with the same name ({name}) already exists. Returning it")
    
    def find_workflow(self, name: str, debug: int=0):
        ## Workflow names are de-capitalized and spaces replaced by underscores
        name = name.lower()
        name = name.replace(" ","_")
        if debug > 1: print(f"BRANCH.FIND_WORKFLOW: Searching Workflow with {name=}:") 
        for wrk in self.workflows:
            if debug > 1: print(f"BRANCH.FIND_WORKFLOW: Comparing with",wrk.name)
            if wrk.name == name: return True, wrk 
        return False, None

    def remove_workflow(self, name: str, debug: int=0):
        found = False
        for idx, wf in enumerate(self.workflows):
            if wf.name == name and not found: found = True; found_idx = idx
        if found: del self.workflows[found_idx]

    ##############
    ### Status ###
    ##############
    def set_status(self, status: str):
        status = status.lower()
        if not hasattr(self._system,"system_path"):
            raise ValueError("BRANCH.SET_STATUS: system doesnt have system_path")
        if status not in ['active','terminated','finished']:
            raise ValueError("BRANCH.SET_STATUS: status should be 'active','terminated' or 'finished'")
        self.status = status

        ## Creates file when necessary to avoid having to load the system
        if status in ['terminated','finished']:
            if   status == 'finished':   filename = f"{self.name}_FINISHED"
            elif status == 'terminated': filename = f"{self.name}_TERMINATED"
            filepath = f"{self._system.system_path}{filename}"
            if not os.path.exists(filepath): open(filepath, "a").close() # Should create an empty file
        return self.status
    
    def read_status(self):
        # There is an alternative to this function that does not require loading the system and branch. 
        # It is in Utils/Run_Job.py as get_status()
        system_path = self._system.system_path
        if   os.path.isfile(f"{system_path}{self.name}_TERMINATED"):     self.status == "terminated"
        elif os.path.isfile(f"{system_path}{self.name}_FINISHED"):       self.status == "finished"
        return self.status

    def clear_status(self):
        system_path = self._system.system_path
        terminated_file = f"{system_path}{self.name}_TERMINATED"
        finished_file   = f"{system_path}{self.name}_FINISHED"
        if os.path.isfile(terminated_file): os.remove(terminated_file)
        if os.path.isfile(finished_file):   os.remove(finished_file)
    
    ####################
    ### Registration ###
    ####################
    def register(self, debug: int=0):
        if debug > 1: print("BRANCH.REGISTER: Registering Branch:", self.name)
        allgood     = True
        allfinished = True
        if len(self.workflows) > 0:
            for wrk in self.workflows:
                if not wrk.isregistered:                 wrk.register(debug=debug)
                if not wrk.isgood:                       allgood     = False
                if not wrk.isfinished:                   allfinished = False
        else: 
            allgood  = False
            allfinished = False
        if allgood:     self.isgood     = True
        if allfinished: self.isfinished = True
        self.isregistered = True
        if debug > 1: print("Registered Branch:", self.name, "[REG, GOOD, FIN]", self.isregistered, self.isgood, self.isfinished)

    ######
    def remove_output_lines(self):
        for wrk in self.workflows:
            for job in wrk.jobs:
                for comp in job.computations:
                    comp.delete_lines()

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
        to_print += f' Num Workflows         = {len(self.workflows)}\n'
        to_print += '\n'
        return to_print

############################
###### WORKFLOW CLASS ######
############################
class Workflow(object):
    def __init__(self, name: str, source: object, _branch: object, debug: int=0) -> None:
        self.object_type      = "workflow"
        self._branch          = _branch
        self.path             = _branch.path
        self.name             = name
        self.source           = source
        self.jobs             = []
        self.isregistered     = False
        self.isgood           = False
        self.isfinished       = False
        self.results          = dict()

    ###############
    ### Results ###
    ###############
    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  self.results[result.key] = result

    def remove_result(self, key: str):
        return self.results.pop(key, None)

    ############
    ### JOBS ### 
    ############
    def add_job(self, job_data, debug: int=0): ## As opposed to add_branch or add_workflow, add_job does not need a name, but a job_data input
        exists, new_job = self.find_job(job_data=job_data, debug=debug)
        if not exists:
            new_job = Job(job_data, _workflow=self)
            self.jobs.append(new_job)
        return new_job 

    def remove_job(self, name=None, hierarchy=None):
        found = False
        if name is None and hierarchy is None: print("WORKFLOW.REMOVE_JOB: Error removing job, please indicate either name, or hierarchy number")
        elif name is not None and hierarchy is not None: print("WORKFLOW.REMOVE_JOB: Error removing job, please indicate only name or hierarchy number, not BOTH")
        else:
            for idx, jb in enumerate(self.jobs):
                if name is not None and hierarchy is None:
                    name = str(name)
                    if jb.name == name and not found: found = True; found_idx = idx
                elif hierarchy is not None and name is None:
                    hierarchy = int(hierarchy)
                    if jb.hierarchy == hierarchy and not found: found = True; found_idx = idx
        if found: del self.jobs[found_idx]

    def find_job(self, name=None, hierarchy=None, job_data=None, debug: int=0):
        if name is None and hierarchy is None and job_data is not None:
            assert hasattr(job_data,"job") and hasattr(job_data,"hierarchy")
            if debug > 1: print(f"WORKFLOW.FIND_JOB: Searching Job with name: '{job_data.job}' and hierarchy '{job_data.hierarchy}'")
            for jb in self.jobs:
                if jb.name == job_data.job and jb.hierarchy == job_data.hierarchy:
                    if debug > 1: print(f"WORKFLOW.FIND_JOB: Job found")
                    return True, jb
        elif name is None and hierarchy is not None and job_data is None:  
            assert type(hierarchy) == int
            if debug > 1: print(f"WORKFLOW.FIND_JOB: Searching Job with and hierarchy '{hierarchy}'")
            for jb in self.jobs:
                if jb.hierarchy == hierarchy:
                    if debug > 1: print(f"WORKFLOW.FIND_JOB: Job found")
                    return True, jb
        elif name is not None and hierarchy is None and job_data is None:  
            assert type(name) == str
            if debug > 1: print(f"WORKFLOW.FIND_JOB: Searching Job with and name '{name}'")
            for jb in self.jobs:
                if jb.name == name:
                    if debug > 1: print(f"WORKFLOW.FIND_JOB: Job found")
                    return True, jb
        elif name is not None and hierarchy is not None and job_data is None:  
            assert type(name) == str and type(hierarchy) == int
            if debug > 1: print(f"WORKFLOW.FIND_JOB: Searching Job with name: '{name}' and hierarchy '{hierarchy}'")
            for jb in self.jobs:
                if jb.name == name and jb.hierarchy == hierarchy:
                    if debug > 1: print(f"WORKFLOW.FIND_JOB: Job found")
                    return True, jb
        return False, None

    ####################
    ### Registration ###
    ####################
    def register(self, debug: int=0):
        if debug > 0: print("WORKFLOW.REGISTER: Registering Workflow:", self.name)
        allgood     = True
        allfinished = True
        if len(self.jobs) > 0:
            for job in self.jobs:
                if not job.isregistered:                 job.register(debug=debug)
                if not job.isgood:                       allgood     = False
                if not job.isfinished:                   allfinished = False
        else: 
            allgood     = False
            allfinished = False
        if allgood:                 self.isgood       = True
        if allfinished:             self.isfinished   = True
        self.isregistered = True
        if debug > 0: print("WORKFLOW.REGISTER: Registered Workflow:", self.name, "[REG, GOOD, FIN]", self.isregistered, self.isgood, self.isfinished)

    def remove_output_lines(self):
        for job in self.jobs:
            for comp in job.computations:
                comp.delete_lines()

    #############
    ### Other ###
    #############
    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> >>> WORKFLOW                                \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Source Name                 = {self.source.name}\n'
        to_print += f' Source Type                 = {self.source.object_type}\n'
        to_print += f' Source sub-Type             = {self.source.object_subtype}\n'

        if hasattr(self.source,"charge"):    
            if self.source.charge is not None: to_print += f' Source Charge               = {self.source.charge}\n'
        if hasattr(self.source,"spin"):        
            if self.source.spin is not None:   to_print += f' Source Spin                 = {self.source.spin}\n'
        if hasattr(self.source,"phase"):       
            if self.source.phase is not None:  to_print += f' Source Phase                = {self.source.phase}\n'

        to_print += f'---------------------------------------------------\n'
        to_print += f' Workflow Name               = {self.name}\n'
        to_print += f' Num Jobs                    = {len(self.jobs)}\n'
        if len(self.jobs) > 0: 
            self.jobs.sort(key=lambda x: x.hierarchy)
            to_print += f'\tLast Job Name        = {self.jobs[-1].name}\n'
            to_print += f'\tLast Job Hierarchy   = {self.jobs[-1].hierarchy}\n'
        to_print += '\n'
        return to_print

#######################
###### JOB CLASS ######
#######################
class Job(object):
    def __init__(self, job_data: object, _workflow: object):        
        self.object_type      = "job"
        self._workflow        = _workflow
        self.source           = _workflow.source
        self.path             = _workflow.path
        self.job_data         = job_data    
        self.computations     = []
        self.isregistered     = False
        self.isgood           = False
        self.isfinished       = False
        
        ## I hate to do this; repeat variables from job_data
        self.name             = job_data.job
        self.type             = job_data.job_type
        self.istate           = job_data.istate
        self.fstate           = job_data.fstate
        self.hierarchy        = int(job_data.hierarchy)
        self.requisites       = job_data.requisites
        self.constrains       = job_data.constrains
        self.job_setup        = job_data.job_setup.lower()
        self.must_be_good     = job_data.must_be_good                

        ## Corrects self.path in case the user forgets to add '/' 
        if self.path[-1] != '/': self.path += '/'
        ## Corrects self.job_setup in case the user forgets to change
        if self.type == 'findiff' or self.type == 'findif': self.job_setup = 'findiff'
        
    def check_job_data(self, inp_path: str, debug: int=0):
        from scope.classes_input import set_input_data
        if debug > 0: print(f"CHECK_JOB_DATA: reading job_data from path: {inp_path}")
        _, _, new_job_data, _ = set_input_data(inp_path, debug=0)
        old_job_data    = self.job_data 
        if new_job_data != old_job_data: 
            print(f"CHECK_JOB_DATA: identified changes in job_data for job.name={self.name}")
            self.update_job_data(old_job_data, new_job_data) 
            for comp in self.computations:
                comp.check_qc_data(inp_path=inp_path, debug=debug)
            return True
        if debug > 0: print(f"CHECK_JOB_DATA: no changes in job_data")
        return False

    def update_job_data(self, old_job_data, new_job_data, debug: int=0):
        self.job_data         = new_job_data
        self.job_data        += old_job_data
        ## This is done to mimic the __init of a job class
        self.name             = new_job_data.job.lower()
        self.type             = new_job_data.job_type.lower()
        self.requisites       = new_job_data.requisites
        self.constrains       = new_job_data.constrains
        self.job_setup        = new_job_data.job_setup.lower()
        self.must_be_good     = new_job_data.must_be_good                
        self.check_requisites()

    def get_max_step(self, debug: int=0):
        max_step = 0
        for idx, comp in enumerate(self.computations):
            if comp.step > max_step: 
                max_step = comp.step
        return max_step

    def check_convergence(self, debug: int=0):
        from scope.other import check_convergence
        self.isconverged = check_convergence(self.energies, None, self.job_data.energy_thres)
        return self.isconverged

    def find_computation(self, keyword: str='', step: int=1, run_number: int=1, debug: int=0):
        for comp in self.computations:
            if not hasattr(comp,"step"): comp.step = 1
            if comp.keyword == keyword and comp.step == step and comp.run_number == run_number: this_comp = comp; return True, this_comp
        return False, None

    def add_computation(self, qc_data: object, step: int=1, path: str='', comp_keyword: str='', is_update: bool=False, debug: int=0):
        ## Name of the computation and file paths are not created automatically. Use self.set_name and self.set_paths
        if path == '': path = self.path
        new_computation = Computation(self, qc_data, step, path, comp_keyword, is_update=is_update, debug=debug)
        self.computations.append(new_computation)
        return new_computation 
    
    def remove_computation(self, comp_keyword=None, comp_step=None, comp_index=None, remove_files: bool=False, debug: int=0):
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
            if remove_files: 
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
        if debug > 1: print("JOB.CHECK_REQUISITES: Requisites", self.requisites, "for job:",self.name)
        if debug > 1: print("JOB.CHECK_REQUISITES: Constrains", self.constrains, "for job:",self.name)
        for idx, job in enumerate(self._workflow.jobs):

            if self != job:
                ## If necessary, it registers any related job
                if debug > 1: print("JOB.CHECK_REQUISITES: Evaluating Job with name:", job.name)
                if debug > 1: print("JOB.CHECK_REQUISITES: Evaluating Job, isregistered:", job.isregistered)
                if debug > 1: print("JOB.CHECK_REQUISITES: Evaluating Job, isgood:", job.isgood)
                if debug > 1: print("JOB.CHECK_REQUISITES: Evaluating Job, isfinished:", job.isfinished)
                if (job.name in self.requisites or job.name in self.constrains) and not job.isregistered: 
                    if debug > 1: print("JOB.CHECK_REQUISITES: Registering Previous Unregistered Job", job.name)
                    job.register(debug=debug)
                    if debug > 1: print("JOB.CHECK_REQUISITES: Registered Job while checking requisites", job.name)
                    if debug > 1: print(job.name, job.isregistered, job.isgood, job.isfinished)

                ## Evaluates Requisites and Constrains
                if job.name in self.requisites and job.isfinished:
                    if job.must_be_good and job.isgood: 
                        requisites_fulfilled[where_in_array(self.requisites,job.name)[0]] = 1
                        if debug > 1: print("JOB.CHECK_REQUISITES: Requisite: ", job.name, "fulfilled case 1")
                    if not job.must_be_good:
                        requisites_fulfilled[where_in_array(self.requisites,job.name)[0]] = 1
                        if debug > 1: print("JOB.CHECK_REQUISITES: Requisite: ", job.name, "fulfilled case 2")
                elif job.name in self.constrains and job.isfinished and job.isgood: 
                    constrains_fulfilled[where_in_array(self.constrains,job.name)[0]] = 1
                    if debug > 1: print("JOB.CHECK_REQUISITES: Constrain: ", job.name, "not fulfilled")
                elif job.name == self.name and 'self' in self.constrains and job.isfinished and job.isgood: 
                    constrains_fulfilled[where_in_array(self.constrains,'self')[0]] = 1
                    if debug > 1: print("JOB.CHECK_REQUISITES: Constrain: ", job.name, "in self not fulfilled")
                elif job.name not in self.requisites and job.name not in self.constrains:
                    if debug > 1: print("JOB.CHECK_REQUISITES: Unrelated or Unregisterd Job with name:", job.name)
                elif job.name in self.requisites and not job.isfinished:
                    if debug > 1: print("JOB.CHECK_REQUISITES: Job in Requisites has not finished:", job.name)

        # Takes Decision
        if all(a == 1 for a in requisites_fulfilled) or len(self.requisites) == 0: self.requisites_fulfilled = True
        if all(b == 0 for b in constrains_fulfilled) or len(self.constrains) == 0: self.constrains_fulfilled = True
        if self.requisites_fulfilled and self.constrains_fulfilled: 
            if debug > 1: print("JOB.CHECK_REQUISITES: Requisites fulfilled")
            return True
        else: 
            if not   self.requisites_fulfilled: 
                if debug > 1: print("JOB.CHECK_REQUISITES: Requisites NOT fulfilled")
            elif not self.constrains_fulfilled: 
                if debug > 1: print("JOB.CHECK_REQUISITES: Constrains NOT fulfilled")
            return False

    def remove_output_lines(self):
        for comp in self.computations:
            comp.delete_lines()

###############################
### Create Computations Set ###
###############################
    def set_computations_from_setup(self, qc_data: object, debug: int=0): 

        #####################
        ## 1- Setup for regular computations: "1 job => 1 computation"
        #####################
        if self.job_setup == "regular" or self.job_setup == "reg":
            exists, new_comp = self.find_computation()
            if not exists: new_comp = self.add_computation(qc_data, 1, self.path, comp_keyword="", is_update=False, debug=debug)
            new_comp.qc_data._add_attr("istate",self.istate)         
            new_comp.qc_data._add_attr("fstate",self.fstate)         
            new_comp.set_name()
            new_comp.set_paths()

        #####################
        ## 2- Setup to create displaced geometries using VNM to escape from local minima
        #####################
        elif self.job_setup == "displacement" or self.job_setup == "disp":
            
            if debug > 0: print(f"SET COMPUTATIONS FROM SETUP: setting displacement starting from {self.istate}")
            exists, initial_state = find_state(self.source, self.istate)
            if not exists: print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' was not found")

            if hasattr(initial_state,"VNMs") and hasattr(initial_state,"coord"):
                #### Only applies to geometries that are not minimum
                if not hasattr(initial_state,"isminimum"): initial_state.check_minimum()

                if initial_state.isminimum: 
                    if debug > 0: print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' is already a minimum")
                else: 
                    ## 0-Checks that the VNM exists, and that they have eigenvalues. If not, registers those eigenvalues.
                    if not hasattr(initial_state.VNMs[0],"atomidxs"):  #Actually only checking for the first one, but its ok
                        print(f"SET COMPUTATIONS FROM SETUP: initial state '{self.istate}' does not have VNMs with eigenvectors.") 
                        print(f"SET COMPUTATIONS FROM SETUP: now searching a job with frequencies in the state class object")
                        found = False
                        for idx, comp in enumerate(initial_state.computations):
                            print(idx, comp._job.name, comp.isgood)
                            if comp._job.type == "freq" and comp.isgood and not found:
                                print(f"SET COMPUTATIONS FROM SETUP: I will try to read the eigenvectors from", idx, comp.out_path)
                                ### 1-Parsing and storage (this block is similar to register_frequencies)
                                if not hasattr(comp,"output"): reg_general(comp)
                                VNMs = comp.output.get_vnms(witheigen=True)
                                if VNMs is not None: 
                                    initial_state.set_VNMs(VNMs)
                                    found = True
                                    print(f"SET COMPUTATIONS FROM SETUP: Worked!")

                    ## 1-Displaces Coordinates Following Negative Freqs ###
                    from scope.vnm_tools import displace_neg_freqs 
                    disp_coord = displace_neg_freqs(initial_state.coord, initial_state.VNMs,debug=debug)
                    exists, displ_state = find_state(self.source, "displaced")          # Checks if state already exists
                    if not exists: displ_state = State(self.source, "displaced")        # If not, creates it
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
                self._workflow.remove_job(name=self.name)    # I'm trying to delete the job when it is not necessary
        
        #####################
        ## 3- Setup for finite Differences
        #####################
        elif self.job_setup == "findiff": 
            print(f"SET COMPUTATIONS FROM SETUP: findiff selected")
            found, initial_state = find_state(self.source, self.istate)
            if not found: print(f"SET COMPUTATIONS FROM SETUP: initial state not found")
            else:
                if hasattr(initial_state,"coord"):
                    from scope.findiff import findiff_displacements
                    #self.path = self.path+"findiff_test4"
                    if not os.path.isdir(self.path): os.makedirs(self.path); print(f"SET COMPUTATIONS FROM SETUP: findiff folder created")
                    geoms, names = findiff_displacements(initial_state.coord)

                    for idx, geo in enumerate(geoms):
                        assert len(geo) == len(initial_state.labels)
                        exists, displ_state = find_state(self.source, names[idx])          # Checks if state already exists
                        if not exists: displ_state = state(self.source, names[idx])        # If not, creates it and adds it to source
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

        else: print(f"SET COMPUTATIONS FROM SETUP: {self.job_setup=} not recognized. No computations were created")

###############################
## Continuation Computations ##
###############################
    def set_continuation_computation(self, comp: object, typ: str, debug: int=0):
        comp.has_update = True
        if debug > 1: print("Set_Continuation_Comp: Creating new computation to continue job. Type:", typ, "Path:", comp.path)
    
        # 0-We make sure that the new_run does not exist. If so, we return it directly:
        exists, new_comp = self.find_computation(keyword=comp.keyword, step=comp.step, run_number=comp.run_number+1)
        if exists: print("Set_Continuation_Comp: Continuation Computation exists"); return new_comp

        ###############################
        ## Operates Depending on typ ##
        ###############################
        if typ == "opt":
            new_comp = self.add_computation(comp.qc_data, comp.step, comp.path, comp_keyword=comp.keyword, is_update=True, debug=debug)
            new_comp.qc_data = deepcopy(comp.qc_data)
            new_comp.set_name()
            new_comp.set_paths()
            update_fstate = True

        elif typ == "ts":
            ## In TS searches, if it doesn't work the first time, we go for calcAll
            new_comp = self.add_computation(comp.qc_data, comp.step, comp.path, comp_keyword=comp.keyword, is_update=True, debug=debug)
            new_comp.qc_data = deepcopy(comp.qc_data)
            new_comp.set_name()
            new_comp.set_paths()
            new_comp.qc_data._mod_attr("fctype","calcall")
            print("JOB.SET_CONTINUATION_COMP: recalcFC changed to calcAll in G16 Input")
            update_fstate = False ## In TS searches, it is preferable to retry from istate, rather than continuing from fstate

        elif typ == "scf":
            new_comp = self.add_computation(comp.qc_data, comp.step, comp.path, comp_keyword=comp.keyword, is_update=True, debug=debug)
            new_comp.qc_data = deepcopy(comp.qc_data)
            new_comp.set_name()
            new_comp.set_paths()
            update_fstate = True
            # If scf has failed with quantum espresso, it is worth trying a different mixing_beta value 
            if new_comp.qc_data.software == "qe":
                import random
                updated = False
                old_value = new_comp.qc_data.mix_beta
                while not updated: 
                    new_value = round(random.uniform(0.2,0.8), 2)
                    if new_value != old_value: updated = True
                new_comp.qc_data._mod_attr("mix_beta",new_value)
                print("JOB.SET_CONTINUATION_COMP: Mixing Beta changed to:", new_comp.qc_data.mix_beta)

        elif typ == "rep_opt":
            new_comp            = self.add_computation(comp.qc_data, comp.step+1, comp.path, comp_keyword=comp.keyword, is_update=False, debug=debug)
            new_comp.qc_data    = deepcopy(comp.qc_data)
            new_comp.step       = comp.step + 1
            new_comp.set_name()
            new_comp.set_paths()
            update_fstate = True
            print(f"JOB.SET_CONTINUATION_COMP: Repetitive Optimization Created with {new_comp.name=}")

        else: raise ValueError(f"JOB.SET_CONTINUATION_COMP: received unknown type of continuation computation: {typ=}")

        print(f"JOB.SET_CONTINUATION_COMP: set new computation with {new_comp.run_number=} and {new_comp.step=}")
    
        ######
        ## Irrespectively of typ, the new computation will continue from the fstate, which should contain the latest available geometry and properties
        #####
        if  hasattr(comp.qc_data,"fstate"):
            iscorrect = comp.verify_state(comp.qc_data.fstate, target='opt')
            if iscorrect and update_fstate:
                new_comp.qc_data._add_attr("istate",comp.qc_data.fstate)
                new_comp.qc_data._add_attr("fstate",comp.qc_data.fstate)
                print("JOB.SET_CONTINUATION_COMP: istate of new computation is modified to", new_comp.qc_data.istate)
                print("JOB.SET_CONTINUATION_COMP: fstate of new computation is modified to", new_comp.qc_data.fstate)
            else: 
                print("JOB.SET_CONTINUATION_COMP: istate of new computation remains as:", new_comp.qc_data.istate)
                print("JOB.SET_CONTINUATION_COMP: fstate of new computation remains as:", new_comp.qc_data.fstate)
    
        return new_comp

####################
### Registration ###
####################
    def register(self, debug: int=0):
        if debug > 1: print("JOB.REGISTER: Registering Job:", self.name)

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
#        if allgood and self.job_setup == 'findiff': 
#            if debug > 1: print("------------------------------------------")
#            if debug > 1: print("Registering Finite Differences of this job")
#            if debug > 1: print("------------------------------------------")
#            from scope.Register_Data import reg_findiff 
#            worked = reg_findiff(self)
#            if not worked: allgood = False

        if allgood:                                           self.isgood       = True
        if allfinished:                                       self.isfinished   = True
        self.isregistered = True
        if debug > 1: print("Registered Job:", self.name, "[REG, GOOD, FIN]", self.isregistered, self.isgood, self.isfinished)

#############
### Other ###
#############
    def __repr__(self):
        if not hasattr(self,"source"): self.source = self._workflow.source 
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> >>> >>> JOB                                 \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Source Type           = {self.source.object_type}\n'
        to_print += f' Source Name           = {self.source.name}\n'
        to_print += f' Branch Name           = {self._workflow._branch.name}\n'
        to_print += f' Workflow Name         = {self._workflow.name}\n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Job path              = {self.path}\n'
        to_print += f' Job name              = {self.name}\n'
        to_print += f' Job type              = {self.type}\n'
        to_print += f' Job hierarchy         = {self.hierarchy}\n'
        to_print += f' Job requisites        = {self.requisites}\n'
        to_print += f' Job constrains        = {self.constrains}\n'
        to_print += f' Job setup             = {self.job_setup}\n'
        to_print += f' Num Computations      = {len(self.computations)}\n'
        if self.job_setup == "rep_opt": to_print += f' Num Steps             = {self.get_max_step()}\n'
        to_print += '----------------------------------------------------\n'
        to_print += f' self.isregistered (Temp) = {self.isregistered}\n'
        to_print += f' self.isgood       (Temp) = {self.isgood}\n' 
        to_print += f' self.isfinished   (Temp) = {self.isfinished}\n' 
        to_print += '\n'
        return to_print

###########################
#### COMPUTATION CLASS ####
###########################
class Computation(object):
    def __init__(self, _job: object, qc_data: object, step: int, path: str, keyword: str, is_update: bool=False, debug: int=0):        
        self.object_type      = "computation"
        self._job             = _job       
        self.qc_data          = qc_data
        self.software         = qc_data.software
        self.type             = qc_data.comp_type     ## Type inherited from the parent job
        self.step             = step
        self.keyword          = keyword               ## Normally blank, but it can be used to identify the computation in jobs with complex setups
        self.path             = path
        self.index            = len(_job.computations)+1
        self.isregistered     = False
        self.has_update       = False
        self.is_update        = is_update
        self.run_number       = self.set_run_number() 
        self.states           = []
        self.source           = _job.source      ## Just to simplify calling this variable 

    #####################
    ### Name of Files ###
    #####################
    def set_file_extension(self):
        if self.software == 'g16':
            inp = ".com"
            out = ".log"
        if self.software == 'qe':
            inp = ".input"
            out = ".out"
        sub = ".sub"
        return inp, out, sub
 
    def get_mod_filename(self, mod_item_vars: list, mod_item_vals: list, debug: int=0):
        if not hasattr(self,"filename"): self.set_filename()
        new_filename = deepcopy(self.filename)
        found = False
        for idx, mod in enumerate(mod_item_vars):
            for jdx, item in enumerate(new_filename.items):
                if debug > 2: print(f"comp.GET_MOD_FILENAME: comparing {item.variable.lower()=} and {mod.lower()=}")
                if item.variable.lower() == mod.lower():
                    item.mod_value(mod_item_vals[idx])
                    if debug > 2: print(f"comp.GET_NAME_FROM_CONFIG: modifying {mod=} in new_filename")
        return new_filename

    def set_filename(self, use_sys_name: bool=True, use_sou_name: bool=True, use_job_name: bool=True, use_step: bool=False, use_run_number: bool=True, use_spin: bool=True, debug: int=0):
        ### Here is the convention I'm using to name files. It is better not to change once computations have been submitted
        ### Uses a filename-class object, as defined below, defined as a sum of items
        if not hasattr(self,"run_number"): self.set_run_number() 
        self.filename = Filename()   ## Class defined at the end of this file
        if use_sys_name:       new_item = Filename_item("sys_name",   self.source._sys.name);  self.filename.add_item(new_item)
        if use_sou_name:       new_item = Filename_item("sou_name",   self.source.name);       self.filename.add_item(new_item)
        if use_job_name:       new_item = Filename_item("job_name",   self._job.name);         self.filename.add_item(new_item)
        if use_step:           # Only step=2 and above are printed in name 
            new_item = filename_item("step",       self.step,'s')
            new_item.set_min_value(int(2))
            self.filename.add_item(new_item)
        if use_run_number:     new_item = Filename_item("run_number", self.run_number,'r'); self.filename.add_item(new_item)
        #if use_spin:           new_item = filename_item("spin",       self.spin);           self.filename.add_item(new_item)
        if self.keyword != '': new_item = Filename_item("keyword",    self.keyword);        self.filename.add_item(new_item)
        return self.filename

    def set_name(self, spacer: str='_', debug: int=0):
        if not hasattr(self,"filename"): 
            if debug > 0: print("COMPUTATION.set_name: filename not found. Creating it")
            if self.step > 1: self.set_filename(use_step=True)
            else:             self.set_filename(use_step=False)
        self.name = self.filename.get_name(spacer=spacer)
        return self.name

    def set_paths(self, new_path: str=None):
        if not hasattr(self,"filename"): self.set_name()
        if new_path is not None: self.path = new_path 
        if self.path[-1] != '/': self.path += '/'
        inp, out, sub = self.set_file_extension()
        # Filenames
        self.inp_name = ''.join([self.name,inp])
        self.out_name = ''.join([self.name,out])
        self.sub_name = ''.join([self.name,sub])
        # Paths
        self.inp_path = self.path+self.inp_name
        self.out_path = self.path+self.out_name
        self.sub_path = self.path+self.sub_name

    def check_files(self) -> None:
        if not hasattr(self,"inp_path"): self.set_paths()
        self.input_exists     = os.path.isfile(self.inp_path)
        self.output_exists    = os.path.isfile(self.out_path)
        self.subfile_exists   = os.path.isfile(self.sub_path)

        if self.input_exists:   self.input_modtime    = os.path.getmtime(self.inp_path)
        else:                   self.input_modtime    = int(0) 
        if self.output_exists:  self.output_modtime   = os.path.getmtime(self.out_path)
        else:                   self.output_modtime   = int(0) 
        if self.subfile_exists: self.subfile_modtime  = os.path.getmtime(self.sub_path)
        else:                   self.subfile_modtime  = int(0) 
        
    ##################################
    #### Update-related functions ####
    ##################################
    def check_updates(self, max_run_number: int=100, debug: int=0) -> int:
        ## Checks for updates in the computation
        self.has_update = False
        if debug > 1: print(f"COMP.CHECK_UPDATES: entering part 1: {self.has_update=}")

        ## 1-Searches in the job it is contained
        for comp in self._job.computations:
            if not hasattr(comp,"step"): comp.step = 1
            if comp.keyword == self.keyword and comp.step == self.step and hasattr(comp,"run_number"):
                if comp.run_number > self.run_number: self.has_update = True

        if debug > 1: print(f"COMP.CHECK_UPDATES: entering part 2: {self.has_update=}")
        ## 2-Checks for newer files with a similar filename in self.path
        if not self.has_update:
            inp, out, sub = self.set_file_extension()
            for rn in range(self.run_number, max_run_number):
                if debug > 1: print(f"COMP.CHECK_UPDATES: in part 2, trying: {rn=}")
                mod_filename = self.get_mod_filename(list(["run_number"]),list([rn]), debug=debug)  ## This creates a new version of the filename
                mod_name     = mod_filename.get_name()
                if debug > 1: print(f"COMP.CHECK_UPDATES: in part 2, searching: {mod_name=}")
                mod_path     = mod_filename.set_path(self.path)
                mod_outfile  = ''.join([mod_path,".out"])
                mod_exists   = os.path.isfile(mod_outfile)
                if mod_exists: 
                    self.has_update = True
                    if debug > 0: print(f"COMP.CHECK_UPDATES: found update with {mod_name=}")
        return self.has_update

    def set_run_number(self) -> int:
        run_number = 0
        for comp in self._job.computations:               # Searches in the job it is contained
            if comp.keyword == self.keyword and comp.step == self.step and hasattr(comp,"isfinished") and hasattr(comp,"run_number"):
                if comp.run_number > run_number: run_number = comp.run_number
        run_number += 1
        return run_number

    ###################################
    #### QC_DATA-related functions ####
    ###################################
    def check_qc_data(self, inp_path: str, debug: int=0):
        from scope.classes_input import set_input_data
        old_qc_data    = deepcopy(self.qc_data)
        _, _, _, new_qc_data = set_input_data(inp_path, debug=0)
        if new_qc_data != old_qc_data: 
            if debug > 1:
                print("COMP.CHECK_QC_DATA: found different qc_data:")
                print("COMP.CHECK_QC_DATA: it may be just either (i) the addition of defaults or (ii) the addition of the states information")
                print("--- OLD QC_DATA ---")
                print(old_qc_data)
                print("--- NEW QC_DATA ---")
                print(new_qc_data)
                print("----")
                print("CHECK_QC_DATA will now update:")
            self.update_qc_data(old_qc_data, new_qc_data) 
            return True
        return False

    def update_qc_data(self, old_qc_data, new_qc_data, debug: int=0):
        ## Updates qc_data
        self.qc_data          = new_qc_data
        self.software         = new_qc_data.software
        ## Adds any old information that is now not present. I'm not sure about this one
        self.qc_data         += old_qc_data
        self.type             = self.qc_data.comp_type
        ## Exceptions
        if self.is_update and self.software == 'qe' and old_qc_data.mix_beta != new_qc_data.mix_beta:
            self.qc_data._mod_attr("mix_beta",old_qc_data.mix_beta)
        return self.qc_data

    ##################################
    #### Output-related functions ####
    ##################################
    def create_output(self, debug: int=0):
        if not hasattr(self,'output_lines'): self.read_lines()
        ## Gaussian Computations
        if   self.software == 'g16': 
            from scope.software.gaussian.g16_output import G16_output
            allowed_types = ['specie']
            assert self.source.object_type in allowed_types
            self.output = G16_output(self.output_lines, self)
        ## Quantum Espresso Computations
        elif self.software == 'qe':  
            from scope.software.quantum_espresso.qe_output import QE_output
            allowed_types = ['specie', 'cell']
            assert self.source.object_type in allowed_types
            self.output = QE_output(self.output_lines, self)
        else: print(f"COMPUTATION.CREATE_OUTPUT: Output of {comp.software} computationss is not implemented."); return None
        return self.output 

    def delete_output(self) -> None:
        if hasattr(self,"output"): delattr(self,"output")

    def read_lines(self, flat: bool=True) -> None:        
        if not hasattr(self,"output_exists"): self.check_files()
        if self.output_exists: self.output_lines = read_lines_file(self.out_path, flat=flat)
        else:                  self.output_lines = []

    def delete_lines(self) -> None:
        if hasattr(self,"output_lines"): delattr(self,"output_lines")
        #self.output_lines = []

    #################################
    #### State-related functions ####
    #################################
    def add_state(self, state):
        if not hasattr(self,'states'): self.states = []
        found = False
        for st in self.states:
            if state.name == st.name: found = True
        if not found: self.states.append(state)

    def verify_state(self, name, target: str='opt'):
        found, state = self.source.find_state(name)
        if not found: return False
        #if not found: return None
        if target == 'opt':
            if hasattr(state,'coord') and hasattr(state,'labels'): return True
        else: 
            print("COMPUTATION.VERIFY STATE: target not implemented")
            return False 

    ######################################
    #### Submission-related functions ####
    ######################################
    def check_submission_status(self, environment: object, debug: int=0) -> None:
        if not hasattr(self,"name"): self.set_name() 
        key = self.name
        self.isrunning = environment.check_submitted(job_name=key, debug=debug)
        
    def add_submission_init(self, nprocs: int, queue: object) -> None:
        self.nprocs                = nprocs
        self.submission_queue      = queue.name
        self.submission_user       = queue._environment.user
        ## self.job_id is retrieved in self.submit

    def add_registration_data(self, user: str=set_user()) -> None:
        self.registration_time     = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        self.registration_user     = user
            
    ###################################
    #### Execute-related functions ####
    ###################################
    def store(self, debug: int=0) -> None:
        ### This function should be connected with the new environment.storage variable
        import shutil
        if self.output_exists: new_out_path= self.out_path+'_part'
        shutil.move(self.out_path, new_out_path)

###########################################
    def submit(self, environment: object, want_job_id: bool=False, debug: int=0) -> None:
        currdir = os.getcwd()
        ## Goes to path. Otherwise subprocess will fail
        os.chdir(self.path)

        ## If job_id is selected. Tries to get it and store it. Else, it just runs
        if not want_job_id: 
            subprocess.run(['bash','-c', f"{environment.commands['submit']} {self.sub_name}"])
            self.job_id = None
        else: 

            ### ALL below should go to environment
            res  = run_command(f"{environment.commands['submit']} {self.sub_name}")
            if not res.ok: raise RuntimeError(f"COMPUTATION.SUBMIT: Submission Failed: \n{res.command}\n{res.stderr}")
            text = res.stdout.rstrip().splitlines()[0]
                
            if environment.management_type == "sge":
                blocks = text.split()
                if len(blocks) == 7 and blocks[2].isdigit() and blocks[6] == 'submitted':   
                    self.job_id = int(blocks[2])
                    print(f"COMPUTATION.SUBMIT: job submitted with job_id: {self.job_id}")
                else:
                    print(f"COMPUTATION.SUBMIT: job submitted with unknown job_id. Blocks: {blocks}")
                    self.job_id = None

            elif environment.management_type == "slurm":
                blocks = text.split()
                if len(blocks) == 4 and blocks[3].isdigit() and blocks[0] == 'Submitted':   
                    self.job_id = int(blocks[3])
                    print(f"COMPUTATION.SUBMIT: job submitted with job_id: {self.job_id}")
                else:
                    print(f"COMPUTATION.SUBMIT: job submitted with unknown job_id. Blocks: {blocks}")
                    self.job_id = None
        os.chdir(currdir)  ## We go back to the original directory
 
###########################################
    def run(self, environment: object, options: object, debug: int=0) -> None:
        from scope.software.quantum_espresso.qe_input    import gen_qe_input, gen_qe_subfile 
        from scope.software.gaussian.g16_input           import gen_g16_input, gen_g16_subfile 

        ## 0-Checks that Resources are available
        if options.want_submit: sent_procs, sent_jobs = environment.get_user_requested(debug=debug)
        else:                   sent_procs = 0; sent_jobs = 0
        if sent_procs >= environment.max_procs or sent_jobs >= environment.max_jobs:
            if debug > 0: print(f"    Over maximum jobs/cores reached")
            return None

        ## 1-Gets Resources
        if options.want_submit:
            askqueue = environment.get_best_queue(debug=debug)
            askprocs = environment.requested_procs 
            ## 1.1-Adds Resources
            self.add_submission_init(nprocs=askprocs, queue=askqueue)
            ## 1.2-Creates Files
            self.check_files()
            if not self.input_exists or options.overwrite_inputs:
                if self.software == 'g16':  gen_g16_input(self, debug=0)
                elif self.software == 'qe': gen_qe_input(self, debug=0)
            if not self.subfile_exists or options.overwrite_inputs:
                if self.software == 'g16':  gen_g16_subfile(self, queue=askqueue, module=environment.g16_module, procs=askprocs, savechk=False)
                elif self.software == 'qe': gen_qe_subfile(self, queue=askqueue, module=environment.qe_module, procs=askprocs)

        ## 2-If output exists, prompts for registration
        if self.output_exists and not self.isregistered:
            print(f"    Output file Found Pending to be REGISTERED")
            print(f"    {self.output_path}")
            print(f"    ")

        ## 3-Evaluates Submission
        if options.want_submit:
            can_submit = True
            ## 3.1-Evaluates if output exists
            if self.output_exists and not options.overwrite_logs:
                can_submit = False
                if debug > 0: print("Output exists and not overwriting logs")

            ## 3.2-Evaluates if output is running, or ignores it if user sets options.ignore_submitted
            if can_submit and not options.ignore_submitted:
                self.check_submission_status(environment)   ### retrieves self.isrunning
                if self.isrunning:
                    can_submit = False
                    if debug > 0: print("Job already running")

            ## 3.3-Submits if possible
            if can_submit:
                self.submit(environment, debug=debug)

###########################################
    def register(self, debug: int=0) -> None:
        if debug > 0: print(f"COMP.REGISTER: Registering Computation with Job Name: {self._job.name}")
       
        ## Checks whether the output file exists:
        if not hasattr(self,"output_exists") or not hasattr(self,"output_modtime"): self.check_files()

        ## If so, it reads the lines under 2 conditions: 
        if self.output_exists:
            # 1-If lines have never been read: 
            if not hasattr(self,"output_lines"):                                                         self.read_lines()
            # 2-If the output file has been modified since job was last registered
            elif hasattr(self,"output_lines") and os.path.getmtime(self.out_path) > self.output_modtime: self.read_lines()

        ## 0-Creates Output, the object that will contain the parsing of data
        self.create_output(debug=debug)

        ## 1-Registration of General Attributes 
        reg_general(self, debug=debug)     # Gives self.isfinished, self.elapsed_time and self.status. It always works
     
        if not self.isfinished or self.status == "aborted":
            ## Retry with reading lines
            self.delete_output()
            self.read_lines()
            self.create_output(debug=debug)
            reg_general(self, debug=debug)  

        ## 2-Registration of Energy 
        worked = reg_energy(self, debug=debug)      # Stores the "last energy of a complete block" to State if it is not None

        ## 3-Registration of Optimization, Frequency and TD/TDA Tasks 
        opt_keywords = ['relax', 'vc-relax', 'opt', 'ts']
        if self.type in opt_keywords:
            worked = reg_optimization(self, debug=debug)
        elif self.type == 'freq': 
            worked = reg_frequencies(self, witheigen=False, debug=debug)
        elif self.type in ['td', 'tda']: 
            worked = reg_excited_states(self, debug=debug)

        ## 4-Wraps Up
        if worked:
            self.isregistered = True
            self.add_registration_data()
            self.delete_lines()   ## Output lines are deleted to save disk
            self.delete_output()  ## Output Object too, since it also stores output lines
        else: 
            print(f"COMP.REGISTER: Registration didn't work for: {self.out_path}")

        return worked

###########################################
    def __repr__(self) -> None:
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> >>> >>> >>> COMPUTATION                     \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Source Type           = {self.source.object_type}\n'
        to_print += f' Source sub-Type       = {self.source.object_subtype}\n'
        to_print += f' Branch Name           = {self._job._workflow._branch.name}\n'
        to_print += f' Workflow Name         = {self._job._workflow.name}\n'
        to_print += f' Job Name              = {self._job.name}\n'
        to_print += f' Job Type              = {self._job.type}\n'
        to_print += f'---------------------------------------------------\n'
        if hasattr(self,"istate"): to_print += f' Initial State         = {self.istate}\n'
        else:                      to_print += f' Initial State         = {self._job.istate}\n'
        if hasattr(self,"fstate"): to_print += f' Final State           = {self.fstate}\n'
        else:                      to_print += f' Final State           = {self._job.fstate}\n'
        to_print += f' Comp software         = {self.software}\n'
        to_print += f' Comp index            = {self.index}\n' 
        to_print += f' Comp step             = {self.step}\n' 
        to_print += f' Comp run_number       = {self.run_number}\n' 
        if self.keyword != '': to_print += f' Comp keyword          = {self.keyword}\n' 
        to_print += f' Comp inp_path         = {self.inp_path}\n' 
        to_print += f' Comp out_path         = {self.out_path}\n' 
        to_print += f' Comp isregistered     = {self.isregistered}\n' 
        if self.isregistered: to_print += f' Comp isgood           = {self.isgood}\n' 
        if self.isregistered: to_print += f' Comp isfinished       = {self.isfinished}\n' 
        if self.isregistered: to_print += f' Comp elapsed_time     = {self.elapsed_time} seconds\n' 
        if hasattr(self,"output"): to_print += f' Has OUTPUT            = YES\n'
        else:                      to_print += f' Has OUTPUT            = NO\n'
        to_print += '\n'
        return to_print

##############################################################
### FILENAME Class to facilitate controling the file names ###
##############################################################
class Filename(object):
    def __init__(self):
        self.typ      = 'name_global'
        self.items    = []
    
    def add_item(self, item):    
        if isinstance(item, Filename_item): self.items.append(item)

    def get_name(self, spacer: str='_', prefix: str='', suffix: str=''):
        self.name     = ''
        for idx, i in enumerate(self.items):
            if idx == 0: self.name += str(prefix)
            self.name += i.format()
            if idx == len(self.items)-1: self.name += str(suffix)
            else:                        
                if i.format() != '': self.name += spacer 
        return self.name

    def set_path(self, path: str=os.getcwd()): 
        if not hasattr(self,"name"): self.get_name()
        if path[-1] != '/': path += '/'
        self.path = path + self.name 
        return self.path

    def __repr__(self) -> None:
        to_print = ''
        for it in self.items:
            to_print += f'{str(it.variable)}: {str(it.value)}, Format: {it.format()}\n'
        return to_print

#######################
class Filename_item(object):
    ## Simple object to create filenames for computation files: input, output and submission
    def __init__(self, variable: str, value, prefix: str=''):
        self.typ              = 'filename_item'
        self.variable         = variable
        self.value            = value
        try:    self.value    = literal_eval(value)
        except: self.value    = value
        self.prefix           = prefix 

    def set_min_value(self, min_value):
        self.min_value        = min_value

    def mod_value(self, new_value):
        try:    self.value    = literal_eval(new_value)
        except: self.value    = new_value
        return self.value

    def format(self):
        to_print = ''
        if (type(self.value) == int or type(self.value) == float) and hasattr(self,"min_value"):
            if self.value >= self.min_value:
                if self.prefix != '': to_print += f'{self.prefix}'
                to_print += f'{str(self.value)}'
        else:
            if self.prefix != '': to_print += f'{self.prefix}'
            to_print += f'{str(self.value)}'
        return to_print

    def __repr__(self) -> None:
        to_print = ''
        if type(self.value) == int or type(self.value) == float:
            if hasattr(self,"min_value"):
                if self.value >= self.min_value:
                    if self.prefix != '': to_print += f'{self.prefix}'
                    to_print += f'{str(self.value)}'
        else:
            if self.prefix != '': to_print += f'{self.prefix}'
            to_print += f'{str(self.value)}'
        return to_print

    def __add__(self,other) -> None:
        if not isinstance(other, type(self)): return self
        return str(self.format()+other.format())
