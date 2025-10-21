from copy import deepcopy
import os
from datetime import datetime
from Scope.Classes_Spin                                 import *
from Scope.Classes_Environment                          import * 
from Scope.Register_Data                                import reg_general, reg_optimization, reg_frequencies, reg_energy
from Scope.Parse_General                                import read_lines_file

#########################
###### COMPUTATION ######
#########################
class computation(object):
    def __init__(self, _job: object, qc_data: object, step: int, path: str, keyword: str, is_update: bool=False, debug: int=0):        
        self.type             = "computation"
        self._job             = _job       
        self.qc_data          = qc_data
        self.software         = qc_data.software
        self.jobtype          = qc_data.jobtype       ## Type of computation (opt, scf, vc-relax)
        self.step             = step
        self.keyword          = keyword               ## Not the software, but a string used to identify the type of job computation
        self.suffix           = _job.suffix           ## A string to add to the file name. This should be eliminated 
        self.path             = path
        self.index            = len(_job.computations)+1
        self.isregistered     = False
        self.has_update       = False
        self.is_update        = is_update
        self.run_number       = self.set_run_number() 
        self.states           = []
        self.source           = _job._recipe.source      ## Just to simplify calling this variable 

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
        from ..Other import where_in_array, extract_from_list
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

    def set_filename(self, use_sys_name: bool=True, use_sou_name: bool=True, use_suffix: bool=True, use_step: bool=False, use_run_number: bool=True, use_spin: bool=True, debug: int=0):
        ### Here is the convention I'm using to name files. It is better not to change once computations have been submitted
        ### Uses a filename-class object, as defined below, defined as a sum of items
        if not hasattr(self,"run_number"): self.set_run_number() 
        self.filename = filename()   ## Class defined at the end of this file
        if use_sys_name:       new_item = filename_item("sys_name",   self.source._sys.name);  self.filename.add_item(new_item)
        if use_sou_name:       new_item = filename_item("sou_name",   self.source.name);       self.filename.add_item(new_item)
        if use_suffix:         new_item = filename_item("suffix",     self._job.suffix);       self.filename.add_item(new_item)
        if use_step:           # Only step=2 and above are printed in name 
            new_item = filename_item("step",       self.step,'s')
            new_item.set_min_value(int(2))
            self.filename.add_item(new_item)
        if use_run_number:     new_item = filename_item("run_number", self.run_number,'r'); self.filename.add_item(new_item)
        #if use_spin:           new_item = filename_item("spin",       self.spin);           self.filename.add_item(new_item)
        if self.keyword != '': new_item = filename_item("keyword",    self.keyword);        self.filename.add_item(new_item)
        return self.filename

    def set_name(self, spacer: str='_', debug: int=0):
        if not hasattr(self,"filename"): 
            if debug > 0: print("COMPUTATION.set_name: filename not found. Creating it")
            if self.step > 1: self.set_filename(use_step=True)
            else:             self.set_filename(use_step=False)
        self.name = self.filename.get_name(spacer=spacer)
        return self.name

    def set_paths(self):
        if not hasattr(self,"filename"): self.set_name()
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
    def check_qc_data(self, job_path: str, debug: int=0):
        from ..Classes_Input import set_qc_data
        old_qc_data    = deepcopy(self.qc_data)
        new_qc_data    = set_qc_data(job_path, section="&qc_data" , debug=0)
        if new_qc_data != old_qc_data: 
            if debug > 1:
                print("CHECK_QC_DATA found different qc_data:")
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
            from Scope.Software.Gaussian.G16_Class_Output import g16_output
            allowed_types = ['specie']
            assert self._job._recipe.source.type in allowed_types
            self.output = g16_output(self.output_lines, self)
        ## Quantum Espresso Computations
        elif self.software == 'qe':  
            from Scope.Software.Quantum_Espresso.QE_Class_Output import qe_output
            allowed_types = ['specie', 'cell']
            assert self._job._recipe.source.type in allowed_types
            self.output = qe_output(self.output_lines, self)
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
        source = self._job._recipe.source
        found, state = source.find_state(name)
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
        ## Goes to path. Otherwise subprocess will fail
        os.chdir(self.path)

        ## If job_id is selected. Tries to get it and store it. Else, it just runs
        if not want_job_id: 
            subprocess.run(['bash','-c', environment.command_submit+' '+self.sub_name])
            self.job_id = None
        else: 

            ### ALL below should go to environment
            raw = subprocess.check_output(['bash','-c', environment.command_submit+' '+self.sub_name])
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")[0]
                
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
 
###########################################
    def run(self, environment: object, options: object, debug: int=0) -> None:
        from Scope.Software.Quantum_Espresso.QE_Input    import gen_QE_input, gen_QE_subfile 
        from Scope.Software.Gaussian.G16_Input           import gen_G16_input, gen_G16_subfile 

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
                if self.software == 'g16':  gen_G16_input(self, debug=0)
                elif self.software == 'qe': gen_QE_input(self, debug=0)
            if not self.subfile_exists or options.overwrite_inputs:
                if self.software == 'g16':  gen_G16_subfile(self, queue=askqueue, module=environment.g16_module, procs=askprocs, savechk=False)
                elif self.software == 'qe': gen_QE_subfile(self, queue=askqueue, module=environment.qe_module, procs=askprocs)

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
     
        ## Retry with reading lines
        if not self.isfinished or self.status == "aborted":
            self.delete_output()
            self.read_lines()
            self.create_output(debug=debug)
            reg_general(self, debug=debug)  

        ## 2-Registration of Energy 
        worked = reg_energy(self, debug=debug)      # Stores the "last energy of a complete block" to State if it is not None

        ## 3-Registration of Optimization of Frequency Tasks 
        if 'opt' in self._job.keyword or 'relax' in self._job.keyword:
            worked = reg_optimization(self, debug=debug)
        elif 'freq' in self._job.keyword: 
            worked = reg_frequencies(self, witheigen=False, debug=debug)

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
        source = self._job._recipe.source
        to_print += f' Source Type           = {source.type}\n'
        to_print += f' Source sub-Type       = {source.subtype}\n'
        to_print += f' Branch Name           = {self._job._recipe._branch.name}\n'
        to_print += f' Recipe Name           = {self._job._recipe.name}\n'
        to_print += f' Job Keyword           = {self._job.keyword}\n'
        to_print += f'---------------------------------------------------\n'
        if hasattr(self,"istate"): to_print += f' Initial State         = {self.istate}\n'
        else:                      to_print += f' Initial State         = {self._job.istate}\n'
        if hasattr(self,"fstate"): to_print += f' Final State           = {self.fstate}\n'
        else:                      to_print += f' Final State           = {self._job.fstate}\n'
        to_print += f' Comp software         = {self.software}\n'
        to_print += f' Comp index            = {self.index}\n' 
        to_print += f' Comp step             = {self.step}\n' 
        to_print += f' Comp run_number       = {self.run_number}\n' 
        #to_print += f' Comp spin             = {self.spin}\n'
        to_print += f' Comp keyword          = {self.keyword}\n' 
        to_print += f' Comp inp_path         = {self.inp_path}\n' 
        to_print += f' Comp out_path         = {self.out_path}\n' 
        to_print += f' Comp isregistered     = {self.isregistered}\n' 
        if self.isregistered: to_print += f' Comp isgood           = {self.isgood}\n' 
        if self.isregistered: to_print += f' Comp isfinished       = {self.isfinished}\n' 
        if self.isregistered: to_print += f' Comp elapsed_time     = {self.elapsed_time}\n' 
        if hasattr(self,"output"): to_print += f'{self.output}\n'
        to_print += '\n'
        return to_print

##############################################################
### FILENAME Class to facilitate controling the file names ###
##############################################################
class filename(object):
    def __init__(self):
        self.typ      = 'name_global'
        self.items    = []
    
    def add_item(self, item):    
        if isinstance(item, filename_item): self.items.append(item)

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
            to_print += f'{str(it.variable)}: {str(it.value)}\n'
        return to_print

#######################
class filename_item(object):
    ## Simple object to create filenames for computation files: input, output and submission
    def __init__(self, variable: str, value, prefix: str=''):
        self.typ      = 'filename_item'
        self.variable = variable
        self.value    = value
        try:    self.value    = literal_eval(value)
        except: self.value    = value
        self.prefix   = prefix 

    def set_min_value(self, min_value):
        self.min_value = min_value

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
