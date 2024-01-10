import sys
import copy
from copy import deepcopy
import os
import numpy as np
from datetime import datetime

from Scope.Classes_Spin import *
#from Scope.Classes_Input import interpret_software
from Scope.Environment import * 
#from Scope.Environment import set_cluster, set_user, check_submitted, check_usage, get_queue_and_procs, send_command
from Scope.Register_Data import reg_general, reg_optimization, reg_frequencies, reg_energy
from Scope.Software.Quantum_Espresso.Write_QE_Inputs import *
from Scope.Software.Gaussian.Write_G16_Inputs import *

#########################
###### COMPUTATION ######
#########################
class computation(object):
    def __init__(self, index: int, keyword: str, qc_data: object, path: str, _job: object, is_update: bool=False, debug: int=0):        
        self.type             = "computation"
        self._job             = _job       
        self.index            = index      
        self.keyword          = keyword  ## Not the software, but a string used to identify the computation
        self.qc_data          = qc_data
        self.software         = qc_data.software
        self.path             = path
        self.refcode          = _job._recipe.subject._sys.refcode
        self.run_number       = self.set_run_number(debug=0)
        self.isregistered     = False
        self.has_update       = False
        self.is_update        = is_update
        self.states           = []

        ############
        ### SPIN ###
        ############
        ## For the moment, self.spin is just a string with HS or LS 
        ## Spin Config can be used to define more complex spin states. 
        ## Particularly in crystals with intermediate spin states, which are not implemented yet
        ############
        if not hasattr(self.qc_data,"spin"): self.spin             = self._job._recipe.subject.spin
        else:                                self.spin             = qc_data.spin
        self.spin_config = get_spin_config(self._job._recipe.subject, self.spin, debug=debug)
        self.suffix           = str("_"+str(_job.suffix)+"_r"+str(self.run_number)+"_"+str(self.spin)+str(self.keyword))

        if self.path[-1] != '/': self.path += '/'

        # Filenames depend on Software
        if self.software == 'g16':
            inp_extension = ".com"
            out_extension = ".log"
        if self.software == 'qe':
            inp_extension = ".input"
            out_extension = ".out"
        sub_extension = ".sub"

        # Filenames
        self.inp_name = ''.join([self.refcode,self.suffix,inp_extension])
        self.out_name = ''.join([self.refcode,self.suffix,out_extension])
        self.sub_name = ''.join([self.refcode,self.suffix,sub_extension])
        # Paths
        self.inp_path = self.path+self.inp_name
        self.out_path = self.path+self.out_name
        self.sub_path = self.path+self.sub_name

    def check_input(self, job_path: str, debug: int=0):
        from Scope.Classes_Input import set_qc_data
        new_qc_data  = set_qc_data(job_path, section="&qc_data" , debug=0)
        old_qc_data  = self.qc_data 
        if new_qc_data != old_qc_data: 
            self.update_qc_data(new_qc_data) 

    def update_qc_data(self, new_qc_data, debug: int=0):
        self.qc_data          = new_qc_data
        self.software         = new_qc_data.software

    def set_run_number(self, debug: int=1) -> int:
        run_number = 0
        if debug > 0: print("SET_RUN_NUMBER: evaluating computation with keyword", self.keyword)
        for idx, comp in enumerate(self._job.computations):  ## Searches in the recipe it is contained
            if debug > 0: print("SET_RUN_NUMBER: in job there is:", comp.keyword, comp.isfinished, comp.run_number)
            if comp.keyword == self.keyword and hasattr(comp,"isfinished") and hasattr(comp,"run_number"):
                if comp.run_number > run_number: run_number = comp.run_number
                #if comp.isfinished and comp.run_number > run_number: run_number = comp.run_number
        run_number += 1
        if not self._job.job_data.must_be_good: run_number = 1
        return run_number

    def check_files(self) -> None:
        self.input_exists     = os.path.isfile(self.inp_path)
        self.output_exists    = os.path.isfile(self.out_path)
        self.subfile_exists   = os.path.isfile(self.sub_path)
        if self.input_exists:   self.input_modtime    = os.path.getmtime(self.inp_path)
        else:                   self.input_modtime    = int(0) 
        if self.output_exists:  self.output_modtime   = os.path.getmtime(self.out_path)
        else:                   self.output_modtime   = int(0) 
        if self.subfile_exists: self.subfile_modtime  = os.path.getmtime(self.sub_path)
        else:                   self.subfile_modtime  = int(0) 
        
    def read_lines(self, flat: bool=True) -> None:        
        if not hasattr(self,"output_exists"): self.check_files()
        if self.output_exists: self.output_lines = read_lines_file(self.out_path, flat=flat)
        else:                  self.output_lines = []

    def delete_lines(self) -> None:
        if hasattr(self,"output_lines"): delattr(self,"output_lines")
        #self.output_lines = []
            
    def check_submission_status(self, environment: object, debug: int=0) -> None:
        key = str(self.refcode+self.suffix)
        self.isrunning = environment.check_submitted(job_name=key, debug=debug)
        
    def add_submission_init(self, nprocs: int, queue: object) -> None:
        self.nprocs                = nprocs
        self.submission_queue      = queue.name
        self.submission_cluster    = queue._environment.cluster
        self.submission_user       = queue._environment.user
        ## self.job_id is retrieved in self.submit

    def add_registration_data(self, cluster: str=set_cluster(), user: str=set_user()) -> None:
        self.registration_time     = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        self.registration_cluster  = cluster
        self.registration_user     = user

    def create_output(self, debug: int=0):
        if not hasattr(self,'output_lines'): self.read_lines()
        if   self.software == 'g16': 
            from Scope.Software.Gaussian.G16_Class_Output import g16_output
            self.output = g16_output(self.output_lines, self)
        elif self.software == 'qe':  
            from Scope.Software.Quantum_Espresso.QE_Class_Output import qe_output
            self.output = qe_output(self.output_lines, self)
        else: print(f"COMPUTATION.CREATE_OUTPUT: Output of {comp.software} computationss is not implemented."); return None
        return self.output 

    def delete_output(self) -> None:
        if hasattr(self,"output"): delattr(self,"output")

    def add_state(self, state):
        if not hasattr(self,'states'): self.states = []
        self.states.append(state)

    def verify_state(self, name, target: str='opt'):
        subject = self._job._recipe.subject
        found, state = find_state(subject, name)
        if not found: return None
        if target == 'opt':
            if hasattr(state,'coord') and hasattr(state,'labels'): return True
        else: 
            print("COMPUTATION.VERIFY STATE: target not implemented")
            return False 
            
###########################################
    def store(self, debug: int=0) -> None:
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

        ## 0-Checks that Resources are available
        if options.want_submit: sent_procs, sent_jobs = environment.get_user_requested(debug=debug)
        else:                   sent_procs = 0; sent_jobs = 0
        if sent_procs >= environment.max_procs or sent_jobs >= environment.max_jobs:
            if debug > 0: print(f"    Over maximum jobs OR cores reached")
            return None

        ## 1-Gets Resources
        if options.want_submit:
            askqueue = environment.get_best_queue(debug=debug)
            askprocs = environment.requested_procs 
            #askqueue, askprocs = get_queue_and_procs(environment=environment)
            ## 1.1-Adds Resources
            self.add_submission_init(nprocs=askprocs, queue=askqueue)
            ## 1.2-Creates Files
            self.check_files()
            if not self.input_exists or options.overwrite_inputs:
                if self.software == 'g16':  gen_G16_input(self, debug=0)
                elif self.software == 'qe': gen_QE_input(self, environment, debug=0)
            if not self.subfile_exists or options.overwrite_inputs:
                if self.software == 'g16':  gen_G16_subfile(self, queue=askqueue, procs=askprocs)
                elif self.software == 'qe': gen_QE_subfile(self, queue=askqueue, procs=askprocs)

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

            ## 3.2-Evaluates if output is running
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
            ## 1-If lines have never been read: 
            if not hasattr(self,"output_lines"):                                                         self.read_lines()
            ## 2-If the output file has been modified since it was read
            elif hasattr(self,"output_lines") and os.path.getmtime(self.out_path) > self.output_modtime: self.read_lines()

        ## 0-Creates Output, the object that will contain the parsing of data
        self.create_output(debug=debug)

        ## 1-Registration of General Attributes 
        reg_general(self, debug=debug)     # Gives self.isfinished, self.elapsed_time and self.status. It always works
     
        ## 2-Registration of Energy 
        reg_energy(self, debug=debug)      # Stores the "last energy of a complete block" to State if it is not None
        #if not worked1: print(f"    COMP.REGISTER: Energy Registration didn't work for: {self.out_path}"); return False

        ## 3-Registration of Optimization of Frequency Tasks 
        if 'opt' in self._job.keyword or 'relax' in self._job.keyword:
        #if ('opt' in self._job.keyword or 'relax' in self._job.keyword) and worked1:
            worked2 = reg_optimization(self, debug=debug)
        elif 'freq' in self._job.keyword: 
        #elif self.isgood and 'freq' in self._job.keyword and worked1:
            worked2 = reg_frequencies(self, witheigen=False, debug=debug)
        else: 
            worked2 = True

        ## 4-Wraps Up
        if worked2:
            self.isregistered = True
            self.add_registration_data()
            self.delete_lines()   ## Output lines are deleted to save disk
            self.delete_output()  ## Output Object too, since it also stores output lines
        else: 
            print(f"    COMP.REGISTER: Opt/Freq Registration didn't work for: {self.out_path}")

        ### 5-Deletes Job_Id from queue pending
        #if hasattr(self,"job_id") and hasattr(self,"submission_queue"):

        return worked2

###########################################
    def __repr__(self) -> None:
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> >>> >>> >>> COMPUTATION                     \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Crystal               = {self._job._recipe.subject.refcode}\n'
        to_print += f' Type of Object        = {self._job._recipe.subject.type}\n'
        to_print += f' Recipe                = {self._job._recipe.keyword}\n'
        to_print += f' Job                   = {self._job.keyword}\n'
        to_print += f'---------------------------------------------------\n'
        if hasattr(self,"istate"): to_print += f' self.istate           = {self.istate}\n'
        else:                      to_print += f' Job Initial State     = {self._job.istate}\n'
        if hasattr(self,"fstate"): to_print += f' self.fstate           = {self.fstate}\n'
        else:                      to_print += f' Job Final State       = {self._job.fstate}\n'
        to_print += f' self.software         = {self.software}\n'
        to_print += f' self.index            = {self.index}\n' 
        to_print += f' self.spin             = {self.spin}\n'
        to_print += f' self.keyword          = {self.keyword}\n' 
        to_print += f' self.inp_path         = {self.inp_path}\n' 
        to_print += f' self.out_path         = {self.out_path}\n' 
        to_print += f' self.isregistered     = {self.isregistered}\n' 
        if self.isregistered: to_print += f' self.isgood           = {self.isgood}\n' 
        if self.isregistered: to_print += f' self.isfinished       = {self.isfinished}\n' 
        if self.isregistered: to_print += f' self.elapsed_time     = {self.elapsed_time}\n' 
        if hasattr(self,"output"): to_print += f'{self.output}\n'
        return to_print
###########################################
