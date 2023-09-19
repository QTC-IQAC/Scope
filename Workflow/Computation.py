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
from Scope.Write_G16_Inputs import *
from Scope.Write_QE_Inputs import *

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
        self.run_number       = self.set_run_number()
        self.isregistered     = False
        self.has_update       = False
        self.is_update        = is_update

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
        if debug > 0: print("SET_RUN_NUMBER: evaluating", self.keyword)
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
        if self.output_exists: self.output_lines = read_lines_file(self.out_path, flat=flat)
        else:                  self.output_lines = []

    def delete_lines(self) -> None:
        self.output_lines = []
            
    def check_submission_status(self, debug: int=0) -> None:
        key = str(self.refcode+self.suffix)
        self.isrunning = check_submitted(key, debug=debug)
        
    def add_submission_init(self, nprocs: str='Unk', queue: str='Unk', cluster: str=set_cluster(), user: str=set_user()) -> None:
        self.nprocs = nprocs
        self.queue = queue
        self.submission_cluster = cluster
        self.submission_user = user

    def add_registration_data(self, cluster: str=set_cluster(), user: str=set_user()) -> None:
        self.registration_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        self.registration_cluster = cluster
        self.registration_user = user

###########################################
    def store(self, debug: int=0) -> None:
        import shutil
        if self.output_exists: new_out_path= self.out_path+'_part'
        shutil.move(self.out_path, new_out_path)

###########################################
    def run(self, environment: object, options: object, debug: int=0) -> None:

        ## 0-Checks that Resources are available
        if options.want_submit: sent_procs, sent_jobs = check_usage()
        else:                   sent_procs = 0; sent_jobs = 0
        if sent_procs >= environment.max_procs or sent_jobs >= environment.max_jobs:
            if debug > 0: print(f"    Over maximum jobs and cores reached")
            return None

        ## 1-Gets Resources
        if options.want_submit:
            askqueue, askprocs = get_queue_and_procs(environment=environment)
            ## 1.1-Adds Resources
            self.add_submission_init(nprocs=askprocs, queue=askqueue)
            ## 1.2-Creates Files
            self.check_files()
            if not self.input_exists or options.overwrite_inputs:
                if self.software == 'g16':  gen_G16_input(self, debug=debug)
                elif self.software == 'qe': gen_QE_input(self, debug=debug)
            if not self.subfile_exists or options.overwrite_inputs:
                if self.software == 'g16':  gen_G16_subfile(self, procs=askprocs, queue=askqueue)
                elif self.software == 'qe': gen_QE_subfile(self, procs=askprocs, queue=askqueue)

        ## 2-If output exists, prompts for registration
        if self.output_exists and not self.isregistered:
            print(f"    Output file Found Pending to be REGISTERED")
            print(f"    {self.output_path}")
            print(f"    ")

        ## 2-Evaluates Submission
        if options.want_submit:
            can_submit = True
            ## 2.1-Evaluates if output exists
            if self.output_exists and not options.overwrite_logs:
                can_submit = False
                if debug > 0: print("Output exists and not overwriting logs")

            ## 2.2-Evaluates if output is running
            if can_submit and not options.ignore_submitted:
                self.check_submission_status()   ### retrieves self.isrunning
                if self.isrunning:
                    can_submit = False
                    if debug > 0: print("Job already running")

            ## 2.3-Submits if possible
            if can_submit:
                os.chdir(self.path)
                send_command("submit", filename=self.sub_name)
                #if debug > 0: print(f"Job {self.out_path} submitted")

###########################################
    def register(self, debug: int=0) -> None:
       
        ## Checks whether the output file exists
        if not hasattr(self,"output_exists") or not hasattr(self,"output_modtime"): self.check_files()

        ## If so, it reads the lines under 2 conditions: 
        if self.output_exists:
            ## 1-If lines have never been read: 
            if not hasattr(self,"output_lines"):                                                         self.read_lines()
            ## 2-If the output file has been modified since it was read
            elif hasattr(self,"output_lines") and os.path.getmtime(self.out_path) > self.output_modtime: self.read_lines()

        ## 1-Registration of General Attributes 
        reg_general(self, debug=debug)  # Gives self.isgood and self.isfinished

        ## 2-Registration of Optimization Tasks 
        if 'opt' in self._job.keyword or 'relax' in self._job.keyword:
            worked = reg_energy(self, debug=debug)
            if worked: worked = reg_optimization(self, debug=debug)
            else: 
                print(f"    WARNING: Update_Registry: Registration didn't work for: {self.out_path}")

        ## 3-Registration of Frequencies
        elif self.isgood and 'freq' in self._job.keyword:
            worked = reg_frequencies(self, debug=debug)
        else:  
            print(f"    WARNING: Update_Registry: Registration didn't work for: {self.out_path}")
            print(f"        -last line:", self.output_lines[-1])
            return False

        ## 4-Wraps Up
        if worked:
            self.isregistered = True
            self.add_registration_data()
        else:  
            print(f"    WARNING: Update_Registry: Registration didn't work for: {self.out_path}")
            if len(self.output_lines) > 0: print(f"        -last line:", self.output_lines[-1])

        return worked

###########################################
    def __repr__(self) -> None:
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> >>> >>> >>> COMPUTATION                     \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Crystal               = {self._job._recipe.subject.refcode}\n'
        to_print += f' Type of Object        = {self._job._recipe.subject.type}\n'
        to_print += f' Recipe                = {self._job._recipe.keyword}\n'
        to_print += f' Job                   = {self._job.keyword}\n'
        to_print += f' Initial State         = {self._job.istate}\n'
        to_print += f' Final State           = {self._job.fstate}\n'
        to_print += f'---------------------------------------------------\n'
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
        to_print += '----------------------------------------------------\n'
        return to_print
###########################################
