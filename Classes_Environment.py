import sys
import os
import copy
import pwd
import grp 
import subprocess
import numpy as np

from Scope.Classes_Queue import queue 
from Scope.Classes_Input import input_data
from Scope import Classes_Input

def set_user():
    return pwd.getpwuid( os.getuid() ).pw_name

def set_group():
    group_id = pwd.getpwnam(set_user()).pw_gid
    return grp.getgrgid(group_id).gr_name

def set_cluster():
    return os.uname()[1]

###############
### CLUSTER ###
###############
class environment(object):
    def __init__(self):
        self.type                   = "environment"
        self.cluster                = set_cluster()
        self.user                   = set_user()
        self.group                  = set_group()
        self.available_queues       = [] 
        self.selected_queues        = [] 
        self.method                 = 'weighted'
        self.get_management_type() 

    def read_user_queue_list(self, line):
        ## Function to digest the queue names assuming that the user will use strange formats
 
        list_of_user_q = line.strip().split(",")
        for user_q in list_of_user_q:
            user_q = user_q.strip()
            for q in self.available_queues:
                if q.name == user_q or q.alter_name == user_q:
                    q.select_queue()
                elif user_q in q.name or q.name in user_q:
                    message = f"Found Similar queue as the one you requested. Is {q.name} your selection: {user_q} Y/N "
                    tmp = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"])
                    if tmp == "Y" or tmp == 'y':
                        q.select_queue()
                        q.alter_name = user_q

    def read_local_environment(self, file_path, debug: int=0):
        from Scope.Classes_Input import set_environment_data
        if not hasattr(self,"available_queues"): print("ENVIRONMENT.ADD: please run 'user_queue_preferences' first"); return None
        local_env = set_environment_data(file_path, debug=debug)
        if debug > 0: print("ENV.READ_USER_SPECS: reading data:", local_env)

        self.added_attr = {}
        for d in dir(local_env):

            ## General Attributes
            if d[0] != '_' and not callable(getattr(local_env,d)) and d != "dct" and d != "type" and d != "queues" and d != "queue":
                at1 = getattr(local_env,d)
                if debug > 0: print(f"ENV.READ_USER_SPECS: adding key={d}, value={at1}")
                self._add_attr(d,at1)
                self.added_attr[d] = at1 
            ## Queues are selected
            elif d == "queues" or d == "queue":
                at1 = getattr(local_env,d)
                self.read_user_queue_list(at1)
        return self.added_attr

    def check_paths(self, debug: int=0):
        if not hasattr(self,"cell2mol_path"): self.set_paths()
        if os.path.isdir(self.cell2mol_path): self.iscell2mol_path   = True 
        else:                                 self.iscell2mol_path   = False
        if os.path.exists(self.calcs_path):   self.iscalcs_path      = True
        else:                                 self.iscalcs_path      = False
        if os.path.exists(self.systems_path): self.issystems_path    = True 
        else:                                 self.issystems_path    = False
        if self.issystems_path and self.iscalcs_path and self.iscell2mol_path: return True
        else:                                                                  return False

    def set_commands(self):
        if self.management_type == "slurm":
            self.command_get_user_usage    = 'squeue -o "%.9P %.50j %.12u %.2t %.12M %.5C %.3D %R"'
            self.command_check_job           = 'squeue -o "%.60j %.12u"'
            self.command_submit              = 'sbatch'
        elif self.management_type == "sge":
            self.command_get_user_usage    = "qstat"
            self.command_get_user_waiting  = "qstat -f | grep ' qw '"
            self.command_check_job           = "qstat -xml | grep JB_name"
            self.command_submit              = "qsub"

    def set_PP_Library(self, debug: int=0):
        if   os.path.isdir(self.scope_home_path+"PP_Library"):       self.PP_Library= self.scope_home_path+"PP_Library/"
        elif os.path.isdir(self.scope_scratch_path+"PP_Library"):    self.PP_Library= self.scope_scratch_path+"PP_Library/"
        else:                                                        self.PP_Library= str(input("Please Specify PP_Library Path: ")) 
        print("PP_Library set to", self.PP_Library)
        return self.PP_Library
  
    def set_paths(self, debug: int=0):
        if 'login' in self.cluster or 'csuc' in self.cluster:
            self.scope_home_path    = f"/home/{self.user}/SCOPE/Database_SCO/" 
            self.scope_scratch_path = f"/scratch/{self.user}/SCOPE/Database_SCO/"
        elif 'lemma' in self.cluster:
            self.scope_home_path    = f"/Users/{self.user}/Documents/SCOPE/Database_SCO/"
            self.scope_scratch_path = None
        elif 'uam' in self.cluster:
            self.scope_home_path    = f"/home/proyectos/{self.group}/SCOPE/"
            self.scope_scratch_path = f"/scratch/{self.group}/"   #does not exist because they use ub100 instead of ub100435
        elif 'portal' in self.cluster or 'node' in self.cluster or 'visual' in self.cluster:
            if 'g2vela' in self.user:
                self.scope_home_path    = f"/home/{self.user}/SCOPE/Database_SCO/"
                self.scope_scratch_path = f"/scratch/{self.user}/SCOPE/Database_SCO/"
            elif 'g4vela' in self.user:
                self.scope_home_path    = f"/home/{self.user}/SCOPE/Database_SCO/"
                self.scope_scratch_path = f"/scratch/{self.user}/SCOPE/Database_SCO/"
        else:
            print(f"Cluster {self.cluster} not recognized")

        keyword = "4-Merged/"
        if   self.scope_home_path is not None and    os.path.isdir(self.scope_home_path+keyword):    self.cell2mol_path = self.scope_home_path+keyword
        elif self.scope_scratch_path is not None and os.path.isdir(self.scope_scratch_path+keyword): self.cell2mol_path = self.scope_scratch_path+keyword
        else:                                                self.cell2mol_path = str(input("Please Specify Cell2mol Path: ")) 

        keyword = "5-Complexes_Iso/"
        if   self.scope_home_path is not None and    os.path.isdir(self.scope_home_path+keyword):    self.calcs_path = self.scope_home_path+keyword
        elif self.scope_scratch_path is not None and os.path.isdir(self.scope_scratch_path+keyword): self.calcs_path = self.scope_scratch_path+keyword
        else:                                                self.calcs_path = str(input("Please Specify Calcs Path: ")) 

        keyword = "6-Systems_V3/"
        if   self.scope_home_path is not None and    os.path.isdir(self.scope_home_path+keyword):    self.systems_path = self.scope_home_path+keyword
        elif self.scope_scratch_path is not None and os.path.isdir(self.scope_scratch_path+keyword): self.systems_path = self.scope_scratch_path+keyword
        else:                                                self.systems_path = str(input("Please Specify Systems Path: ")) 
        
        if self.cell2mol_path[-1]   != '/': self.cell2mol_path += '/'
        if self.calcs_path[-1]      != '/': self.calcs_path    += '/'
        if self.systems_path[-1]    != '/': self.systems_path  += '/'

    def get_management_type(self, debug: int=0):
        self.management_type = "None"
        
        ### Sun Grid Engine, SGE ###
        worked_sge = False
        try: 
            sge_queues = subprocess.check_output(['bash','-c', "qconf -sql"], stderr=subprocess.DEVNULL)
            worked_sge = True
        except: pass
  
        ### Slurm Manager ###
        worked_slurm = False
        try: 
            slurm_queues = subprocess.check_output(['bash','-c', "sinfo"], stderr=subprocess.DEVNULL)
            worked_slurm = True
        except: pass

        if worked_sge and not worked_slurm:   self.management_type = "sge"
        elif not worked_sge and worked_slurm: self.management_type = "slurm"
        else: 
            print("GET_MANAGEMENT_TYPE: Confict with the recognition of the Queue System")
            print("GET_MANAGEMENT_TYPE: SGE:", worked_sge)
            print("GET_MANAGEMENT_TYPE: Slurm:", worked_slurm)
            print("GET_MANAGEMENT_TYPE: Assuming this is a local computer")
        return self.management_type

    def get_mqueues(self, debug: int=0):
        self.mqueues = []

        ##  SGE  ##
        if not hasattr(self,"management_type"): self.get_management_type()
        if self.management_type == "sge": 
            raw = subprocess.check_output(['bash','-c', "qconf -sql"]) 
            #raw = subprocess.check_output(['bash','-c', "qstat -g c"]) ## this could also be interesting 
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
            for idx, line in enumerate(text):
                name = line.split()[0] 
                raw2 = subprocess.check_output(['bash','-c', f"qstat -f | grep '{name}'"])
                dec2 = raw2.decode("utf-8")
                text2 = dec2.rstrip().split("\n")
                # Add queue
                new_queue = queue(name, self)
                self.add_mqueue(new_queue) 

        ## SLURM ##
        elif self.management_type == "slurm": 
            raw = subprocess.check_output(['bash','-c', "sinfo"]) 
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
            for idx, line in enumerate(text):
                if idx > 0 and len(line.split()) == 6: 
                    name, avail, time_limit, num_nodes, state, name_nodes = line.split()
                    new_queue = queue(name, self, avail=avail, time_limit=time_limit, state=state)
                    self.add_mqueue(new_queue) 
        return self.mqueues

    def add_mqueue(self, new_queue: object):
        if new_queue.name not in list(q.name for q in self.mqueues): 
            self.mqueues.append(new_queue)
            #if new_queue.num_nodes > 0: self.mqueues.append(new_queue)
        else: 
            for q in self.mqueues:
                if q.name == new_queue.name: q += new_queue

    def find_queue(self, queue_name: str):
        found = False
        for q in self.mqueues:
            if q.name == queue_name or q.name_alter == queue_name:
                return True, q
        return False, None

    def make_queue_available(self, queue_name: str):
        found, q = self.find_queue(queue_name)
        if found and q not in self.available_queues: self.available_queues.append(q)
        if not found: print(f"ENV.MAKE_QUEUE_AVAILABLE: queue {queue_name} not found")

    def make_queue_selected(self, queue_name: str):
        found, q = self.find_queue(queue_name)
        if found and q not in self.selected_queues: self.selected_queues.append(q)
        if not found: print(f"ENV.MAKE_QUEUE_AVAILABLE: queue {queue_name} not found")

    def save(self, filepath: str):
        self.filepath    = filepath
        from Scope.Read_Write import save_binary
        save_binary(self, filepath)

#####################################
###  Connection with Execute_Job  ###
#####################################
    def get_user_requested(self, debug: int=0):
        if not hasattr(self,"user_waiting_cpus"): self.get_user_waiting(debug=debug)
        if not hasattr(self,"user_running_cpus"): self.get_user_running(debug=debug)
        self.user_requested_cpus = self.user_waiting_cpus + self.user_running_cpus
        self.user_requested_jobs = self.user_waiting_jobs + self.user_running_jobs
        return self.user_requested_cpus, self.user_requested_jobs

    def get_user_waiting(self, debug: int=0):
        if not hasattr(self,"command_get_user_waiting"): self.set_commands()
        self.user_waiting_cpus = int(0)
        self.user_waiting_jobs = int(0)

        try: 
            raw = subprocess.check_output(['bash','-c', self.command_get_user_waiting])
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
        except Exception as exc:
            self.user_waiting_cpus = int(0)
            self.user_waiting_jobs = int(0)
            return self.user_waiting_cpus, self.user_waiting_jobs

        ## SGE clusters
        if self.management_type == "sge":
            try:
                for line in text:
                    blocks = line.split()
                    if len(blocks) == 8:
                        if blocks[4] == 'qw':
                            self.user_waiting_cpus  += int(blocks[7])
                            self.user_waiting_jobs  += 1
            except Exception as exc:
                self.user_waiting_cpus = int(0)
                self.user_waiting_jobs = int(0)
        return self.user_waiting_cpus, self.user_waiting_jobs

    def get_user_running(self, method: str='direct', debug: int=0):
        if not hasattr(self,"command_get_user_waiting"): self.set_commands()
        if len(self.available_queues) == 0: 
            print("CHECK_USER_USAGE: please set available queues first by running the method 'set_queues'") 
            return None
        self.user_running_cpus = 0
        self.user_running_jobs = 0

        ## Method 1 - Queue by Queue
        if method == 'queues': 
            for idx, q in enumerate(self.available_queues):
                cpus, jobs     = q.get_user_running(debug=debug)
                self.user_running_cpus += cpus
                self.user_running_jobs += jobs
            ## Finally, we add the waiting jobs
            #wcpus, wjobs = self.get_user_waiting(debug=debug)
            #self.user_running_cpus += wcpus
            #self.user_running_jobs += wjobs

        ## Method 2 - Directly
        elif method == 'direct':
            raw = subprocess.check_output(['bash','-c', self.command_get_user_usage])
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
            ## SGE clusters
            if self.management_type == "sge":
                try:
                    for line in text:           
                        blocks = line.split()
                        if len(blocks) == 9:
                            if blocks[4] == 'r':
                                self.user_running_cpus  += int(blocks[8])
                                self.user_running_jobs  += 1
                except Exception as exc:
                    self.user_running_cpus += int(0)
                    self.user_running_jobs += int(0)
                    print("ENVIRONMENT.CHECK_USER_USAGE: exception:", exc)
        
            ## SLURM clusters
            elif self.management_type == "slurm":
                try:
                    for line in text:
                        if q.name in line or q.alter_name in line:
                            blocks = line.split()
                            self.user_running_cpus      += int(blocks[5])
                            self.user_running_jobs      += 1
                except Exception as exc:
                    self.user_running_cpus = int(0)
                    self.user_running_jobs = int(0)
                    print("ENVIRONMENT.CHECK_USER_USAGE: exception:", exc)
        return self.user_running_cpus, self.user_running_jobs
         
    def get_best_queue(self, autoselect: bool=False, debug: int=0):
   
        ## If user has not selected queues yet. Then we take all available queues
        if len(self.selected_queues) == 0: 
            print(f"GET_BEST_QUEUE: queues have not been selected. Taking all available queues")
            target_queues = getattr(self,"available_queues")
            if autoselect:
                for aq in self.available_queues:
                    aq.select_queue()
        elif len(self.available_queues) > 0:
            target_queues = getattr(self,"selected_queues")
        else:
            print(f"GET_BEST_QUEUE: neither selected nor available queues were found")
            return None

        # Retrieve score
        scores = []
        for idx, q in enumerate(target_queues):
            tmp = q.get_queue_score(method=self.method)
            scores.append(tmp)
            if debug > 0: print(f"GET_BEST_QUEUE: evaluated {q.name} with score {tmp}")

        # Select best
        best_idx = np.argmax(scores) 
        if debug > 0: print(f"GET_BEST_QUEUE: selected best_idx={best_idx} of {scores}")
        if debug > 0: print(f"GET_BEST_QUEUE: returning {target_queues[best_idx].name}")
        return target_queues[best_idx]

    def check_submitted(self, job_name: str, debug: int=0):
        if not hasattr(self,"command_check_job"): self.set_commands()
        if not hasattr(self,"management_type"): self.get_management_type()
        if self.management_type is "None": return False 

        raw  = subprocess.check_output(['bash','-c', self.command_check_job])
        dec  = raw.decode("utf-8")
        flat = dec.replace("\n", "")
        if self.management_type == "sge": flat = flat.replace("<","").replace(">","").replace("/","").replace("JB_name","").split()
        if str(job_name) in flat: found = True
        else:                     found = False
        return found

#########################
###  User Interaction ###
#########################
    def set_queues(self, debug: int=0):
        from Scope.Read_Write import read_user_input
        if not hasattr(self,"mqueues"):   self.get_mqueues()

        suggested_names = list(q.name for q in self.mqueues if q.available)
        print(f"Setting Queues For Environment")

        if len(suggested_names) > 0: 
            print(f"-------------------------------------------------")
            print(f"         The following Queues were found:        ")
            print(f"-------------------------------------------------")
            for n in suggested_names:
                print(n)
            print(f"-------------------------------------------------")
            message = "Do you want to keep ALL suggested Queues? Y/N "
            print(" ")
            keep = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"]) 

            if keep == "Y" or keep == 'y':        
                user_q = list(q.name for q in self.mqueues) 
                verify  = False
                correct = True
            elif keep == "N" or keep == 'n':      
                verify = True
                user_q = []
                message = "Write the Queues/Partitions that will be available to SCOPE. Write them one by one: Enter to Stop  "
                finished = False
                while not finished: 
                    tmp = read_user_input(message=message, rtype=True, rtype_options=[str]) 
                    if tmp != '':    
                        interpreted_q = tmp.strip().replace(",","").split(" ")
                        for iq in interpreted_q:
                            if iq != '': user_q.append(iq.strip())
                    else: finished = True
            else: correct = False
                    
        elif len(suggested_names) == 0: 
            verify = True
            user_q = []
            message = "Write the Queues/Partitions that will be available to SCOPE. Write them one by one: Enter to Stop  "
            finished = False
            while not finished: 
                tmp = read_user_input(message=message, rtype=True, rtype_options=[str]) 
                if tmp != '':    
                    interpreted_q = tmp.strip().replace(",","").split(" ")
                    for iq in interpreted_q:
                        if iq != '': user_q.append(iq.strip())
                else: finished = True
                            
        if len(user_q) > 0 and verify:
            print(f" ")
            print(f"The Interpreted Queues are:") #{','.join(user_q)}")
            for q in user_q:
                print(f"{q}")
            
            message = "Is that correct? Y/N "
            tmp = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"]) 

            if tmp == "Y" or tmp == 'y':      correct = True
            elif tmp == "N" or tmp == 'n':    correct = False
            else:                             correct = False

        if correct:
            for uq in user_q:
                found = False
                for mq in self.mqueues:
                    if uq == mq.name or uq == mq.alter_name and not found: 
                        found = True
                        mq.set_nodes()
                        self.available_queues.append(mq)
                    elif uq in mq.name or mq.name in uq and not found:
                        message = f"I found a similar queue as the one you requested. Is {mq.name} your selection: {uq}? Y/N "
                        tmp = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"])
                        if tmp == "Y" or tmp == 'y':
                            found = True
                            mq.set_nodes()
                            mq.alter_name = uq
                            self.available_queues.append(mq)
                if not found: print(f"USER_QUEUES: could not find {uq} in the system queues")
        else:
            return self.available_queues

        ####################
        ## Queue Priority ##
        ####################
        if correct and len(self.available_queues) > 0:
            print(" ")
            print(f"---------------------------------------------------------------------------------------- ")
            print("Do you want to set manual priorities for the queues?                                      ")
            print("If so, the value you introduce will weight the ratio of free-CPU/total-CPU among queues   ")
            print("Otherwise, all queues will be valued equally, and will be chosen based on the above ratio ")
            print(f"---------------------------------------------------------------------------------------- ")
            print(" ")

            message = " Y/N "
            prio = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"]) 

            if prio == "Y" or prio == "y":
                print("Setting Priorities. Please type integer or float")
                for q in self.available_queues:
                    message = f"Set priority for queue={q.name}: "
                    q_prio = read_user_input(message=message, rtype=True, rtype_options=[int, float], debug=0) 
                    q.set_priority(prio=q_prio)
            elif prio == "N" or prio == "n":
                for q in self.available_queues:
                    q.set_priority(prio=int(1))
            return self.available_queues

#########################
###  Hardcoded Paths  ###
#########################
    def set_PP_Library(self, debug: int=0):
        if   os.path.isdir(self.scope_home_path+"PP_Library"):       self.PP_Library= self.scope_home_path+"PP_Library/"
        elif os.path.isdir(self.scope_scratch_path+"PP_Library"):    self.PP_Library= self.scope_scratch_path+"PP_Library/"
        else:                                                        self.PP_Library= str(input("Please Specify PP_Library Path: "))
        print("PP_Library set to", self.PP_Library)
        return self.PP_Library

    def set_paths(self, debug: int=0):
        if 'login' in self.cluster or 'csuc' in self.cluster:
            self.scope_home_path    = f"/home/{self.user}/SCOPE/Database_SCO/"
            self.scope_scratch_path = f"/scratch/{self.user}/SCOPE/Database_SCO/"
        elif 'lemma' in self.cluster:
            self.scope_home_path    = f"/Users/{self.user}/Documents/SCOPE/Database_SCO/"
            self.scope_scratch_path = None
        elif 'uam' in self.cluster:
            self.scope_home_path    = f"/home/proyectos/{self.group}/SCOPE/"
            self.scope_scratch_path = f"/scratch/{self.group}/"   #does not exist because they use ub100 instead of ub100435
        elif 'portal' in self.cluster or 'node' in self.cluster or 'visual' in self.cluster:
            if 'g2vela' in self.user:
                self.scope_home_path    = f"/home/{self.user}/SCOPE/Database_SCO/"
                self.scope_scratch_path = f"/scratch/{self.user}/SCOPE/Database_SCO/"
            elif 'g4vela' in self.user:
                self.scope_home_path    = f"/home/{self.user}/SCOPE/Database_SCO/"
                self.scope_scratch_path = f"/scratch/{self.user}/SCOPE/Database_SCO/"
        else:
            print(f"Cluster {self.cluster} not recognized")

        keyword = "4-Merged/"
        if   self.scope_home_path is not None and    os.path.isdir(self.scope_home_path+keyword):    self.cell2mol_path = self.scope_home_path+keyword
        elif self.scope_scratch_path is not None and os.path.isdir(self.scope_scratch_path+keyword): self.cell2mol_path = self.scope_scratch_path+keyword
        else:                                                self.cell2mol_path = str(input("Please Specify Cell2mol Path: "))

        keyword = "5-Complexes_Iso/"
        if   self.scope_home_path is not None and    os.path.isdir(self.scope_home_path+keyword):    self.calcs_path = self.scope_home_path+keyword
        elif self.scope_scratch_path is not None and os.path.isdir(self.scope_scratch_path+keyword): self.calcs_path = self.scope_scratch_path+keyword
        else:                                                self.calcs_path = str(input("Please Specify Calcs Path: "))

        keyword = "6-Systems_V3/"
        if   self.scope_home_path is not None and    os.path.isdir(self.scope_home_path+keyword):    self.systems_path = self.scope_home_path+keyword
        elif self.scope_scratch_path is not None and os.path.isdir(self.scope_scratch_path+keyword): self.systems_path = self.scope_scratch_path+keyword
        else:                                                self.systems_path = str(input("Please Specify Systems Path: "))

        if self.cell2mol_path[-1]   != '/': self.cell2mol_path += '/'
        if self.calcs_path[-1]      != '/': self.calcs_path    += '/'
        if self.systems_path[-1]    != '/': self.systems_path  += '/'

########################
###  Dunder Methods  ###
########################
    def __repr__(self) -> None:
        to_print  = f'\n-----------------------------------------------------------\n'
        to_print += f' Formatted input interpretation of Environment Class Object()\n'
        to_print += f'-------------------------------------------------------------\n'
        to_print += f' Cluster               = {self.cluster}\n'
        to_print += f' User                  = {self.user}\n'
        to_print += f' Group                 = {self.group}\n'
        to_print += f'\n'
        to_print += f' Queue System          = {self.management_type}\n'
        if hasattr(self,"method"): to_print += f' Method of Queue Sel   = {self.method}\n'
        if hasattr(self,"filepath"): to_print += f' Path of saved file    = {self.filepath}\n'
        if hasattr(self,"available_queues"): 
            to_print += f' Number of available queues = {len(self.available_queues)}:\n'
            for aq in self.available_queues:
                to_print += f'    {aq.name}\n'
        if hasattr(self,"selected_queues"):  
            to_print += f' Number selected queues = {len(self.selected_queues)}:\n'
            for sq in self.selected_queues:
                to_print += f'    {sq.name}\n'

        ## Prints the options added as local environment
        if hasattr(self,"added_attr"):
            if len(self.added_attr) > 0:  
                to_print += f'---------------------------------------------------\n'
                to_print += f' LOCAL ENVIRONMENT OPTIONS ADDED:\n'
                to_print += f'-------------------------------------------------------------\n'
                string    = 'self.{:15}| {:20}| {:10}\n'
                to_print += string.format('Key', 'Data Type', 'Value')
                for key in self.added_attr.keys():
                    val = self.added_attr[key]
                    to_print += string.format(key, str(type(val)), str(val))
                to_print += f'-------------------------------------------------------------\n'
        return to_print

    def _add_attr(self, key: str, value):
        try:      attr = literal_eval(value)
        except:   attr = value
        if hasattr(self,"dct"): self.dct[key] = attr
        setattr(self, key, attr)

    def _mod_attr(self, key: str, value):           ## Same as above
        try:      attr = literal_eval(value)
        except:   attr = value
        if hasattr(self,"dct"): self.dct[key] = attr
        setattr(self, key, attr)
