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
    release = os.uname()[2]
    if   release == "4.18.0-513.11.1.el8_9.x86_64": cluster = "csuc3"
    elif release == "3.10.0-693.5.2.el7.x86_6":     cluster = "csuc2"
    elif release == "4.19.0-10-amd64":              cluster = "portal"
    elif release == "4.18.0-305.3.1.el8_4.x86_64":  cluster = "cesga"
    elif release == "23.6.0":                       cluster = "duke"
    return cluster

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
        self.set_commands()

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

    def set_commands(self):
        if self.management_type == "slurm":
            string_squeue = str('"%.10i %.9P %.50j %.12u %.2t %.12M %.5C %.3D %R"')
            self.command_get_user_usage      = f"squeue -o {string_squeue}"
            self.command_get_user_waiting    = f"squeue -o {string_squeue} | grep {self.user} | grep ' PD '"
            self.command_check_job           = 'squeue -o "%.60j %.12u"'
            self.command_submit              = 'sbatch'
        elif self.management_type == "sge":
            self.command_get_user_usage      = "qstat"
            self.command_get_user_waiting    = "qstat -f | grep ' qw '"
            self.command_check_job           = "qstat -xml | grep JB_name"
            self.command_submit              = "qsub"

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
           ## It should be possible to get the queue memory or mem-per-cpu doing: "sinfo -o "%15N %10c %10m  %25f %10G""
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
            if q.name == queue_name or q.alter_name == queue_name:
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

    def save(self, filepath=None):
        if filepath is None and hasattr(self,"filepath"):            pass
        elif filepath is not None and hasattr(self,"filepath"):      self.filepath = filepath
        elif filepath is not None and not hasattr(self,"filepath"):  self.filepath = filepath
        elif filepath is None and not hasattr(self,"filepath"): 
            print("ENVIRONMENT.SAVE: please re-run and provide filepath")
            return None
        from Scope.Read_Write import save_binary
        save_binary(self, self.filepath)

#####################################
###  Connection with Execute_Job  ###
#####################################
    def get_user_requested(self, debug: int=0):
        #if not hasattr(self,"user_waiting_cpus"): self.get_user_waiting(debug=debug)
        #if not hasattr(self,"user_running_cpus"): self.get_user_running(debug=debug)
        self.get_user_waiting(debug=debug)
        self.get_user_running(debug=debug)
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
        try:
            if self.management_type == "sge":
                for line in text:
                    blocks = line.split()
                    if len(blocks) == 8:
                        if blocks[4] == 'qw':
                            self.user_waiting_cpus  += int(blocks[7])
                            self.user_waiting_jobs  += 1
            elif self.management_type == "slurm":
                for line in text:
                    blocks = line.split()
                    if len(blocks) == 9:
                        if blocks[4] == 'PD':
                            self.user_waiting_cpus  += int(blocks[6])
                            self.user_waiting_jobs  += 1
        except Exception as exc:
            self.user_waiting_cpus = int(0)
            self.user_waiting_jobs = int(0)

        return self.user_waiting_cpus, self.user_waiting_jobs

############## NOW ###
    def assign_waiting_jobs(self, debug: int=0):
        if not hasattr(self,"command_get_user_waiting"): self.set_commands()
        self.jobs_assigned = []
        self.jobs_pending  = []
        ### Retrieve all waiting jobs
        try:
            raw = subprocess.check_output(['bash','-c', self.command_get_user_waiting], stderr=subprocess.DEVNULL)
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
        except Exception as exc:
            print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: exception retrieving waiting jobs")
            print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: likely due to no jobs in queue")
            text = []

        ### Save all their job_ids
        all_job_id = []   # This will store all job_ids that are going to be parsed now

        ## SGE clusters
        if self.management_type == "sge":
            for line in text:
                blocks = line.split()
                if len(blocks) == 8:
                    if blocks[4] == 'qw':
                        job_id = int(blocks[0])
                        if debug > 0: print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: job {job_id} might go to pending")
                        all_job_id.append(job_id)
                        if job_id not in self.jobs_assigned and job_id not in self.jobs_pending: 
                            self.jobs_pending.append(job_id)
                            if debug > 0: print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: pending job: {job_id}")

        ## SLURM clusters
        elif self.management_type == "slurm":
            for line in text:
                blocks = line.split()
                if len(blocks) == 9:
                    if blocks[3] == self.user and blocks[4] == 'PD':
                        job_id = int(blocks[0])
                        if debug > 0: print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: job {job_id} might go to pending")
                        all_job_id.append(job_id)
                        if job_id not in self.jobs_assigned and job_id not in self.jobs_pending: 
                            self.jobs_pending.append(job_id)
                            if debug > 0: print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: pending job: {job_id}")

        ### Assign pending jobs
        for jp in self.jobs_pending[:]: 
            if debug > 0: print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: evaluating pending job: {jp}")

            if self.management_type == "sge":
                # Find Queue
                raw2   = subprocess.check_output(['bash','-c', f"qstat -j {jp} | grep 'hard_queue_list'"])
                dec2   = raw2.decode("utf-8")
                text2  = dec2.rstrip().split("\n")[0]
                blocks = text2.split()
                if len(blocks) == 2:
                    q_name = text2.split()[1]
                    for aq in self.available_queues: 
                        if q_name == aq.name or q_name == aq.alter_name:
                            if debug > 0: print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: pending job: {jp}. queue found: {aq.name}")
                            if not hasattr(aq,"waiting_cpus"): aq.waiting_cpus = 0
                            if not hasattr(aq,"waiting_jobs"): aq.waiting_jobs = []
                            queue_to_assign = aq
                else: print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: parsing of hard_queue_list failed. Blocks:{blocks}")

                # Find number of requested processors
                raw3   = subprocess.check_output(['bash','-c', f"qstat -j {jp} | grep 'parallel environment'"])
                dec3   = raw3.decode("utf-8")
                text3  = dec3.rstrip().split("\n")[0]
                blocks = text3.split()
                if len(blocks) == 5 and blocks[4].isdigit():
                    queue_to_assign.waiting_cpus += int(blocks[4])
                    queue_to_assign.waiting_jobs.append(jp)
                    self.jobs_assigned.append(jp) 
                    self.jobs_pending.remove(jp) 
                    if debug > 0: print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: job: {jp} assigned to {aq.name}")
                    if debug > 0: print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: queue has now {queue_to_assign.waiting_cpus} cpus waiting")
                else: print(f"ENVIRONMENT.ASSIGN_WAITING_JOBS: parsing of requested processors failed. Blocks:{blocks}")

            elif self.management_type == "slurm": 
                pass

            ### Flushes out old assigned job ids, so that one of the checks above, doesnt get too costly
            if len(self.jobs_assigned) > 0:
                min_job_id = np.min(all_job_id)
                for ja in self.jobs_assigned[:]:
                    if ja < (min_job_id + 5000): self.jobs_assigned.remove(ja)

############## NOW ###

    def get_user_running(self, method: str='direct', debug: int=0):
        if not hasattr(self,"command_get_user_waiting"): self.set_commands()
        self.user_running_cpus = 0
        self.user_running_jobs = 0
        if len(self.available_queues) == 0: 
            print("CHECK_USER_RUNNING: please set available queues first by running the method 'set_queues'") 
            return None

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
                    print("ENVIRONMENT.CHECK_USER_RUNNING: exception:", exc)
        
            ## SLURM clusters
            elif self.management_type == "slurm":
                try:
                    for line in text:
                        blocks = line.split()
                        if len(blocks) == 9:
                            if blocks[3] == self.user and blocks[4] == 'R':
                                self.user_running_cpus      += int(blocks[6])
                                self.user_running_jobs      += 1
                except Exception as exc:
                    self.user_running_cpus = int(0)
                    self.user_running_jobs = int(0)
                    print("ENVIRONMENT.CHECK_USER_RUNNING: exception:", exc)
        return self.user_running_cpus, self.user_running_jobs
         
    def get_best_queue(self, autoselect: bool=False, debug: int=0):
        ## This method gives back the best queue to submit a computation.
        ## Basically, it computes a score for each queue that is either available (self.available_queues)
        ## or selected (self.selected_queues) by the user. 
        ## The score is computed for each queue separately using:
        ## 1-the current queue availability (accounting for all users)
        ## 2-the currently pending jobs     (accounting for one user, i.e. the one using scope)
        ##
        ## 1) is computed once every 60 seconds, as it is demanding
        ## 2) is computed every time
   
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

        # Assign Waiting Jobs. Point (2) in intro
        self.assign_waiting_jobs()

        # Retrieve score. Point (1) in intro
        scores = []
        for idx, q in enumerate(target_queues):
            tmp = q.get_queue_score(method=self.method)
            scores.append(tmp)
            if debug > 0: print(f"GET_BEST_QUEUE: evaluated {q.name} with score {tmp:8.4f}")

        # Select best
        best_idx = np.argmax(scores) 
        if debug > 0: print(f"GET_BEST_QUEUE: selected best_idx={best_idx} of {scores}")
        if debug > 0: print(f"GET_BEST_QUEUE: returning {target_queues[best_idx].name}")
        return target_queues[best_idx]

    def check_submitted(self, job_name: str, debug: int=0):
        if not hasattr(self,"command_check_job"): self.set_commands()
        if not hasattr(self,"management_type"): self.get_management_type()
        if self.management_type == "None": return False 

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
        ## If it is not a computation cluster, returns an empty list
        if self.management_type == "None": return []

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
                        if mq.name not in list(aq.name for aq in self.available_queues):   ## maybe .alter_name should also be checked 
                            mq.set_nodes()
                            self.available_queues.append(mq)

                    elif uq in mq.name or mq.name in uq and not found:
                        message = f"I found a similar queue as the one you requested. Is {mq.name} your selection: {uq}? Y/N "
                        tmp = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"])
                        if tmp == "Y" or tmp == 'y':
                            found = True
                            if mq.name not in list(aq.name for aq in self.available_queues):   ## maybe .alter_name should also be checked 
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

###############
###  Paths  ###
###############
    def set_storage_path(self, debug: int=0):
        if   self.cluster == "csuc3" : self.storage_path = f"/data/{self.group}/{self.user}"
        elif self.cluster == "csuc2" : self.storage_path = f"/scratch/{self.user}/"
        elif self.cluster == "duke"  : self.storage_path = f"/scratch/{self.user}/"
        elif self.cluster == "portal": self.storage_path = f"/scratch/{self.user}/"
        else: 
            print(f"Cluster {self.cluster} not implemented")
            self.storage_path = str(input("Please Specify Path of Storage Folder (e.g. user scratch):"))
            if self.storage_path[-1] != '/': self.storage_path += '/'
        return self.storage_path

    def set_scope_main_path(self, debug: int=0):
        if   self.cluster == "csuc3" : self.scope_main_path = f"/home/{self.user}/SCOPE/Database_SCO/"
        elif self.cluster == "csuc2" : self.scope_main_path = f"/home/{self.user}/SCOPE/Database_SCO/" 
        elif self.cluster == "duke"  : self.scope_main_path = f"/Users/{self.user}/Documents/SCOPE/Database_SCO/"
        elif self.cluster == "portal": self.scope_main_path = f"/home/{self.user}/SCOPE/Database_SCO/"
        else: 
            print(f"Cluster {self.cluster} not implemented")
            self.scope_main_path = str(input("Please Specify Main Scope Folder:"))
            if self.scope_main_path[-1] != '/': self.scope_main_path += '/'
        return self.scope_main_path

    def set_PP_Library(self, debug: int=0):
        if not hasattr(self,"storage_path"): self.set_storage_path()
        if not hasattr(self,"scope_main_path"): self.set_scope_main_path()
        if   os.path.isdir(self.scope_main_path+"PP_Library"):   self.PP_Library= self.scope_main_path+"PP_Library/"
        elif os.path.isdir(self.storage_path+"PP_Library"):      self.PP_Library= self.storage_path+"PP_Library/"
        else:                                                    self.PP_Library= str(input("Please Specify PP_Library Path: ")) 
        # Corrects PP_Library path if necessary
        if self.PP_Library[-1] != '/': self.PP_Library += '/'
        print("PP_Library set to", self.PP_Library)
        return self.PP_Library

    def set_paths(self, debug: int=0):
        if not hasattr(self,"storage_path"): self.set_storage_path()
        if not hasattr(self,"scope_main_path"): self.set_scope_main_path()

        keyword = "4-Merged/"
        if   self.scope_main_path is not None and os.path.isdir(self.scope_main_path+keyword): self.cell2mol_path = self.scope_main_path+keyword
        elif self.storage_path is not None and    os.path.isdir(self.storage_path+keyword):    self.cell2mol_path = self.storage_path+keyword
        else:                                                                                  self.cell2mol_path = str(input("Please Specify Cell2mol Path: "))

        keyword = "5-Complexes_Iso/"
        if   self.scope_main_path is not None and os.path.isdir(self.scope_main_path+keyword): self.calcs_path = self.scope_main_path+keyword
        elif self.storage_path is not None and    os.path.isdir(self.storage_path+keyword):    self.calcs_path = self.storage_path+keyword
        else:                                                                                  self.calcs_path = str(input("Please Specify Calcs Path: "))

        keyword = "6-Systems_V3/"
        if   self.scope_main_path is not None and os.path.isdir(self.scope_main_path+keyword): self.sys_path = self.scope_main_path+keyword
        elif self.storage_path is not None and    os.path.isdir(self.storage_path+keyword):    self.sys_path = self.storage_path+keyword
        else:                                                                                  self.sys_path = str(input("Please Specify Systems Path: "))

        if self.cell2mol_path[-1]   != '/': self.cell2mol_path += '/'
        if self.calcs_path[-1]      != '/': self.calcs_path    += '/'
        if self.sys_path[-1]    != '/': self.sys_path  += '/'

    def check_paths(self, debug: int=0):
        if not hasattr(self,"cell2mol_path"): self.set_paths()
        if os.path.isdir(self.cell2mol_path): self.iscell2mol_path   = True 
        else:                                 self.iscell2mol_path   = False
        if os.path.exists(self.calcs_path):   self.iscalcs_path      = True
        else:                                 self.iscalcs_path      = False
        if os.path.exists(self.sys_path):     self.issys_path        = True 
        else:                                 self.issys_path        = False
        if debug > 0 and not os.path.isdir(self.cell2mol_path): print(f"ENVIRONMENT.CHECK_PATHS: {self.cell2mol_path} does not exist")
        if debug > 0 and not os.path.isdir(self.sys_path):      print(f"ENVIRONMENT.CHECK_PATHS: {self.sys_path} does not exist")
        if debug > 0 and not os.path.isdir(self.calcs_path):    print(f"ENVIRONMENT.CHECK_PATHS: {self.calcs_path} does not exist")
        if self.issys_path and self.iscalcs_path and self.iscell2mol_path: return True
        else:                                                              return False

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
        if hasattr(self,"scope_main_path"):  
            to_print += f' Paths:\n'
            to_print += f'     Scope Main       = {self.scope_main_path}\n'
            to_print += f'     Cell2mol         = {self.cell2mol_path}\n'
            to_print += f'     Computations     = {self.calcs_path}\n'
            to_print += f'     Systems          = {self.sys_path}\n'
            to_print += f'     Storage          = {self.storage_path}\n'
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
