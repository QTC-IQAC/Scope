import sys
import os
import copy
import pwd
import grp 
import subprocess

from Scope.Classes_Queue import queue 
from Scope.Classes_Input import input_data
from Scope import Classes_Input
from Scope.Read_Write import read_user_input

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
        self.type        = "environment"
        self.cluster     = set_cluster()
        self.user        = set_user()
        self.group       = set_group()
        self.get_management_type() 

    def __repr__(self) -> None:
        to_print  = f'\n-----------------------------------------------------------\n'
        to_print += f' Formatted input interpretation of Environment Class Object()\n'
        to_print += f'-------------------------------------------------------------\n'

        for d in dir(self):
            if d[0] != '_' and not callable(getattr(self,d)) and d != "dct" and d != 'mqueues' and d != 'queues':
                val = getattr(self,d)
                to_print += f' {d:20} = {val} \n'

        if hasattr(self,"queues"):   
            if len(self.queues) > 0: to_print += f' Selected Queues: \n'
            for q in self.queues: 
                if hasattr(q,"selected"):
                    if q.selected: to_print += f'     Queue Name             = {q.name}\n'
        to_print += f'---------------------------------------------------\n'
        return to_print

    def save(self, filepath: str):
        self.filepath    = filepath
        from Scope.Read_Write import save_binary
        save_binary(self, filepath)

    def __add__(self, other):   ## environment-class can be enriched using the input_data-class
        if not hasattr(other,"type"): return self
        if other.type != "input_data": print("Not input data but:", type(other)); return self
        for d in dir(other):
            if '_' not in d and not callable(getattr(other,d)) and d != "dct" and d != "type":
                at1 = getattr(other,d)
                self._add_attr(d,at1)
        return self

    def _add_attr(self, key: str, value):
        try:      attr = literal_eval(value)
        except:   attr = value
        if hasattr(self,"dct"): self.dct[key] = attr
        setattr(self, key, attr)

    def read_user_specs_from_file(self, path_to_add, debug: int=0):
        from Classes_Input import set_environment_data
        new_data = set_environment_data(path_to_add, section="&environment", debug=debug)
        for d in dir(new_data):
            if '_' not in d and not callable(getattr(new_data,d)) and d != "dct" and d != "type":
                at1 = getattr(new_data,d)
                self._add_attr(d,at1)
        return self

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

    def set_PP_Library(self, debug: int=0):
        if   os.path.isdir(self.scope_home_path+"PP_Library"):       self.PP_Library= self.scope_home_path+"PP_Library/"
        elif os.path.isdir(self.scope_scratch_path+"PP_Library"):    self.PP_Library= self.scope_scratch_path+"PP_Library/"
        else:                                                        self.PP_Library= str(input("Please Specify PP_Library Path:")) 
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
        else:                                                self.cell2mol_path = str(input("Please Specify Cell2mol Path:")) 

        keyword = "5-Complexes_Iso/"
        if   self.scope_home_path is not None and    os.path.isdir(self.scope_home_path+keyword):    self.calcs_path = self.scope_home_path+keyword
        elif self.scope_scratch_path is not None and os.path.isdir(self.scope_scratch_path+keyword): self.calcs_path = self.scope_scratch_path+keyword
        else:                                                self.calcs_path = str(input("Please Specify Calcs Path:")) 

        keyword = "6-Systems_V3/"
        if   self.scope_home_path is not None and    os.path.isdir(self.scope_home_path+keyword):    self.systems_path = self.scope_home_path+keyword
        elif self.scope_scratch_path is not None and os.path.isdir(self.scope_scratch_path+keyword): self.systems_path = self.scope_scratch_path+keyword
        else:                                                self.systems_path = str(input("Please Specify Systems Path:")) 
        
        if self.cell2mol_path[-1]   != '/': self.cell2mol_path += '/'
        if self.calcs_path[-1]      != '/': self.calcs_path    += '/'
        if self.systems_path[-1]    != '/': self.systems_path  += '/'

    def get_management_type(self, debug: int=0):
        self.management_type = "unknown"
        
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
            print("Confict with the recognition of the Queue System")
            print("SGE:", worked_sge)
            print("Slurm:", worked_slurm)
        return self.management_type

    def get_mqueues(self, debug: int=0):
        self.mqueues = []

        ##  SGE  ##
        if not hasattr(self,"management_type"): self.get_management_type()
        if self.management_type == "sge": 
            raw = subprocess.check_output(['bash','-c', "qconf -sql"]) 
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

####################################################
    def user_queues_preferences(self, debug: int=0):
        if not hasattr(self,"mqueues"):   self.get_mqueues()

        suggested_names = list(q.name for q in self.mqueues if q.available)
        print(f"Setting Queues For Environment")

        if len(suggested_names) > 0: 
            print(f"-------------------------------------------------")
            print(f"         We found the following Queues:          ")
            print(f"-------------------------------------------------")
            #print(f"Suggested queues are: {','.join(suggested_names)}")
            for n in suggested_names:
                print(n)
            print(f"-------------------------------------------------")
            message = "Do you want to keep all suggested Queues? Y/N "
            print(" ")
            keep = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"]) 

            if keep == "Y" or keep == 'y':        
                #user_q = self.mqueues.copy()
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
                    if tmp != '':    user_q.append(tmp.strip())
                    else: finished = True
                    
        elif len(suggested_names) == 0: 
            verify = True
            user_q = []
            message = "Write the Queues/Partitions that will be available to SCOPE. Write them one by one: Enter to Stop  "
            finished = False
            while not finished: 
                tmp = read_user_input(message=message, rtype=True, rtype_options=[str]) 
                if tmp != '':    user_q.append(tmp.strip())
                else: finished = True
                            
        if len(user_q) > 0 and verify:
            print(f"The Interpreted Queues are:") #{','.join(user_q)}")
            for q in user_q:
                print(f"{q}")
            
            message = "Is that correct? Y/N "
            tmp = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"]) 

            if tmp == "Y" or tmp == 'y':      correct = True
            elif tmp == "N" or tmp == 'n':    correct = False
            else:                             correct = False

        self.queues = []
        if correct:
            for uq in user_q:
                found = False
                for mq in self.mqueues:
                    if uq == mq.name or uq in mq.name:
                        mq.select_queue()
                        mq.set_nodes()
                        self.queues.append(mq)
                        found = True
                if not found: print(f"USER_QUEUES: could not find {uq} in the system queues")
                
        else:
            return self.queues


        ####################
        ## Queue Priority ##
        ####################
        if correct and len(self.queues) > 0:
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
                for q in self.queues:
                    message = f"Set priority for queue={q.name}: "
                    q_prio = read_user_input(message=message, rtype=True, rtype_options=[int, float], debug=0) 
                    q.set_priority(prio=q_prio)
            elif prio == "N" or prio == "n":
                for q in self.queues:
                    q.set_priority(prio=int(1))
            return self.queues

