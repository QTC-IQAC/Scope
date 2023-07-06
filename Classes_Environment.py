import sys
import os
import copy
from Scope.Classes_Queue import queue 


def set_user():
    return pwd.getpwuid( os.getuid() ).pw_name

def set_cluster():
    return os.uname()[1]

###############
### CLUSTER ###
###############
class environment(object):
    def __init__(self):
        self.type    = "environment"
        self.cluster = set_cluster()
        self.user    = set_user()

    def read_user_specs(self, to_add, debug: int=0):
        if isinstance(to_add, input_data): self += to_add
        elif type(to_add) == str: 
            environment = set_environment_data(to_add, section="&environment", debug=debug)
            self += to_add

    def check_paths(self, debug: int=0):
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
            self.scope_home_path    = f"/home/{user}/SCOPE/Database_SCO/" 
            self.scope_scratch_path = f"/scratch/{user}/SCOPE/Database_SCO/"
        elif 'lemma' in self.cluster:
            self.scope_home_path    = f"/Users/{user}/Documents/SCOPE/Database_SCO/"
            self.scope_scratch_path = None
        elif 'uam' in self.cluster:
            self.scope_home_path    = f"/home/proyectos/{user}/SCOPE/"
            self.scope_scratch_path = None
        elif 'portal' in self.cluster or 'node' in self.cluster or 'visual' in self.cluster:
            if 'g2vela' in user:
                self.scope_home_path    = f"/home/{user}/SCOPE/Database_SCO/"
                self.scope_scratch_path = f"/scratch/{user}/SCOPE/Database_SCO/"
            elif 'g4vela' in user:
                self.scope_home_path    = f"/home/{user}/SCOPE/Database_SCO/"
                self.scope_scratch_path = f"/scratch/{user}/SCOPE/Database_SCO/"
        else:
            print(f"Cluster {self.cluster} not recognized")

        keyword = "4-Merged/"
        if   os.path.isdir(self.scope_home_path+keyword):    self.cell2mol_path = self.scope_home_path+keyword
        elif os.path.isdir(self.scope_scratch_path+keyword): self.cell2mol_path = self.scope_scratch_path+keyword
        else:                                                self.cell2mol_path = str(input("Please Specify Cell2mol Path:")) 

        keyword = "5-Complexes_Iso/"
        if   os.path.isdir(self.scope_home_path+keyword):    self.calcs_path = self.scope_home_path+keyword
        elif os.path.isdir(self.scope_scratch_path+keyword): self.calcs_path = self.scope_scratch_path+keyword
        else:                                                self.calcs_path = str(input("Please Specify Calcs Path:")) 

        keyword = "6-Systems_V3/"
        if   os.path.isdir(self.scope_home_path+keyword):    self.systems_path = self.scope_home_path+keyword
        elif os.path.isdir(self.scope_scratch_path+keyword): self.systems_path = self.scope_scratch_path+keyword
        else:                                                self.systems_path = str(input("Please Specify Systems Path:")) 
        
        if self.cell2mol_path[-1]   != '/': self.cell2mol_path += '/'
        if self.calcs_path[-1]      != '/': self.calcs_path    += '/'
        if self.systems_path[-1]    != '/': self.systems_path  += '/'

    def __add__(self, other):
        if not isinstance(other, input_data): return self
        for d in dir(other):
            if '_' not in d and not callable(getattr(other,d)) and d != "dct":
                at1 = getattr(other,d)
                self._add_attr(d,at1)
        return self

    def get_management_type(self, debug: int=0):
        self.management_type = "unknown"
        
        ### Sun Grid Engine, SGE ###
        worked_sge = False
        try: 
            sge_queues = subprocess.check_output(['bash','-c', "qconf -sql"])
            worked_sge = True
        except: pass
  
        ### Slurm Manager ###
        worked_slurm = False
        try: 
            slurm_queues = subprocess.check_output(['bash','-c', "sinfo"])
            worked_slurm = True
        except: pass

        if worked_sge and not worked_slurm:   self.management_type = "sge"
        elif not worked_sge and worked_slurm: self.management_type = "slurm"
        else: 
            print("Confict with the recognition of the Queue System")
            print("SGE:", worked_sge)
            print("Slurm:", worked_slurm)
        return self.management_type

    def get_all_queues(self, debug: int=0):
        self.queues  = []
        if not hasattr(self,"management_type"): self.get_management_type()
        if self.management_type == "sge": 
            raw = subprocess.check_output(['bash','-c', "qconf -sql"]) 

        elif self.management_type == "slurm": 
            raw = subprocess.check_output(['bash','-c', "sinfo"]) 
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
            for idx, line in enumerate(text):
                if idx > 0: 
                    name, avail, time_limit, num_nodes, state, name_nodes = line.split()
                    new_queue = queue(, name, avail, time_limit, num_nodes, state, name_nodes)
                    self.queues.append(new_queue)
                
                new_queue = queue(name, time_limit, num_nodes)
                 
####################################################
    def user_queues_preferences(self, debug: int=0):
        if not hasattr(self,"management_queues"):   self.get_management_queue_info()
        if   'login'  in self.cluster or 'csuc' in self.cluster: suggested = list(["std"])
        elif 'uam'    in self.cluster:                           suggested = list(["class_a"])
        elif 'portal' in self.cluster:                        
            suggested = []
            list_of_exceptions = [1, 3, 5, 7]
            for i in range(1,11):
                if i in list_of_exceptions: pass
                elif i < 10:                suggested.append(str("iqtc0"+str(i)))
                else:                       suggested.append(str("iqtc"+str(i))) 
        else: suggested = []

        print(f"Setting Queues For Environment")
        if len(suggested) > 0: print("Suggested queues are: {','.join(suggested)}")
        
        attempt = 0
        correct = False
        while attempt < 3 and not correct:
            if len(suggested) > 0: 
                keep = input("Do you want to keep the suggested ones? Y/N ")
                if keep == "Y" or keep == 'y':  
                    self.list_q = suggested.copy()
                elif keep == "N" or keep == 'n':
                    self.list_q = []
                    tmp = input("Write the available Queues/Partitions to submit jobs as a comma-separated list ")
                    tmp = tmp.strip().split(',')
                    for t in tmp:
                        self.list_q.append(str(t.strip())) 
                elif keep != "N" and keep != "Y" and keep != 'n' and keep != 'y':
                    print("I could not understand our option. Please Try Again ")
                    
            elif len(suggested) == 0: 
                self.list_q = []
                tmp = input("Write the available Queues/Partitions to submit jobs as a comma-separated list ")
                tmp = tmp.strip().split(',')
                for t in tmp:
                    self.list_q.append(str(t.strip())) 
                            
            #if len(self.list_q) > 0:
            print(f"The Interpreted Queues are: {','.join(self.list_q)}")
            for q in self.list_q:
                print(f"{q}")
            
            understood = False
            while not understood:
                tmp = input("Is that correct? Y/N ")
                if tmp != "N" and tmp != "Y" and tmp != 'n' and tmp != 'y': understood = False
                else:                                                       understood = True
            
            if understood and (tmp == "Y" or tmp == 'y'): 
                correct = True
                return self.list_q
            elif understood and (tmp == "N" or tmp == 'n'): 
                correct = False
                self.list_q = []
