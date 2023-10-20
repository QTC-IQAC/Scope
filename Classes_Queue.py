import subprocess
import numpy as np

from Scope.Parse_General import slurm_time_to_seconds

#############
### QUEUE ###
#############
class queue(object):
    def __init__(self, name: str, _environment: object, avail: str='up', time_limit: str='infinite', state: str='idle'):
        if avail == "up":    self.available  = True
        else:                self.available  = False
        self.name                = name
        self.time_limit_plain    = time_limit
        self.time_limit          = slurm_time_to_seconds(time_limit)
        self.state               = state
        self.selected            = False
        self._environment        = _environment

        self.set_commands()

    def set_nodes(self):
        self.nodes = []
        self.max_cpu_x_node = 0 

        if not hasattr(self,"command_check_nodes_state"): self.set_commands()
        try:
            raw = subprocess.check_output(['bash','-c', self.command_check_queue_state])
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
        except: text = ""

        if self._environment.management_type == 'slurm':
            for line in text:
                blocks = line.replace('/',' ').split()
                if len(blocks) == 6:
                    node_name = str(blocks[0])
                    queue     = str(blocks[1])
                    total     = int(blocks[5])
                    if total > self.max_cpu_x_node: self.max_cpu_x_node = total
                    newnode = node(node_name, self) 
                    self.nodes.append(newnode)
                else:
                    print("SET_NODES: Unexpected length of block list", blocks)

        elif self._environment.management_type == 'sge':
            for idx, line in enumerate(text):
                blocks = line.replace("@"," ").replace("/"," ").split()
                if len(blocks) == 8 or len(blocks) == 9:
                    blocks = line.replace("@"," ").replace("/"," ").split()
                    queue           = str(blocks[0])
                    node_name       = str(blocks[1])
                    total           = int(blocks[5])
                    if total > self.max_cpu_x_node: self.max_cpu_x_node = total
                    newnode = node(node_name, self) 
                    self.nodes.append(newnode)
                else:
                    print("SET_NODES: Unexpected length of block list", blocks)
            
        self.num_nodes = len(self.nodes)
        self.nodes.sort(key= lambda x: x.name)
        return self.nodes

    def set_priority(self, prio: int=1):
        self.priority = prio
        return self.priority

    def select_queue(self):
        self.selected = True
        return self.selected

    def set_commands(self):
        if self._environment.management_type == "slurm":
            self.command_check_user_usage    = 'squeue -o "%.9P %.50j %.12u %.2t %.12M %.5C %.3D %R" | grep '+self.name  ## The rest shouldnt be necessary
            self.command_check_queue_state   = 'sinfo -o "%N %P %C" | grep '+self.name 
            self.command_job_count           = 'squeue | grep '+self.name+' | wc -l'
        elif self._environment.management_type == "sge":
            self.command_check_user_usage    = "qstat | grep "+self.name
            self.command_check_queue_state   = "qstat -f | grep "+self.name
            self.command_job_count           = "qstat | grep "+self.name+" | wc -l"

    def get_queue_score(self, method: str='weighted'):
        if not hasattr(self,"allocated"): self.check_usage(user='all')
        if not hasattr(self,"priority"):  self.set_priority()
        if self.total_cpus == 0: self.score = float(0.0)
        else:                  
            if   method == "weighted": self.score = float(self.free)*float(self.priority)/float(self.total)
            elif method == "total":    self.score = self.free
            elif method == "ratio":    self.score = self.free/self.total 
        return self.score

    def check_user_usage(self):
        if not hasattr(self, "nodes"): self.set_nodes()
        self.user_cpus = 0
        self.user_jobs = 0

        raw = subprocess.check_output(['bash','-c', self.command_check_user_usage])
        dec = raw.decode("utf-8")
        text = dec.rstrip().split("\n")
        ## SGE clusters
        if self._environment.management_type == "sge":
            try:
                for line in text:
                    if self.name in line:
                        blocks = line.split()
                        self.user_cpus  += int(blocks[8])
                        self.user_jobs  += 1
            except Exception as exc:
                self.user_cpus = int(0)
                self.user_jobs = int(0)
                print("QUEUE.CHECK_USER_USAGE: exception:", exc)

        ## SLURM clusters
        elif self._environment.management_type == "slurm":
            try:
                for line in text:
                    blocks = line.split()
                    self.user_cpus      += int(blocks[5])
                    self.user_jobs      += 1
            except Exception as exc:
                self.user_cpus = int(0)
                self.user_jobs = int(0)
                print("QUEUE.CHECK_USER_USAGE: exception:", exc)

        return self.user_cpus, self.user_jobs

    def check_overall_usage(self, debug: int=0):
        if not hasattr(self, "nodes"): self.set_nodes()
        self.allocated    = 0
        self.free         = 0
        self.max_free     = 0
        self.total        = 0
        for node in self.nodes:
            node.check_overall_usage()
            self.allocated   += node.allocated 
            self.free        += node.free
            self.total       += node.total
            if node.free > self.max_free: self.max_free = node.free

##### MOVED TO ENVIRONMENT LEVEL

#    def check_submitted(self, job_name=None, job_id=None, debug: int=0):
#        if job_name is not None and job_id is None:
#            raw  = subprocess.check_output(['bash','-c', self.command_check_job])
#            dec  = raw.decode("utf-8")
#            flat = dec.replace("\n", "")
#            if name in flat: found = True
#            else:            found = False
#        elif job_name is None and job_id is not None:
#            print("QUEUE.CHECK_SUBMITTED: this is not yet implemented")
#            return None
#        elif job_name is None and job_id is None:
#            print("QUEUE.CHECK_SUBMITTED: I need or job_name or job_id to find job")
#            return None

    #    ######################
    #    ### NUMBER of JOBS ###        
    #    ######################
    #    try: 
    #        raw = subprocess.check_output(['bash','-c', self.command_job_count])
    #        dec = raw.decode("utf-8")
    #        text = dec.rstrip().split("\n")
    #    except: 
    #        text = '0'
    #    self.num_jobs  = int(text[0])

##########################
######## DUNDER ##########
##########################
    def __repr__(self) -> None:
        to_print  = f'\n-------------------------------------------------------\n'
        to_print += f' Formatted input interpretation of Queue Class Object()\n'
        to_print += f'-------------------------------------------------------\n'
        to_print += f' Name                  = {self.name}\n'
        to_print += f' Time Limit Plain      = {self.time_limit_plain}\n'
        to_print += f' Time Limit (s)        = {self.time_limit} seconds\n'
        if hasattr(self,"num_nodes"):  to_print += f' Number of Nodes       = {self.num_nodes}\n'
        if hasattr(self,"priority"):   to_print += f' Priority              = {self.priority}\n'
        if hasattr(self,"selected"):   to_print += f' Selected              = {self.selected}\n'
        if hasattr(self,"allocated"):  to_print += f' Allocated CPUs        = {self.allocated}\n'
        if hasattr(self,"free"):  to_print += f' Available CPUs        = {self.free}\n'
        if hasattr(self,"total"):      to_print += f' Total CPUs            = {self.total}\n'
        #if hasattr(self,"num_jobs"):   to_print += f' Num Jobs              = {self.num_jobs}\n'        ### Removed
        to_print += f'---------------------------------------------------\n'
        return to_print

    def __add__(self, other):
        if not isinstance(other, type(self)): return self
        if self.name == other.name and self.free and other.free:
            if hasattr(other,"nodes"):
                if type(other.nodes) == list:
                    for on in other.nodes:
                        if on.name not in [n.name for n in self.nodes]:
                            self.nodes.append(on)
            #    else: print(f"QUEUE: {other.name} does not have a list of nodes: {type(other.nodes)}")
            #else: print(f"QUEUE: {other.name} has no variable other.nodes")
        return self


############
### NODE ###
############
class node(object):
    def __init__(self, name: str, _queue: object): 
        self.name                = name
        self._queue              = _queue
        self.set_commands()

    def set_occupation(self, total: int, allocated: int, free: int):
        self.total               = total
        self.allocated           = allocated
        self.free           = free
        self.other               = total - allocated - free

    def set_commands(self):
        if self._queue._environment.management_type == "slurm":
            self.command_check_state = 'sinfo -o "%n %P %C" | grep '+self._queue.name+' | grep '+self.name 
        elif self._queue._environment.management_type == "sge":
            self.command_check_state = 'qstat -f | grep '+self._queue.name+' | grep '+self.name 

    #def check_user_usage(self, user: str='all'):

    def check_overall_usage(self, user: str='all'):
        if not hasattr(self,"command_check_state"): self.set_commands()
        raw = subprocess.check_output(['bash','-c', self.command_check_state])
        dec = raw.decode("utf-8")
        text = dec.rstrip().split("\n")

        if self._queue._environment.management_type == 'slurm':
            for line in text:
                blocks = line.replace('/',' ').split()
                if len(blocks) == 6:
                    total      = int(blocks[5])
                    allocated  = int(blocks[2])
                    free       = total - allocated
                    self.set_occupation(total, allocated, free)
                else:
                    print("NODE.CHECK_USAGE: Unexpected length of block list", blocks)

        elif self._queue._environment.management_type == 'sge':
            for line in text:
                blocks = line.replace("@"," ").replace("/"," ").split()
                if len(blocks) == 8 or len(blocks) == 9:
                    free       = int(blocks[4])
                    total      = int(blocks[5])
                    allocated  = total - free
                    self.set_occupation(total, allocated, free)
                else:
                    print("NODE.CHECK_USAGE: Unexpected length of block list", blocks)

    def __repr__(self) -> None:
        to_print  = f'{self._queue.name}:    {self.name}    usage: {self.allocated}/{self.total}'
        return to_print

