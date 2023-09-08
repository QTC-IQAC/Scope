import subprocess
import numpy as np

from Scope.Parse_General import slurm_time_to_seconds

#############
### QUEUE ###
#############
class queue(object):
    def __init__(self, name: str, avail: str, time_limit: str, num_nodes: int, state: str, name_nodes: str, _environment: object): 
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
            raw = subprocess.check_output(['bash','-c', self.command_check_nodes_state])
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
        except: text = ""

        for line in text:
            blocks = line.replace('/',' ').split()
            if len(blocks) == 6:
                name      = str(blocks[0])
                queue     = str(blocks[1])
                allocated = int(blocks[2])
                idle      = int(blocks[3])
                other     = int(blocks[4])
                total     = int(blocks[5])
                #if queue != self.name: print(queue, self.name)
                if total > self.max_cpu_x_node: self.max_cpu_x_node = total
                newnode = node(name, allocated, idle, other, total, self) 
                self.nodes.append(newnode)

        self.num_nodes = len(self.nodes)
        self.nodes.sort(key= lambda x: x.name)

        return self.nodes

    def set_priority(self, prio: int=1):
        self.priority = prio
        return self.priority

    def select_queue(self):
        self.selected = True
        return self.selected

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
        if hasattr(self,"used_cpus"):  to_print += f' Used CPUs             = {self.used_cpus}\n'
        if hasattr(self,"free_cpus"):  to_print += f' Free CPUs             = {self.free_cpus}\n' 
        if hasattr(self,"total_cpus"): to_print += f' Total CPUs            = {self.total_cpus}\n'
        if hasattr(self,"num_jobs"):   to_print += f' Num Jobs              = {self.num_jobs}\n'
        to_print += f'---------------------------------------------------\n'
        return to_print

    def __add__(self, other):
        if not isinstance(other, type(self)): return self
        if self.name == other.name and self.available and other.available: 
            if hasattr(other,"nodes"):
                if type(other.nodes) == list:
                    for on in other.nodes:
                        if on.name not in [n.name for n in self.nodes]:
                            self.nodes.append(on)
            #    else: print(f"QUEUE: {other.name} does not have a list of nodes: {type(other.nodes)}")
            #else: print(f"QUEUE: {other.name} has no variable other.nodes")
        return self
      
    def set_commands(self):
        if self._environment.management_type == "slurm":
            self.command_check_cluster_state = 'squeue -o "%.9P %.50j %.12u %.2t %.12M %.5C %.3D %R" | grep '+str(self._environment.user)
            self.command_check_queue_state   = 'sinfo -o "%N %P %C" | grep '+self.name 
            self.command_check_nodes_state   = 'sinfo -o "%n %P %C" | grep '+self.name  ## Not being used yet 
            self.command_check_job           = 'squeue -o "%.60j %.12u"' 
            self.command_job_count           = 'squeue | grep '+self.name+' | wc -l'
            self.command_submit              = 'sbatch' 
        elif self._environment.management_type == "sge":
            self.command_check_cluster_state = "qstat"
            self.command_check_queue_state   = "qstat -f | grep "+self.name
            #self.command_check_nodes_state   = ''  # Not ready
            self.command_check_job           = "qstat -xml"
            self.command_job_count           = "qstat | grep "+self.name+" | wc -l"
            self.command_submit              = "qsub"

    def get_queue_score(self):
        if not hasattr(self,"used_cpus"): self.check_usage()
        if not hasattr(self,"priority"):  self.set_priority()
        if self.total_cpus > 0: self.score = float(self.free_cpus)*float(self.priority)/float(self.total_cpus)
        else:                   self.scope = 0.0

    def check_usage(self, debug: int=0):
        self.used_cpus    = 0
        self.free_cpus    = 0
        self.total_cpus   = 0
        self.num_jobs     = 0

        ############
        ### CPUs ###        
        ############
        try: 
            raw = subprocess.check_output(['bash','-c', self.command_check_queue_state])
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
        except: text = ""

        if self._environment.management_type == "slurm":
            try:
                for line in text:
                    blocks = line.replace('/',' ').split()
                    self.used_cpus += int(blocks[2])
                    self.free_cpus += int(blocks[3])
                    self.total_cpus += int(blocks[5])
            except:
                self.used_cpus  = int(0)
                self.free_cpus  = int(0)
                self.total_cpus = int(0)
                print(f"CHECK_USAGE: Exception checking usage of queue={self.name}:", exc)

        elif self._environment.management_type == "sge":
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
            try:
                for line in text:
                    if q in line:
                        blocks = line.replace('/',' ').replace('@',' ').split()
                        self.used_cpus  += int(blocks[3])
                        self.free_cpus  += int(blocks[2])
                        self.total_cpus += int(blocks[4])
                        self.num_jobs   += 1
            except Exception as exc:
                self.used_cpus  = int(0)
                self.free_cpus  = int(0)
                self.total_cpus = int(0)
                print(f"CHECK_USAGE: Exception checking usage of queue={self.name}:", exc)

        ######################
        ### NUMBER of JOBS ###        
        ######################
        try: 
            raw = subprocess.check_output(['bash','-c', self.command_job_count])
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
        except: 
            text = '0'
        self.num_jobs  = int(text[0])

############
### NODE ###
############
class node(object):
    def __init__(self, name: str, allocated: int, idle: int, other: int, total: int, _queue: object): 
        self.name                = name
        self.allocated           = allocated
        self.idle                = idle
        self.other               = other
        self.total               = total
        self._queue              = _queue

    def __repr__(self) -> None:
        to_print  = f'{self._queue.name}:    {self.name}    usage: {self.allocated}/{self.total}'
        return to_print
