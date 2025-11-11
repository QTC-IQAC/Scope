from datetime import datetime, timedelta
import subprocess
from scope.parse_general import slurm_time_to_minutes

#############
### QUEUE ###
#############
class Queue(object):
    def __init__(self, name: str, _environment: object, avail: str='up', time_limit: str=None, state: str='idle'):
        if avail == "up":  self.available  = True
        else:              self.available  = False
        self.name                = fix_partition_name(name)
        self.alter_name          = name
        self.time_limit_plain    = time_limit
        self.time_limit          = slurm_time_to_minutes(time_limit)
        self.state               = state
        self.selected            = False
        self._environment        = _environment
        self.set_commands()
        #self.set_max_mem()      # Not ready

    ######
    def set_commands(self):
        if self._environment.management_type == "slurm":
            self.command_get_user_usage      = 'squeue -o "%.9P %.50j %.12u %.2t %.12M %.5C %.3D %R" | grep '+self.name  ## The rest shouldnt be necessary
            self.command_check_queue_state   = 'sinfo -o "%n %P %C" | grep '+self.name 
            self.command_job_count           = 'squeue | grep '+self.name+' | wc -l'
            self.command_get_max_mem         = 'sinfo -o "%15N %10c %10m  %25f %10G" | grep '+self.name
        elif self._environment.management_type == "sge":
            self.command_get_user_usage      = "qstat | grep "+self.name
            self.command_check_queue_state   = "qstat -f | grep "+self.name
            self.command_job_count           = "qstat | grep "+self.name+" | wc -l"

    ######
    def set_nodes(self):
        self.nodes = []
        self.max_cpu_x_node = 0 

        if not hasattr(self,"self.command_check_queue_state"): self.set_commands()
        try:
            raw = subprocess.run(['bash', '-c', self.command_check_queue_state], capture_output=True).stdout
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
                    newnode = Node(node_name, self) 
                    self.nodes.append(newnode)
                else:
                    print("SET_NODES: Unexpected length of block list when retrieving Nodes with SLURM", blocks)

        elif self._environment.management_type == 'sge':
            for idx, line in enumerate(text):
                blocks = line.replace("@"," ").replace("/"," ").split()
                if len(blocks) == 8 or len(blocks) == 9:
                    blocks = line.replace("@"," ").replace("/"," ").split()
                    queue           = str(blocks[0])
                    node_name       = str(blocks[1])
                    total           = int(blocks[5])
                    if total > self.max_cpu_x_node: self.max_cpu_x_node = total
                    newnode = Node(node_name, self) 
                    self.nodes.append(newnode)
                else:
                    print("SET_NODES: Unexpected length of block list when retrieving Nodes with SLURM", blocks)
            
        self.num_nodes = len(self.nodes)
        self.nodes.sort(key= lambda x: x.name)
        return self.nodes

    ######
    def set_priority(self, prio: int=1):
        self.priority = prio
        return self.priority

    ######
    def select_queue(self):
        if not hasattr(self._environment,"selected_queues"): self._environment.selected_queues = []
        if self not in self._environment.selected_queues: 
            self._environment.selected_queues.append(self)
        self.selected = True
        return self.selected

    ######
    def get_queue_score(self, force_update: bool=False, method: str='weighted'):
        if not hasattr(self,"priority"):      self.set_priority(prio=int(1))   ## If user has set a priority. Then 1
        if not hasattr(self,"free"):          self.get_overall_usage(force_run=force_update)

        ## If queue has no computation capacity. Score = 0 
        if self.total == 0: self.score = float(0.0); return self.score

        ## Else. First, we account for jobs that are scheduled to run on this queue. 
        if hasattr(self,"waiting_cpus"):         tmp = self.free - self.waiting_cpus 
        else:                                    tmp = self.free 

        ## Finally, we compute the score
        if   method == "weighted": self.score = float(tmp)*float(self.priority)/float(self.total)
        elif method == "total":    self.score = tmp
        elif method == "ratio":    self.score = tmp/self.total 
        else: print("GET_QUEUE_SCORE: unknown method", method); self.score = float(1.0)
        #if self.score < float(0.0): self.score = float(0.0)
        return self.score

    ######
    def get_user_running(self, debug: int=0):
        ### Evaluating the usage by queues, has the risk that pending jobs will not appear. 
        ### Thus, this function is meant to be called at the environment level, were we can also run the "get_user_waiting" method

        self.user_running_cpus = 0
        self.user_running_jobs = 0
        try: 
            #raw = subprocess.check_output(['bash','-c', self.command_get_user_usage]) ## raises error when node is not active
            raw = subprocess.run(['bash', '-c', self.command_get_user_usage], capture_output=True).stdout
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
        except Exception as exc:
            text = []

        ## SGE clusters
        if self._environment.management_type == "sge":
            try:
                for line in text:
                    if self.name in line or self.alter_name in line:
                        blocks = line.split()
                        self.user_running_cpus  += int(blocks[8])
                        self.user_running_jobs  += 1
            except Exception as exc:
                self.user_running_cpus += int(0)
                self.user_running_jobs += int(0)
                print("QUEUE.CHECK_USER_USAGE: exception:", exc)

        ## SLURM clusters
        elif self._environment.management_type == "slurm":
            try:
                for line in text:
                    if self.name in line or self.alter_name in line:
                        blocks = line.split()
                        self.user_running_cpus      += int(blocks[5])
                        self.user_running_jobs      += 1
            except Exception as exc:
                self.user_running_cpus = int(0)
                self.user_running_jobs = int(0)
                print("QUEUE.CHECK_USER_USAGE: exception:", exc)
        return self.user_running_cpus, self.user_running_jobs

    ######
    def get_overall_usage(self, force_run: bool=False, debug: int=0):
        ### For the overall usage, we have to go node by node
        ### We do not count pending jobs. This is done at the environment level
        if not hasattr(self,"nodes"):            self.set_nodes()
        if not hasattr(self,"allocated"):        self.allocated    = 0
        if not hasattr(self,"free"):             self.free         = 0
        if not hasattr(self,"max_total"):        self.max_total    = 0
        if not hasattr(self,"max_free"):         self.max_free     = 0
        if not hasattr(self,"total"):            self.total        = 0

        ### The function only updates every 60 seconds, unless force_run
        if hasattr(self,"last_usage_check"): 
            current_time = datetime.now()
            time_gap = current_time - self.last_usage_check
            if time_gap.total_seconds() > 60: proceed = True 
            else:                             proceed = False
        else:                                 proceed = True

        if proceed or force_run: 
            for node in self.nodes:
                node.get_overall_usage()
                self.allocated   += node.allocated 
                self.free        += node.free
                self.total       += node.total
                if node.total > self.max_total: self.max_total = node.total
                if node.free > self.max_free:   self.max_free = node.free
            self.last_usage_check = datetime.now()

##########################
######## DUNDER ##########
##########################
    def __repr__(self) -> None:
        to_print  = f'-------------------------------------------------------\n'
        to_print += f' Formatted input interpretation of Queue Class Object()\n'
        to_print += f'-------------------------------------------------------\n'
        to_print += f' Name                  = {self.name}\n'
        if self.alter_name != self.name:  to_print += f' Alternative Name      = {self.alter_name}\n'
        to_print += f' Time Limit Plain      = {self.time_limit_plain}\n'
        to_print += f' Time Limit (min)      = {self.time_limit} minutes\n'
        if hasattr(self,"num_nodes"):  to_print += f' Number of Nodes       = {self.num_nodes}\n'
        if hasattr(self,"priority"):   to_print += f' Priority              = {self.priority}\n'
        if hasattr(self,"score"):      to_print += f' Score                 = {self.score}\n'
        if hasattr(self,"selected"):   to_print += f' Selected              = {self.selected}\n'
        if hasattr(self,"allocated"):  to_print += f' Allocated CPUs        = {self.allocated}\n'
        if hasattr(self,"free"):       to_print += f' Available CPUs        = {self.free}\n'
        if hasattr(self,"total"):      to_print += f' Total CPUs            = {self.total}\n'
        if hasattr(self,"waiting_cpus"): to_print += f' Waiting CPUs          = {self.waiting_cpus}\n'
        if hasattr(self,"max_total"):  to_print += f' Max CPUs              = {self.max_total}\n'
        if hasattr(self,"user_running_jobs"):  to_print += f' Num Jobs {self._environment.user}       = {self.user_running_jobs}\n'        
        if hasattr(self,"user_running_cpus"):  to_print += f' Num CPUs {self._environment.user}       = {self.user_running_cpus}\n'        
        to_print += f'---------------------------------------------------\n'
        to_print += f'\n'
        return to_print

    def __add__(self, other):
        ## Adds new nodes in other to self
        if not isinstance(other, type(self)): return self
        if self.name == other.name or self.alter_name == other.alter_name:
            if hasattr(other,"nodes"):
                if type(other.nodes) == list:
                    for on in other.nodes:
                        if on.name not in [n.name for n in self.nodes]:
                            self.nodes.append(on)
        return self

############
### NODE ###
############
class Node(object):
    def __init__(self, name: str, _queue: object): 
        self.name                = name
        self._queue              = _queue
        self.set_commands()

    def set_occupation(self, total: int, allocated: int, free: int):
        self.total               = total
        self.allocated           = allocated
        self.free                = free
        self.other               = total - allocated - free

    def set_commands(self):
        if self._queue._environment.management_type == "slurm":
            self.command_check_state = 'sinfo -o "%n %P %C" | grep '+self._queue.name+' | grep '+self.name 
        elif self._queue._environment.management_type == "sge":
            self.command_check_state = 'qstat -f | grep '+self._queue.name+' | grep '+self.name 

    def get_overall_usage(self): # , user: str='all'):
        if not hasattr(self,"command_check_state"): self.set_commands()
        #raw = subprocess.check_output(['bash','-c', self.command_check_state]) ## raises error when node is not active
        raw = subprocess.run(['bash', '-c', self.command_check_state], capture_output=True).stdout
        dec = raw.decode("utf-8")
        text = dec.rstrip().split("\n")

        if dec == '':  ## It means the node is not active
            self.set_occupation(0, 0, 0)
        else: 
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
        to_print  = f'------------------------------------------------------\n'
        to_print += f' Formatted input interpretation of NODE Class Object()\n'
        to_print += f'------------------------------------------------------\n'
        if hasattr(self._queue,"name"):   to_print += f' Queue Name            = {self._queue.name}\n'
        if hasattr(self,"name"):          to_print += f' Node Name             = {self.name}\n'
        if hasattr(self,"allocated"):     to_print += f' Allocated CPUS        = {self.allocated}\n'
        if hasattr(self,"total"):         to_print += f' Total CPUS            = {self.total}\n'
        to_print += f'------------------------------------------------------\n'
        to_print += f'\n'
        return to_print

def fix_partition_name(name):
    invalid_chars = ['*', '?', '[', ']', '{', '}', '!', '@', '#', '$', '%', '^', '&', '(', ')', '+', '=', ':', ';', '"', "'", ',', '<', '>', '/', '\\', '|', '`', '~']
    for ch in invalid_chars:
        name = name.replace(ch, '')
    return name
