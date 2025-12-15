import os
import pwd
import grp 
import subprocess
import numpy as np
import glob
import readline
from scope.classes_queue import Queue 
from scope.read_write    import read_user_input
from dataclasses         import dataclass

def set_user():
    return pwd.getpwuid( os.getuid() ).pw_name

def set_group():
    group_id = pwd.getpwnam(set_user()).pw_gid
    return grp.getgrgid(group_id).gr_name

###############
### CLUSTER ###
###############
class Environment(object):
    """
    The `environment` class controls the computational environment of SCOPE, for job submission and resource allocation.
    It should be ready to work in SGE and Slurm, although extensive testing encompassing different version has not been done.

    Attributes:
        type (str):                     Type of the environment object.
        name (str):                     Name of the environment.
        user (str):                     User name, set via `set_user()`.
        group (str):                    Group name, set via `set_group()`.
        available_queues (list):        List of available queue objects.
        selected_queues (list):         List of selected queue objects.
        method (str):                   Method for queue selection (default: 'weighted').
    Methods:
        __init__(name):                 Initializes the environment object.
        set_scheduler():                Detects and sets the job scheduler.
        set_commands():                 Sets command-line instructions for scheduler.
        read_user_queue_list():         Processes user-provided queue names and selects corresponding queues.
        read_job_specs():               Reads and applies user-specific environment settings from a local file.
        -----------------
        set_queues():                   Interactively sets available queues for the environment.
        get_mqueues():                  Retrieves and initializes queues based on scheduler.
        add_mqueue():                   Adds a new queue to the environment.
        find_queue():                   Finds a queue by name or alternate name.
        make_queue_available():         Makes a queue available for selection.
        select_queue():                 Selects a queue for job submission.
        -----------------
        save(filepath):                 Saves the environment object to a file.
        save_config(filepath):          Saves a JSON config file in the user's config dir.
        -----------------
        set_software():                 Sets software modules for Gaussian16 and Quantum Espresso.
        set_storage_path():             Sets the storage path with tab completion.
        set_scope_program():            Sets the main scope program path with tab completion.
        set_paths():                    Sets paths for sources, calculations, and systems.
        check_paths():                  Checks if specified paths exist.
        -----------------
        get_user_requested():           Gets total CPUs and jobs requested by the user.
        get_user_waiting():             Gets number of waiting CPUs and jobs for the user.
        get_user_running():             Gets number of running CPUs and jobs for the user.
        get_best_queue():               Returns the best queue for job submission based on scoring.
        check_submitted():              Checks if a job has been submitted.
        assign_waiting_jobs():          Assigns waiting jobs to queues.
        -----------------
        __repr__():                     Returns a formatted string representation of the environment object.
        _add_attr(key, value):          Adds an attribute to the environment.
        _mod_attr(key, value):          Modifies an attribute in the environment.
    Usage:
        1) Instantiate the class giving a name and follow prompts.
        2) If in a computation cluster, run self.set_queues(), to specify which queues are available
    """
    def __init__(self, name: str):
        self.type                   = "environment"
        self.name                   = name
        self.user                   = set_user()
        self.group                  = set_group()
        self.available_queues       = [] 
        self.selected_queues        = [] 
        self.method                 = 'weighted'

        self.set_scheduler() 
        self.set_commands()

    ######
    def set_scheduler(self, debug: int=0):
        """
        Determines the type of job scheduler available on the host machine.
        The method checks for the presence of Sun Grid Engine (SGE) and Slurm by attempting to
        execute their respective queue listing commands. 
        Args:
            debug (int, optional): Debug level (currently unused). Defaults to 0.
        Returns:
            str: The detected scheduler ("sge", "slurm", or "local" if none is detected).
        """
        ### Sun Grid Engine, SGE ###
        worked_sge = False
        res = run_command("qconf -sql")
        if not res.ok: worked_sge = False 
        else:          worked_sge = True

        ### Slurm ###
        worked_slurm = False
        res = run_command("sinfo")
        if not res.ok: worked_slurm = False 
        else:          worked_slurm = True

        ### Decision ###
        if worked_sge and not worked_slurm:   self.scheduler = "sge"
        elif not worked_sge and worked_slurm: self.scheduler = "slurm"
        elif worked_sge and worked_slurm:     raise ValueError("ENV.SET_SCHEDULER: Confict with the recognition of the Queue System")
        else: 
            print("ENVIRONMENT.SET_SCHEDULER: Could not recognise the Scheduler")
            print("ENVIRONMENT.SET_SCHEDULER: Assuming this is a local computer")
            print("ENVIRONMENT.SET_SCHEDULER: If this is a computation cluster with SLURM or SGE, please report bug")
            self.scheduler = "local"
        return self.scheduler

    ######
    def set_commands(self):
        """
        Sets the command-line instructions for interacting with the scheduler
        """
        if self.scheduler == "slurm":
            string_squeue = str('"%.10i %.9P %.50j %.12u %.2t %.12M %.5C %.3D %R"')
            self.commands = {
                "get_user_usage"   : f'squeue -o {string_squeue}',
                "get_user_waiting" : f'squeue -o {string_squeue} | grep {self.user} | grep " PD "',
                "check_job"        : 'squeue -o "%.60j %.12u"',
                "submit"           : 'sbatch',
                "get_queues"       : 'sinfo',}
        elif self.scheduler == "sge":
            self.commands = {
                "get_user_usage"     : 'qstat',
                "get_user_waiting"   : 'qstat -f | grep " qw "',
                "check_job"          : 'qstat -xml | grep JB_name',
                "submit"             : 'qsub',
                "get_queues"         : 'qconf -sql',}
                #"get_queues"        : 'qstat -g c',} ## this could also be interesting for the future
        else:
            self.commands = {
                "get_user_usage"     : None,
                "get_user_waiting"   : None,
                "check_job"          : None,
                "submit"             : None,
                "get_queues"         : None,}

    ######
    def test_scheduler(self, debug: int=0):
        if not hasattr(self,"scheduler"): self.set_scheduler(debug=debug) 
        if not hasattr(self,"commands"):  self.set_commands()
        if   self.scheduler == 'slurm': 
            checks = {"squeue": 'squeue --version', "sbatch": 'sbatch --help', "squeue_format": 'squeue -o "%.10i %.9P %.50j"'}
        elif self.scheduler == 'sge': 
            checks = {"qstat_version": "qstat -help","qhost": "qhost","qstat_basic": "qstat"}
        else: return None

        ## Tests generic commands
        for name, cmd in checks.items():
            res = run_command(cmd)
            if not res.ok: print(f"ENV.TEST_SCHEDULER: Check ({name}): FAILED")
            else:          print(f"ENV.TEST_SCHEDULER: Check ({name}): OK")

        ## Tests actual commands
        for name, cmd in self.commands.items():
            if name == 'submit': continue
            if name == 'get_user_waiting': continue
            res = run_command(cmd)
            if not res.ok: print(f"ENV.TEST_SCHEDULER: Check ({name}): FAILED")
            else:          print(f"ENV.TEST_SCHEDULER: Check ({name}): OK")

    ######
    def read_user_queue_list(self, line, debug: int=0):
        """
        Processes a comma-separated string of queue names provided by the user, matches them against available queues,
        and selects the corresponding queues. Handles variations in queue name formats and prompts the user to confirm
        selection if a similar queue name is found.

        Args:
            line (str): A comma-separated string containing queue names input by the user.

        Behavior:
            - Strips whitespace and splits the input string into individual queue names.
            - For each queue name, checks for an exact match or an alternate name in the available queues.
            - If a similar queue name is found (partial match), prompts the user for confirmation before selecting.
            - Updates the alternate name of the queue if the user confirms selection of a similar queue.
        """
        ## Function to digest the queue names assuming that the user could use strange formats
        list_of_user_q = line.strip().split(",")
        if debug > 0: print(f"ENV.READ_USER_QUEUE_LIST: received {list_of_user_q}")
        for user_q in list_of_user_q:
            if debug > 0: print(f"ENV.READ_USER_QUEUE_LIST: trying to make '{user_q}' available. Searching in list of env.available_queues")
            user_q = user_q.strip()
            found = False
            for q in self.available_queues:
                if q.name == user_q or q.alter_name == user_q:
                    q.select_queue()
                    found = True
                elif user_q in q.name or q.name in user_q:
                    message = f"Found Similar queue as the one you requested. Is {q.name} your selection: {user_q} Y/N "
                    tmp = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"])
                    if tmp == "Y" or tmp == 'y':
                        q.select_queue(debug=debug)
                        q.alter_name = user_q
                        found = True
            if not found: print(f"ENV.READ_USER_QUEUE_LIST: requested queue '{user_q}' was not found")

    ######
    def read_job_specs(self, content: str, isfile: bool=True, debug: int=0):
        from scope.classes_input import set_environment_data
        if not hasattr(self,"available_queues"): print("ENV.READ_JOB_SPECS: please run 'user_queue_preferences' first in the environment class"); return None
        local_env = set_environment_data(content, isfile=isfile, debug=debug)
        if debug > 0: print("ENV.READ_JOB_SPECS: reading data:", local_env)
        self.added_attr = {}
        for d in dir(local_env):
            ## General Attributes
            if d[0] != '_' and not callable(getattr(local_env,d)) and d != "dct" and d != "type" and d != "queues" and d != "queue":
                at1 = getattr(local_env,d)
                if debug > 0: print(f"ENV.READ_JOB_SPECS: adding key={d}, value={at1}")
                self._add_attr(d,at1)
                self.added_attr[d] = at1 
            ## Queues are selected
            elif d == "queues" or d == "queue":
                at1 = getattr(local_env,d)
                self.read_user_queue_list(at1, debug=debug)
        return self.added_attr

    ######
    def get_mqueues(self, debug: int=0):
        """
        Retrieves and initializes queues based on the detected cluster scheduler (SGE or SLURM).

        Args:
            debug (int, optional): Debug level (currently unused). Defaults to 0.

        Returns:
            self.mqueues(list): A list of initialized queue objects corresponding to the available queues.
        """
        self.mqueues = []

        ## Execute Command
        res = run_command(self.commands["get_queues"])
        if not res.ok: raise RuntimeError(f"ENV.GET_MQUEUES: {res.command} failed:\n{res.stderr}")
        text = res.stdout.splitlines()

        ## Parse for SGE ##
        if not hasattr(self,"scheduler"): self.set_scheduler()
        if self.scheduler == "sge": 
            for idx, line in enumerate(text):
                name = line.split()[0] 
                res2 = run_command("qstat -f | grep '{name}'")
                if not res2.ok: raise RuntimeError(f"ENV.GET_MQUEUES: {res2.command} failed:\n{res2.stderr}")
                text2 = res2.stdout.rstrip().splitlines()
                # Add queue
                new_queue = Queue(name, self)
                self.add_mqueue(new_queue) 

        ## Parse for SLURM ##
        elif self.scheduler == "slurm": 
            for idx, line in enumerate(text):
                if idx > 0 and len(line.split()) == 6: 
                    name, avail, time_limit, num_nodes, state, name_nodes = line.split()
                    new_queue = Queue(name, self, avail=avail, time_limit=time_limit, state=state)
                    self.add_mqueue(new_queue) 
                    ## Note: It should be possible to get the queue memory or mem-per-cpu doing: "sinfo -o "%15N %10c %10m  %25f %10G""
        else: print("ENV.GET_MQUEUES: SLURM or SGE were not detected. Assuming this is a local computer")
        return self.mqueues

    def add_mqueue(self, new_queue: object):
        if new_queue.name not in list(q.name for q in self.mqueues): 
            self.mqueues.append(new_queue)
        else: 
            for q in self.mqueues:
                if q.name == new_queue.name: q += new_queue

    def find_queue(self, queue_name: str):
        if not hasattr(self,"mqueues"): raise ValueError(f"Queues have not been initialized for this environment. Run self.get_mqueues()") 
        for q in self.mqueues:
            if q.name == queue_name or q.alter_name == queue_name:
                return True, q
        return False, None

    def make_queue_available(self, queue_name: str):
        found, q = self.find_queue(queue_name)
        if found and q not in self.available_queues: self.available_queues.append(q)
        if not found: print(f"ENV.MAKE_QUEUE_AVAILABLE: queue {queue_name} not found")

    def select_queue(self, queue_name: str):
        found, q = self.find_queue(queue_name)
        if found and q not in self.selected_queues: self.selected_queues.append(q)
        if not found: print(f"ENV.SELECT_QUEUE: queue '{queue_name}' not found")

    def save(self, filepath=None):
        from scope.read_write import save_binary
        if      filepath is None and hasattr(self,"filepath"):          pass
        elif    filepath is not None and hasattr(self,"filepath"):      self.filepath = filepath
        elif    filepath is not None and not hasattr(self,"filepath"):  self.filepath = filepath
        else:   self.filepath = os.path.abspath(str(f"./scope_env_{self.name}.npy"))
        save_binary(self, self.filepath)

#######################################
###  In case a JSON should be saved ###
#######################################
    def save_config(self, filepath=None):
        from platformdirs import user_config_dir
        from scope.read_write import save_json

        if      filepath is None and hasattr(self,"filepath"):          pass
        elif    filepath is not None and hasattr(self,"filepath"):      self.filepath = filepath
        elif    filepath is not None and not hasattr(self,"filepath"):  self.filepath = filepath
        else:   self.filepath = os.path.abspath(str(f"./scope_env_{self.name}.npy"))
        
        config_dir = user_config_dir("scope")
        self.config_path = os.path.join(config_dir, f"config_{self.name}.json")
        #os.makedirs(config_dir, exist_ok=True)
        config_dict = {f"scope_env_{self.name}_filepath": self.filepath}
        save_json(config_dict, self.config_path)
        return self.config_path

    def load_config(self):
        from scope.read_write import load_json
        config_dict = load_json(self.config_path)
        return config_dict

#####################################
###  Connection with Execute_Job  ###
#####################################
    def get_user_requested(self, debug: int=0):
        self.get_user_waiting(debug=debug)
        self.get_user_running(debug=debug)
        self.user_requested_cpus = self.user_waiting_cpus + self.user_running_cpus
        self.user_requested_jobs = self.user_waiting_jobs + self.user_running_jobs
        return self.user_requested_cpus, self.user_requested_jobs

    def get_user_waiting(self, debug: int=0):
        if not hasattr(self,"commands"): self.set_commands()
        self.user_waiting_cpus = int(0)
        self.user_waiting_jobs = int(0)

        res = run_command(self.commands["get_user_waiting"])
        if res.ok: text = res.stdout.rstrip().splitlines()
        else:      return self.user_waiting_cpus, self.user_waiting_jobs

        ## SGE clusters
        try:
            if self.scheduler == "sge":
                for line in text:
                    blocks = line.split()
                    if len(blocks) == 8:
                        if blocks[4] == 'qw':
                            self.user_waiting_cpus  += int(blocks[7])
                            self.user_waiting_jobs  += 1
            elif self.scheduler == "slurm":
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

    ######
    def assign_waiting_jobs(self, debug: int=0):
        if not hasattr(self,"commands"): self.set_commands()
        self.jobs_assigned = []
        self.jobs_pending  = []
        ### Retrieve all waiting jobs
        res = run_command(self.command["get_user_waiting"])
        if not res.ok: 
            print(f"ENV.ASSIGN_WAITING_JOBS: failure retrieving waiting jobs")
            print(f"ENV.ASSIGN_WAITING_JOBS: likely due to no jobs in queue")
            text = []
        else: text = res.stdout

        ### Save all their job_ids
        all_job_id = []   # This will store all job_ids that are going to be parsed now

        ## SGE clusters
        if self.scheduler == "sge":
            for line in text:
                blocks = line.split()
                if len(blocks) == 8:
                    if blocks[4] == 'qw':
                        job_id = int(blocks[0])
                        if debug > 0: print(f"ENV.ASSIGN_WAITING_JOBS: job {job_id} might go to pending")
                        all_job_id.append(job_id)
                        if job_id not in self.jobs_assigned and job_id not in self.jobs_pending: 
                            self.jobs_pending.append(job_id)
                            if debug > 0: print(f"ENV.ASSIGN_WAITING_JOBS: pending job: {job_id}")

        ## SLURM clusters
        elif self.scheduler == "slurm":
            for line in text:
                blocks = line.split()
                if len(blocks) == 9:
                    if blocks[3] == self.user and blocks[4] == 'PD':
                        job_id = int(blocks[0])
                        if debug > 0: print(f"ENV.ASSIGN_WAITING_JOBS: job {job_id} might go to pending")
                        all_job_id.append(job_id)
                        if job_id not in self.jobs_assigned and job_id not in self.jobs_pending: 
                            self.jobs_pending.append(job_id)
                            if debug > 0: print(f"ENV.ASSIGN_WAITING_JOBS: pending job: {job_id}")

        ### Assign pending jobs
        for jp in self.jobs_pending[:]: 
            if debug > 0: print(f"ENV.ASSIGN_WAITING_JOBS: evaluating pending job: {jp}")

            if self.scheduler == "sge":
                # Find Queue
                res2   = run_command(f"qstat -j {jp} | grep 'hard_queue_list'")
                if not res2.ok: raise RuntimeError(f"ENV.ASSIGN_WAITING_JOBS: {res2.command} failed:\n{res2.stderr}")
                text2  = res2.stdout.rstrip().splitlines()[0]
                blocks = text2.split()
                if len(blocks) == 2:
                    q_name = text2.split()[1]
                    for aq in self.available_queues: 
                        if q_name == aq.name or q_name == aq.alter_name:
                            if debug > 0: print(f"ENV.ASSIGN_WAITING_JOBS: pending job: {jp}. queue found: {aq.name}")
                            if not hasattr(aq,"waiting_cpus"): aq.waiting_cpus = 0
                            if not hasattr(aq,"waiting_jobs"): aq.waiting_jobs = []
                            queue_to_assign = aq
                else: print(f"ENV.ASSIGN_WAITING_JOBS: parsing of hard_queue_list failed. Blocks:{blocks}")

                # Find number of requested processors
                res3   = run_command(f"qstat -j {jp} | grep 'parallel environment'")
                if not res3.ok: raise RuntimeError(f"ENV.ASSIGN_WAITING_JOBS: {res3.command} failed:\n{res3.stderr}")
                text3  = dec3.stdout.rstrip().splitlines()[0]
                blocks = text3.split()
                if len(blocks) == 5 and blocks[4].isdigit():
                    queue_to_assign.waiting_cpus += int(blocks[4])
                    queue_to_assign.waiting_jobs.append(jp)
                    self.jobs_assigned.append(jp) 
                    self.jobs_pending.remove(jp) 
                    if debug > 0: print(f"ENV.ASSIGN_WAITING_JOBS: job: {jp} assigned to {aq.name}")
                    if debug > 0: print(f"ENV.ASSIGN_WAITING_JOBS: queue has now {queue_to_assign.waiting_cpus} cpus waiting")
                else: print(f"ENV.ASSIGN_WAITING_JOBS: parsing of requested processors failed. Blocks:{blocks}")

            elif self.scheduler == "slurm": 
                pass

            ### Flushes out old assigned job ids, so that one of the checks above, doesnt get too costly
            if len(self.jobs_assigned) > 0:
                min_job_id = np.min(all_job_id)
                for ja in self.jobs_assigned[:]:
                    if ja < (min_job_id + 5000): self.jobs_assigned.remove(ja)

    ######
    def get_user_running(self, method: str='direct', debug: int=0):
        if not hasattr(self,"commands"): self.set_commands()
        self.user_running_cpus = 0
        self.user_running_jobs = 0
        if len(self.available_queues) == 0: 
            print("ENV.CHECK_USER_RUNNING: please set available queues first by running the method 'set_queues'") 
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
            res = run_command(self.commands["get_user_usage"])
            if not res.ok: raise RuntimeError(f"ENV.CHECK_USER_RUNNING: {res.command} failed:\n{res.stderr}")
            text = res.stdout.rstrip().splitlines()

            ## SGE clusters
            if self.scheduler == "sge":
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
                    print("ENV.CHECK_USER_RUNNING: exception:", exc)
        
            ## SLURM clusters
            elif self.scheduler == "slurm":
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
                    print("ENV.CHECK_USER_RUNNING: exception:", exc)
        return self.user_running_cpus, self.user_running_jobs
         
    def get_best_queue(self, autoselect: bool=False, debug: int=0):
        """
        Determines and returns the best queue for submitting a computation based on queue availability and pending jobs.
        The method computes a score for each queue, considering:
            1. Current queue availability (across all users), updated every 60 seconds.
            2. Number of currently pending jobs (for the current user), updated every call.
        If no queues have been selected by the user, all available queues in the cluster are considered.
        The queue with the highest score is returned.
        Args:
            autoselect (bool, optional): If True, automatically selects all available queues if none are selected. Defaults to False.
            debug (int, optional): Debug level for verbose output. Defaults to 0.
        Returns:
            Queue: The queue object with the highest computed score, or None if no queues are available.
        """
   
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
        if not hasattr(self,"commands"):   self.set_commands()
        if not hasattr(self,"scheduler"):  self.set_scheduler()
        if self.scheduler == "local": return False 

        res = run_command(self.commands["check_job"])
        if not res.ok: raise RuntimeError(f"ENV.CHECK_SUBMITTED: {res.command} failed:\n{res.stderr}")

        flat = res.stdout.replace("\n", "")
        if self.scheduler == "sge": flat = flat.replace("<","").replace(">","").replace("/","").replace("JB_name","").split()
        if str(job_name) in flat: found = True
        else:                     found = False
        return found

#########################
###  User Interaction ###
#########################
    def set_queues(self, debug: int=0):
        ## If it is not a computation cluster, returns an empty list
        if self.scheduler == "local": 
            print("ENV.SET_QUEUES. Identified local computer, returning empty list of available queues")
            return self.available_queues

        if not hasattr(self,"mqueues"):   self.get_mqueues()

        suggested_names = list(q.name for q in self.mqueues if q.available)
        print(f"\tSetting Queues For Environment")

        if len(suggested_names) > 0: 
            print(f"\t-------------------------------------------------")
            print(f"\t         The following Queues were found:        ")
            print(f"\t-------------------------------------------------")
            for n in suggested_names:
                print(f"\t{n}")
            print(f"\t-------------------------------------------------")
            message = "\tDo you want to keep ALL suggested Queues? Y/N "
            print(" ")
            keep = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"]) 

            if keep == "Y" or keep == 'y':        
                user_q = list(q.name for q in self.mqueues) 
                verify  = False
                correct = True
            elif keep == "N" or keep == 'n':      
                verify = True
                user_q = []
                message = "\tWrite the Queues/Partitions that will be available to SCOPE. Write them one by one: Enter to Stop  "
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
            message = "\tWrite the Queues/Partitions that will be available to SCOPE. Write them one by one: Enter to Stop  "
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
            print(f"\tThe Interpreted Queues are:") #{','.join(user_q)}")
            for q in user_q:
                print(f"\t{q}")
            
            message = f"\tIs that correct? Y/N "
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
                        message = f"\tI found a similar queue as the one you requested. Is {mq.name} your selection: {uq}? Y/N "
                        tmp = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"])
                        if tmp == "Y" or tmp == 'y':
                            found = True
                            if mq.name not in list(aq.name for aq in self.available_queues):   ## maybe .alter_name should also be checked 
                                mq.set_nodes()
                                mq.alter_name = uq
                                self.available_queues.append(mq)
                if not found: print(f"USER_QUEUES: could not find {uq} in the system queues")
        else:
            print(f"USER_QUEUES: You said the interpreted queues are not correct. Please restart function to repeat the process")
            return self.available_queues

        ####################
        ## Queue Priority ##
        ####################
        if correct and len(self.available_queues) > 0:
            print("")
            print(f"\t---------------------------------------------------------------------------------------- ")
            print(f"\tDo you want to set manual priorities for the queues?                                      ")
            print(f"\tYES: the value you introduce will weight the ratio of free-CPU/total-CPU among queues")
            print(f"\t     so Occupancy will still determine the best queue")
            print(f"\tNO:  all queues will be valued equally, and will be chosen based on the above ratio ")
            print(f"\t---------------------------------------------------------------------------------------- ")
            print(" ")

            message = "\tY/N "
            prio = read_user_input(message=message, rtext=True, rtext_options=["Y", "N", "y", "n"]) 

            if prio == "Y" or prio == "y":
                print("\tSetting Priorities. Please type integer or float")
                for q in self.available_queues:
                    message = f"\tSet priority for queue={q.name}: "
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
        from scope.read_write import complete_path
        # Configure readline to use tab completion
        readline.set_completer_delims(' \t\n;')
        readline.parse_and_bind("tab: complete")
        readline.set_completer(complete_path)
        self.storage_path = os.path.abspath(str(input("\tPlease specify path of storage folder (with autocomplete): ")))
        if not self.storage_path.endswith("/"):
            self.storage_path += "/"
        return self.storage_path

    #def set_storage_path(self, debug: int=0):
    #    self.storage_path = str(input("Please Specify Path of Storage Folder (e.g. user scratch):"))
    #    if self.storage_path[-1] != '/': self.storage_path += '/'
    #    return self.storage_path

    def set_scope_program(self, debug: int=0):
        from scope.read_write import complete_path
        readline.set_completer_delims(' \t\n;')
        readline.parse_and_bind("tab: complete")
        readline.set_completer(complete_path)
        self.scope_program = os.path.abspath(str(input("\tPlease Specify Main scope Folder (with autocomplete):")))
        if self.scope_program[-1] != '/': self.scope_program += '/'
        return self.scope_program

    def set_paths(self, create_folders: bool=True, debug: int=0):
        from scope.read_write import complete_path
        print("\t--------------------------------------------------------------------------------------------------------------")
        print("\tSCOPE connects a list of sources (molecules/cells), with their computations, and analyses")
        print("\tThe data is stored in system files. Please define the GENERAL paths where these 3 elements will be stored.")
        print("\t")
        print("\t                          SOURCE <--> COMPUTATION <--> SYSTEM")
        print("\t")
        print("\t !! Notice that each SYSTEM will have its own subfolder inside those paths. !!")
        print("\t--------------------------------------------------------------------------------------------------------------")

        readline.set_completer_delims(' \t\n;')
        readline.parse_and_bind("tab: complete")
        readline.set_completer(complete_path)
        self.sources_path       = os.path.abspath(str(input("\tPlease Specify Sources Path (with autocomplete): ")).strip())
        self.systems_path       = os.path.abspath(str(input("\tPlease Specify Systems Path (with autocomplete): ")).strip())
        self.computations_path  = os.path.abspath(str(input("\tPlease Specify Computations Path (with autocomplete): ")).strip())
        if self.sources_path[-1]       != '/': self.sources_path      += '/'
        if self.systems_path[-1]       != '/': self.systems_path      += '/'
        if self.computations_path[-1]  != '/': self.computations_path += '/'

        ## Create Folders if necessary:
        if create_folders and debug > 0: print("\tFolders will be created if necessary")
        if not os.path.isdir(self.sources_path)      and create_folders: os.makedirs(self.sources_path, exist_ok=True)
        if not os.path.isdir(self.systems_path)      and create_folders: os.makedirs(self.systems_path, exist_ok=True)
        if not os.path.isdir(self.computations_path) and create_folders: os.makedirs(self.computations_path, exist_ok=True)

        print("\t--------------------------------------------------------------------------------------------------------------")
        print("\tAdditionally, you can specify: (1) a storage (scratch or data) folder --> Run environment.set_storage_path()")
        print("\t                               (2) the folder where the program is    --> Run environment.set_scope_program()")
        print("\t--------------------------------------------------------------------------------------------------------------")
        print("")

    def check_paths(self, debug: int=0):
        if not hasattr(self,"sources_path"):  self.set_paths()
        self.issources_path       = os.path.isdir(self.sources_path)
        self.issystems_path       = os.path.isdir(self.systems_path)
        self.iscomputations_path  = os.path.isdir(self.computations_path)
        if not os.path.isdir(self.sources_path):      print(f"ENVIRONMENT.CHECK_PATHS: {self.sources_path} does not exist")
        if not os.path.isdir(self.systems_path):      print(f"ENVIRONMENT.CHECK_PATHS: {self.systems_path} does not exist")
        if not os.path.isdir(self.computations_path): print(f"ENVIRONMENT.CHECK_PATHS: {self.computations_path} does not exist")
        if self.issystems_path and self.iscomputations_path and self.issources_path:  return True
        else:                                                                         return False

###############
### Project ###
###############
    def get_all_systems(self):
        from scope.read_write import load_binary
        systems = []
        print(f"ENV.GET_ALL_SYSTEMS: loading all systems from {self.systems_path}")
        for sys_folder in sorted(os.listdir(self.systems_path)):
            if os.path.isdir(f"{self.systems_path}{sys_folder}") and not sys_folder.startswith("."):
                for file in sorted(os.listdir(f"{self.systems_path}{sys_folder}")):
                    if file.endswith(".npy"):
                        file_path = f"{self.systems_path}{sys_folder}/{file}" 
                        try:
                            sys = load_binary(file_path)
                            systems.append(sys)
                        except Exception as exc:
                            print(exc)
                            #pass
        return systems

################
### Software ###
################
    def set_software(self):
        if self.scheduler != "local": 
            print("\t-------------------------------------------------------------------------------------")
            print("\tSCOPE expects computations to be run with either Gaussian16 or Quantum Espresso")
            print("\tPlease introduce the modules that should be called for these two codes")
            print("\tAlternatively, modify the functions gen_QE_subfile and gen_G16_subfile to your liking")
            print("\t-------------------------------------------------------------------------------------")
            message = "\tPlease, introduce the module to run GAUSSIAN16 in this cluster (Skip if G16 is not available): "
            self.g16_module = read_user_input(message=message, rtext=False)
            message = "\tNow introduce the module to run QUANTUM ESPRESSO in this cluster (Skip if QE is not available): "
            self.qe_module = read_user_input(message=message, rtext=False)
            print("\t-------------------------------------------------------------------------------------")
            print("")

########################
###  Dunder Methods  ###
########################
    def __repr__(self) -> None:
        to_print  = f'\n'
        to_print += f'\t-------------------------------------------------------\n'
        to_print += f'\t                  SCOPE Environment \n'
        to_print += f'\t-------------------------------------------------------\n'
        to_print += f'\t User                  = {self.user}\n'
        to_print += f'\t Group                 = {self.group}\n'
        to_print += f'\n'
        if hasattr(self,"sources_path"):  
            to_print += f'\tPaths:\n'
            to_print += f'\t    Sources          = {self.sources_path}\n'
            to_print += f'\t    Systems          = {self.systems_path}\n'
            to_print += f'\t    Computations     = {self.computations_path}\n'
        if hasattr(self,"storage_path"):  to_print += f'\t    Storage Path     = {self.storage_path}\n'
        if hasattr(self,"scope_program"): to_print += f'\t    Scope Program    = {self.scope_program}\n'
        to_print += f'\n'
        if hasattr(self,"qe_module") or hasattr(self,"g16_module"): 
            to_print += f'\tAvailable Software:\n'
        if hasattr(self,"g16_module"): to_print += f'\t Module of G16         = {self.g16_module}\n'
        if hasattr(self,"qe_module"):  to_print += f'\t Module of QE          = {self.qe_module}\n'
        to_print += f'\n'
        to_print += f'\tScheduler             = {self.scheduler}\n'
        if hasattr(self,"method"):     to_print += f'\t Method of Queue Sel   = {self.method}\n'
        if hasattr(self,"filepath"):   to_print += f'\t Path of saved file    = {self.filepath}\n'
        if hasattr(self,"available_queues"): 
            to_print += f'\t Number of available queues = {len(self.available_queues)}:\n'
            for idx, aq in enumerate(self.available_queues):
                to_print += f'\t  {idx+1}: Name: {aq.name} Time_limit: {aq.time_limit} minutes #Nodes: {len(aq.nodes)} Max_CPU_per_node: {aq.max_cpu_x_node}\n'

        ## Prints the options added as local environment
        if hasattr(self,"added_attr"):
            if len(self.added_attr) > 0:  
                to_print += f'\n'
                to_print += f'\t---------------------------------------------------\n'
                to_print += f'\t LOCAL ENVIRONMENT OPTIONS ADDED:\n'
                to_print += f'\t---------------------------------------------------\n'
                string    = '\t{:15}| {:20}| {:10}\n'
                to_print += string.format('Key', 'Data Type', 'Value')
                for key in self.added_attr.keys():
                    val = self.added_attr[key]
                    to_print += f"{string.format(key, str(type(val)), str(val))}"
                to_print += f'\t---------------------------------------------------\n'
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


##############
## Commands ##
##############

####
@dataclass
class CommandResult: 
    command: str
    returncode: int
    stdout: str
    stderr: str

    @property
    def ok(self) -> bool:
        return self.returncode == 0

####
def run_command(cmd: str, timeout: int = 10) -> CommandResult:
    try:
        completed = subprocess.run(cmd,shell=True,capture_output=True,text=True,timeout=timeout)
        return CommandResult(command=cmd,returncode=completed.returncode,stdout=completed.stdout.strip(),stderr=completed.stderr.strip())
    except subprocess.TimeoutExpired:
        return CommandResult(cmd, -2, "", "Command timed out")
    except FileNotFoundError as e:
        return CommandResult(cmd, -1, "", str(e))

  
