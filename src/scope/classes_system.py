#####################################
##### Contains the SYSTEM Class #####
#####################################
import os
from scope.read_write import save_binary
from scope.classes_workflow import Branch

class System(object):
    def __init__(self, name: str) -> None:
        self.version              = "1.0"
        self.type                 = "system"
        self.subtype              = "system"
        self.origin               = "created"
        self.name                 = name
        self.sources              = []

        ## Connection with Computational Workflokw
        self.results              = dict()
        self.branches             = []
        self.states               = []

    ######
    def __repr__(self, indirect: bool=False):
        to_print  = ''
        if not indirect: to_print += '---------------------------------\n'
        if not indirect: to_print += '   >>> SCOPE System Object >>>   \n'
        if not indirect: to_print += '---------------------------------\n'
        to_print += f' Name                  = {self.name}\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        to_print += f' Subtype               = {self.subtype}\n'
        if hasattr(self,"sources_path"):      to_print += f' Source Path           = {self.sources_path}\n'  ## Path where files with molecular or cell structures are stored
        if hasattr(self,"computations_path"): to_print += f' Computations Path     = {self.computations_path}\n'    ## Path where folders with calculations will be stored
        if hasattr(self,"system_path"):       to_print += f' System File Path      = {self.system_path}\n'      ## Path where the system object is stored
        if hasattr(self,"system_file"):       to_print += f' System File Name      = {self.system_file}\n'      ## Full system path
        if len(self.sources) > 0:
            to_print += '\n'
            to_print += f' # of Sources          = {len(self.sources)}\n'
            to_print += f'     idx: type, name, formula               \n'
            for idx, spec in enumerate(self.sources):
                to_print += f'     {idx}: {spec.type}, {spec.name}, {spec.formula} \n'
        if not indirect: to_print += '\n'
        return to_print

    ######
    def create_folders(self):
        """Create system and calculations folders if they do not exist."""
        for path in [self.system_path, self.computations_path]:
            try:
                os.makedirs(path, exist_ok=True)
            except Exception as exc:
                print(f"Error creating folder: {path}")
                print(exc)

    ######
    def save(self, filepath: str=None):
        if filepath is None: filepath = self.system_file
        save_binary(self, filepath)

    ######
    def find_source(self, name: str, debug: int=0):
        if debug > 0: 
            print(f"FIND_SOURCE. Searching for source with {name} in system with {len(self.sources)} sources:")
            for sour in self.sources:
                print(f"    {sour.name} {sour.type}")
        if len(self.sources) == 0: return False, None
        for sour in self.sources:
            if sour.name.lower() == name.lower(): return True, sour 
        return False, None

    ######
    def add_source(self, name: str, new_source: object, overwrite: bool=False, debug: int=0):
        ## Sources Must Have a Name
        if not hasattr(new_source,"name"): new_source.name = name
        ## Search if source with the same name already exists
        found, old_source = self.find_source(name, debug=debug)
        ## If not, it is added
        if not found: 
            new_source._sys = self ## Links the system to the source 
            self.sources.append(new_source)
        ## If it exists, it is overwritten if specified 
        elif found and overwrite: 
            new_source._sys = self ## Links the system to the source 
            self.sources = [s for s in self.sources if s.name.lower() != name.lower()]
            self.sources.append(new_source)
        else: 
            print(f"ADD_SOURCE: Source with name '{new_source.name}' already exists in system") 
            print(f"If you would like to Overwrite, specify overwrite=True")
        return self.sources

    #############
    ### Paths ###
    #############
    def check_paths(self, debug: int=0) -> bool:
        if not os.path.isfile(self.system_file) or not os.path.isdir(self.system_path) or not os.path.isdir(self.computations_path) or not os.path.isdir(self.sources_path):  
            if debug > 0: 
                if not os.path.isfile(self.system_file):       print(f"SYSTEM.CHECK_PATHS: WARNING: System FILE does not exist {self.system_file=}")
                if not os.path.exists(self.system_path):       print(f"SYSTEM.CHECK_PATHS: WARNING: Systems Folder does not exist {self.system_path=}")
                if not os.path.exists(self.computations_path): print(f"SYSTEM.CHECK_PATHS: WARNING: Computations Folder does not exist {self.computations_path=}")
                if not os.path.exists(self.sources_path):      print(f"SYSTEM.CHECK_PATHS: WARNING: Sources Folder does not exist {self.sources_path=}")
            return False
        return True

    ######
    def set_paths(self, create_folders: bool=True, debug: int=0) -> None: 
        from scope.read_write import complete_path
        """
        Modifies the paths associated with the system, as well as the branches, workflows, jobs, and computation files of a system
        Args:
            debug (int, optional): Debug level. Defaults to 0.
        Returns: None
        """
        import readline
        # Set up autocomplete
        readline.set_completer_delims(' \t\n;')
        readline.parse_and_bind("tab: complete")
        readline.set_completer(complete_path)

        ## Reads User Choice for Paths
        self.system_path       = os.path.abspath(str(input("\tPlease Specify Systems Path for System (with autocomplete): ")).strip())
        self.sources_path      = os.path.abspath(str(input("\tPlease Specify Sources Path for System (with autocomplete): ")).strip())
        self.computations_path = os.path.abspath(str(input("\tPlease Specify Computations Path for System (with autocomplete): ")).strip())
        if self.system_path[-1]        != '/': self.system_path      += '/'
        if self.sources_path[-1]       != '/': self.sources_path  += '/'
        if self.computations_path[-1]  != '/': self.computations_path    += '/'
        ## Sets Default system_file name
        self.system_file       = f"{self.system_path}{self.name}.npy"
        ## Create Folders if necessary:
        if not os.path.isdir(self.system_path)        and create_folders: os.makedirs(self.system_path)
        if not os.path.isdir(self.computations_path)  and create_folders: os.makedirs(self.computations_path)
        if not os.path.isdir(self.sources_path)       and create_folders: os.makedirs(self.sources_path)
        if debug > 0: 
            print(f"SYSTEM.SET_PATHS: new paths:")
            print(f"Source path: {self.sources_path}")
            print(f"Comps path:  {self.computations_path}")
            print(f"System path: {self.system_path}")
            print(f"System file: {self.system_file}")
        ## Chance all paths
        self.set_paths_down_hierarchy(debug=debug)

    ######
    def set_paths_from_environment(self, environment: object, create_folders: bool=True, debug: int=0) -> None: 
        """
        Modifies the paths associated with the system, as well as the branches, workflows, jobs, and computation files of a system
        Based on the paths stored in the environment object
        Args:
            environment (object): The environment to which the system paths must be updated
            debug (int, optional): Debug level. Defaults to 0.
        Returns: None
        """
        ## Fix for older versions
        if not hasattr(environment,"sources_path"):
            if hasattr(environment,"cell2mol_path"):  environment.sources_path = environment.cell2mol_path
            else: print("SYSTEM.SET_PATHS_FROM_ENV: WARNING: no sources_path in environment"); return False

        ## Environment Sends the Global Paths
        self.sources_path         = f"{environment.sources_path}{self.name}/"
        self.computations_path    = f"{environment.computations_path}{self.name}/"
        self.system_path          = f"{environment.system_path}{self.name}/"
        self.system_file          = f"{environment.system_path}{self.name}/{self.name}.npy"

        ## Create Folders if necessary:
        if not os.path.isdir(self.system_path)       and create_folders: os.makedirs(self.system_path)
        if not os.path.isdir(self.computations_path) and create_folders: os.makedirs(self.computations_path)
        if not os.path.isdir(self.sources_path)      and create_folders: os.makedirs(self.sources_path)
        if not create_folders and (not os.path.isdir(self.system_path) or not os.path.isdir(self.computations_path) or not os.path.isdir(self.sources_path)):
            print(f"SYSTEM.SET_PATHS_FROM_ENV: New Folders Could not be created. Please Use create_folders=True")

        if debug > 0: 
            print(f"SYSTEM.SET_PATHS_FROM_ENV: new paths:")
            print(f"Source path: {self.sources_path}")
            print(f"Comps path:  {self.computations_path}")
            print(f"System path: {self.system_path}")
            print(f"System file: {self.system_file}")
        self.set_paths_down_hierarchy(debug=debug)
        return True

    ######
    def set_paths_down_hierarchy(self, create_folders: bool=True, debug: int=0) -> None:
        """
        Modifies the paths associated with all well branches, workflows, jobs, and computation files of a system
        Based on the paths stored in this system-class object
        Args:
            debug (int, optional): Debug level. Defaults to 0.
        Returns: None
        """
        for br in self.branches:
            br.path = self.computations_path+br.keyword+'/'
            if create_folders: os.makedirs(br.path, exist_ok=True)
            if debug > 0: print(f"SET_PATHS_DOWN_HIERARCHY: new branch path: {br.path}")
            for wrk in br.workflows:
                wrk.path = br.path              ## Recipe Path is the same as branch, and is a folder
                for job in wrk.jobs:
                    job.path = wrk.path         ## Job Path is the same as branch, and is a folder. Except for finite differences jobs
                    if job.setup == "findiff": 
                        job.path = job.path+"findiff/"
                        if create_folders: os.makedirs(wrk.path, exist_ok=True)
                    for comp in job.computations:
                        comp.path = job.path                     ## Computation paths is the same as branch folder
                        comp.inp_path = comp.path+comp.inp_name  ## inp_path is a file
                        comp.out_path = comp.path+comp.out_name  ## out_path is a file
                        comp.sub_path = comp.path+comp.sub_name  ## sub_path is a file
                        comp.check_files()
                        if debug > 0: print(f"SET_PATHS_DOWN_HIERARCHY: new computation path: {comp.inp_path}")

    #########################################################
    ### Functions to Interact with Computational Workflow ###
    #########################################################
    def add_branch(self, name: str, debug: int=0):
        new_branch = Branch(self.computations_path+name, name, self, debug=debug)
        if not os.path.isdir(self.computations_path+name): 
            if debug > 0: print(f"SYSTEM.ADD_BRANCH: creating branch in {self.computations_path}{name}")
            os.makedirs(self.computations_path+name, exist_ok=True)
        self.branches.append(new_branch)
        return new_branch

    ######
    def remove_branch(self, br_name=None):
        found = False
        for idx, br in enumerate(self.branches):
            if br.name == str(br_name): found = True; found_idx = idx
        if found:
            to_delete = self.branches[found_idx]
            del self.branches[found_idx]

    ######
    def remove_all_branches(self) -> None:
        if hasattr(self,"branches"): delattr(self,"branches"); setattr(self,"branches",[])
        return self.branches

    ######
    def find_branch(self, name: str, debug: int=0):
        name = name.lower()
        if debug > 1: print(f"FIND_BRANCH. Finding branch with name:", name)
        if debug > 1: print(f"FIND_BRANCH. There are {len(self.branches)} branches in system")
        if len(self.branches) == 0: return False, None
        for idx, br in enumerate(self.branches):
            if debug > 1: print(f"FIND_BRANCH. Evaluating branch {idx} with name: {br.name} and path: {br.path}")
            if br.name.lower() == name.lower():
                if debug > 1: print(f"FIND_BRANCH. Branch was found. Checking path...")
                if not os.path.isdir(br.path) and debug > 0: 
                    print(f"WARNING: The path associated with this branch (below) does not exist. Loading the branch anyway")
                    print(f"WARNING: {br.path=}")
                return True, br
        return False, None

    ######
    def find_computation(self, branch_keyword: str, workflow_keyword: str, job_keyword: str, comp_keyword: str='', comp_step: int=1, comp_run_number: int=1, debug: int=0):
        if len(self.branches) == 0: return False, None
        for br in self.branches:
            if br.keyword.lower() == branch_keyword.lower():
                if len(br.workflows) == 0: return False, None
                for wrk in br.workflows:
                    if wrk.source.spin.lower() == workflows_keyword.lower():
                        if len(wrk.jobs) == 0: return False, None
                        for job in wrk.jobs:
                            if job.keyword.lower() == job_keyword.lower():
                                if len(job.computations) == 0: return False, None
                                for idx, comp in enumerate(job.computations):
                                    if comp.run_number == comp_run_number and comp.step == comp_step and comp.keyword.lower() == comp_keyword.lower(): return True, comp
        return False, None

    #########################################
    ### Functions to Interact with States ###
    #########################################
    def add_state(self, state: object, debug: int=0):
        ## Verifies that a state with the same name does not exist already. If not, appends it
        for my_state in self.states:
            if my_state.name == state.name: 
                found = True
                print("SYSTEM.ADD_STATE: you're trying to create a state that already exists in _source.")
                print("SYSTEM.ADD_STATE: use function called 'find_state' instead to retrieve existing state")
        if not found: self.states.append(state)

    ############################################
    ### Functions to Load Simple XYZ Sources ###
    ############################################
    def load_single_xyz(self, filepath: str, overwrite: bool=False, debug: int=0):
        ## xyz files can only be of species, not unit cells. Unit Cells require the unit cell vectors or parameters
        from scope.classes_specie import Specie
        from scope.read_write import read_xyz
        dir  = os.path.dirname(filepath)
        file = os.path.basename(filepath)
        name = file.split('.')[0]
        labels, coord = read_xyz(filepath)
        new_specie = Specie(labels, coord)
        ## Creates the Initial State
        ini_state = new_specie.add_state("initial")
        ini_state.set_geometry(new_specie.labels, new_specie.coord)
        ## Adds the Specie as a Source of System
        self.add_source(name, new_specie, overwrite=overwrite)

    def load_multiple_xyz(self, folder: str, overwrite: bool=False, debug: int=0):
        if folder[-1] != '/': folder += '/'
        if not os.path.isdir(folder): return None
        for file in sorted(os.listdir(os.path.abspath(folder))):
            if os.path.isfile(folder+file):
                self.load_single_xyz(folder+file, overwrite=overwrite, debug=debug)
        return self
