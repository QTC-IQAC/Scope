from Scope.Workflow.Job import *

####################
###### RECIPE ######
####################
class recipe(object):
    def __init__(self, name: str, source: object, _branch: object, debug: int=0) -> None:
        self.type             = "recipe"
        self._branch          = _branch
        self.path             = _branch.path
        self.name             = name
        self.source           = source
        self.jobs             = []
        self.isregistered     = False
        self.isgood           = False
        self.isfinished       = False
        self.results          = dict()

    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  self.results[result.key] = result

    #####################################
    ### Add // Remove // ------- Jobs ###
    #####################################
    def add_job(self, job_data):                               ## As opposed to add_branch or add_recipe, add_job does not need a name, but a job_data input that will include a keyword
        new_job               = job(job_data, _recipe=self)
        self.jobs.append(new_job)
        return new_job 
    
    def remove_job(self, keyword=None, hierarchy=None):
        found = False
        if keyword is None and hierarchy is None: print("Error removing job, please indicate either keyword, or hierarchy number")
        elif keyword is not None and hierarchy is not None: print("Error removing job, please indicate only keyword or hierarchy number, not BOTH")
        else:
            for idx, jb in enumerate(self.jobs):
                if keyword is not None and hierarchy is None:
                    keyword = str(keyword)
                    if jb.keyword == keyword and not found: found = True; found_idx = idx
                elif hierarchy is not None and keyword is None:
                    hierarchy = int(hierarchy)
                    if jb.hierarchy == hierarchy and not found: found = True; found_idx = idx
        if found: del self.jobs[found_idx]

    def remove_output_lines(self):
        for job in self.jobs:
            for comp in job.computations:
                comp.delete_lines()

#########################
    def find_job(self, keyword=None, hierarchy=None, job_data=None, debug: int=0):
        found_job = False

        if keyword is None and hierarchy is None and job_data is not None:
            assert hasattr(job_data,"keyword") and hasattr(job_data,"hierarchy")
            if debug > 1: print(f"Searching Job with keyword: '{job_data.keyword}' and hierarchy '{job_data.hierarchy}'")
            for idx, jb in enumerate(self.jobs):
                if jb.keyword == job_data.keyword and jb.hierarchy == job_data.hierarchy and not found_job:
                    this_job = jb
                    found_job = True
                    if debug > 1: print(f"Job found")

        elif keyword is None and hierarchy is not None and job_data is None:  
            assert type(hierarchy) == int
            if debug > 1: print(f"Searching Job with and hierarchy '{hierarchy}'")
            for idx, jb in enumerate(self.jobs):
                if jb.hierarchy == hierarchy and not found_job:
                    this_job = jb
                    found_job = True
                    if debug > 1: print(f"Job found")

        elif keyword is not None and hierarchy is None and job_data is None:  
            assert type(keyword) == str
            if debug > 1: print(f"Searching Job with and keyword '{keyword}'")
            for idx, jb in enumerate(self.jobs):
                if jb.keyword == keyword and not found_job:
                    this_job = jb
                    found_job = True
                    if debug > 1: print(f"Job found")

        elif keyword is not None and hierarchy is not None and job_data is None:  
            assert type(keyword) == str and type(hierarchy) == int
            if debug > 1: print(f"Searching Job with keyword: '{keyword}' and hierarchy '{hierarchy}'")
            for jb in self.jobs:
                if jb.keyword == keyword and jb.hierarchy == hierarchy and not found_job:
                    this_job = jb
                    found_job = True
                    if debug > 1: print(f"Job found")

        if found_job: return found_job, this_job
        else: return found_job, None

####################
### Registration ###
####################
    def register(self, debug: int=0):
        if debug > 1: print("Registering Recipe:", self.name)
        allgood     = True
        allfinished = True
        if len(self.jobs) > 0:
            for job in self.jobs:
                if not job.isregistered:                 job.register(debug=debug)
                if not job.isgood:                       allgood     = False
                if not job.isfinished:                   allfinished = False
        else: 
            allgood     = False
            allfinished = False
        if allgood:                 self.isgood       = True
        if allfinished:             self.isfinished   = True
        self.isregistered = True
        if debug > 1: print("Registered Recipe:", self.name, "[REG, GOOD, FIN]", self.isregistered, self.isgood, self.isfinished)

#############
### Other ###
#############
    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> >>> RECIPE                                  \n'
        to_print += f'---------------------------------------------------\n'
        if hasattr(self.source,"name"):      to_print += f' Source Name                 = {self.source.name}\n'
        if hasattr(self.source,"spin"):      to_print += f' Source Spin                 = {self.source.spin}\n'
        if hasattr(self.source,"phase"):     to_print += f' Source Phase                = {self.source.phase}\n'
        to_print += f' Source Type                  = {self.source.type}\n'
        to_print += f' Source sub-Type              = {self.source.subtype}\n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Recipe Name (from Source)    = {self.name}\n'
        to_print += f' Num Jobs                     = {len(self.jobs)}\n'
        if len(self.jobs) > 0: 
            self.jobs.sort(key=lambda x: x.hierarchy)
            to_print += f'\tLast Job Keyword      = {self.jobs[-1].keyword}\n'
            to_print += f'\tLast Job Hierarchy    = {self.jobs[-1].hierarchy}\n'
        to_print += '\n'
        return to_print
