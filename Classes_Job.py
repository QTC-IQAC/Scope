import sys
import copy
from copy import deepcopy
import os
import numpy as np

from Scope.Parse_General import search_string, read_lines_file 
from Scope.Parse_G16_outputs import G16_get_last_geom, G16_time_to_sec
from Scope.Gmol_ops import gmol_update_geom
from Scope.Other import where_in_array

from cell2mol.tmcharge_common import labels2formula
from cell2mol.elementdata import ElementData
elemdatabase = ElementData()

##################################
##### RECIPE and JOB CLASSES #####
##################################
class recipe(object):
    def __init__(self, obj: object, path: str, keyword: str, debug: int=0) -> None:
        self.object = obj
        self.path = path
        self.keyword = keyword
        self.jobs = []

    def update_registry(self, code_restriction: str="None", debug: int=0):
        last_run_number = 0
        for idx, jb in enumerate(self.jobs):
            if os.path.isfile(jb.output_path):
                if not hasattr(jb,"isregistered"):
                    if code_restriction == "None":
                        lines = read_lines_file(jb.output_path, flat=True)
                        jb.register_job(lines, print_output=False, debug=debug)
                        if debug > 0: print(f"     Update_Registry: {jb.job_code} {jb.run_number} registered")
                    if code_restriction != "None":
                        if jb.code == code_restriction:
                            lines = read_lines_file(jb.output_path, flat=True)
                            jb.register_job(lines, print_output=False, debug=debug)
                            if debug > 0: print(f"     Update_Registry: {jb.job_code} {jb.run_number} registered")
            else:
                if debug > 0: print(f"     Update_Registry: deleting job with: {jb.job_code} {jb.run_number}")
                self.jobs.remove(jb)


class job(object):
    def __init__(self, hierarchy_number: int, run_number: int, input_path: str, output_path: str, subfile_path: str, software: str, code: str='', debug: int=0) -> None:
        self.hierarchy_number = int(hierarchy_number)
        self.run_number = int(run_number)
        self.input_path = input_path
        self.output_path = output_path
        self.subfile_path = subfile_path
        self.software = software
        self.code = code
        self.requisites = []

    def add_requisite(self, requisite) -> None:
        self.requisites.append(requisite)

    def add_submission_init(self, nprocs: str='Unk', queue: str='Unk', issubmitted: bool=False) -> None:
        self.nprocs = nprocs
        self.queue = queue
        self.issubmitted = issubmitted

    def register_job(self, lines: str, print_output: bool=True, debug: int=0):
        if os.path.isfile(self.output_path): self.isfinished = True
        else: self.isfinished = False

        if not self.isfinished: print(f"    REGISTRY of Job: {self.code} is impossible, since LOG:{self.output_path} does not exist")
        else:
            ###################
            ### Gaussian 16 ###
            ###################
            if self.software == 'G16':
                line_num, found_good = search_string("Normal termination", lines, typ='all')
                line_time, found_time = search_string("Elapsed time:", lines, typ='last')
                if found_time:
                    elapsed_time_list = lines[line_time].split()[2:]
                    self.elapsed_time = G16_time_to_sec(elapsed_time_list)
                else: elapsed_time = float(0)
            else:
                found_good = False
                print(f"Registry of this software: {self.software} is not implemented. See register_jobs in Control_Jobs.py")
            ###################
    
            if found_good: self.isgood = True
            else: self.isgood = False
    
        if self.isfinished: self.isregistered = True
        if not self.isfinished: self.isregistered = False
    
        if print_output or debug > 0:
            print(f"    REGISTERED for Job: {self.code}. LOG:{self.output_path}")
            if self.issubmitted: print(f"     -> Submitted")
            if self.isfinished:  print(f"     -> Finished")
            if self.isgood:      print(f"     -> Good")
            else:                    print(f"     -> BAD !! ")
            if self.isfinished:  print(f"     -> Elapsed Time: {self.elapsed_time}")


################################
##### ASSOCIATED FUNCTIONS #####
################################
def check_recipe_requisites(recipe: object, code: str, requisites: list=[], constrains: list=[], debug: int=0) -> bool:
    requisites_fulfilled = np.zeros((len(requisites)))  ## To be correct, all must be 1
    constrains_fulfilled = np.zeros((len(constrains)))  ## To be correct, all must be 0
    for idx, job in enumerate(recipe.jobs):
        if hasattr(job,"isfinished") and hasattr(job,"isregistered"):
            if job.code in requisites and job.isfinished and job.isgood: requisites_fulfilled[where_in_array(requisites,job.code)[0]] = 1
            elif job.code in constrains and job.isfinished and job.isgood: constrains_fulfilled[where_in_array(constrains,job.code)[0]] = 1
            elif job.code == code and 'self' in constrains and job.isfinished and job.isgood: constrains_fulfilled[where_in_array(constrains,'self')[0]] = 1
    if all(a == 1 for a in requisites_fulfilled) or len(requisites) == 0:
        correct = True
        if all(b == 0 for b in constrains_fulfilled) or len(constrains) == 0: correct = True
        else: correct = False
    else: correct = False
    return correct


def find_job(recipe: object, code: str, run_number: int):
    found_job = False
    for idx, jb in enumerate(recipe.jobs):
        if jb.code == code and jb.run_number == run_number and not found_job:
            this_job = jb
            found_job = True
    if found_job: return found_job, this_job
    else: return found_job, None


def get_last_run_number(recipe: object, code: str) -> int:
    run_number = 0
    for idx, jb in enumerate(recipe.jobs):
        if jb.code == code and hasattr(jb,"isregistered"):
            if jb.isfinished and jb.run_number > run_number: run_number = jb.run_number
    return int(run_number)
