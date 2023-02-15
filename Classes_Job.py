import sys
import copy
from copy import deepcopy
import os
import numpy as np
from datetime import datetime

from Scope.Parse_General import search_string, read_lines_file 
from Scope.Parse_G16_outputs import G16_get_last_geom, G16_time_to_sec
from Scope.Parse_QE_outputs import *
from Scope.Gmol_ops import gmol_update_geom
from Scope.Other import where_in_array
from Scope.Control_Jobs import set_cluster, set_user, check_submitted_job 

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

#    def update_registry(self, code_restriction: str="None", debug: int=0):
#        last_run_number = 0
#        for idx, jb in enumerate(self.jobs):
#            if os.path.isfile(jb.output_path):
#                if not hasattr(jb,"isregistered"):
#                    if code_restriction == "None":
#                        lines = read_lines_file(jb.output_path, flat=True)
#                        jb.register_job(lines, print_output=False, debug=debug)
#                        if debug > 0: print(f"     Update_Registry: {jb.job_code} {jb.run_number} registered")
#                    if code_restriction != "None":
#                        if jb.code == code_restriction:
#                            lines = read_lines_file(jb.output_path, flat=True)
#                            jb.register_job(lines, print_output=False, debug=debug)
#                            if debug > 0: print(f"     Update_Registry: {jb.job_code} {jb.run_number} registered")
#            else:
#                if debug > 0: print(f"     Update_Registry: deleting job with: {jb.job_code} {jb.run_number}")
#                self.jobs.remove(jb)


class job(object):
    def __init__(self, name: str, hierarchy_number: int, run_number: int, input_path: str, output_path: str, subfile_path: str, software: str, code: str='', debug: int=0) -> None:

        self.name = name
        self.code = code
        self.hierarchy_number = int(hierarchy_number)
        self.run_number = int(run_number)

        self.input_path = input_path
        self.output_path = output_path
        self.subfile_path = subfile_path
        self.software = software

        self.requisites = []
        self.isregistered = False

    def add_requisite(self, requisite) -> None:
        self.requisites.append(requisite)

    def add_submission_init(self, nprocs: str='Unk', queue: str='Unk', cluster: str=set_cluster(), user: str=set_user()) -> None:
        self.input_exists   = os.path.isfile(self.input_path)
        self.subfile_exists = os.path.isfile(self.subfile_path)
        self.nprocs = nprocs
        self.queue = queue
        self.submission_cluster = cluster 
        self.submission_user = user 

    def add_submission_end(self) -> None:
        self.output_exists = os.path.isfile(self.output_path)
        if self.output_exists: self.output_lines = read_lines_file(self.output_path, flat=True)
        #else:                  self.output_lines = '' 

    def check_submission_status(self) -> None:
        if self.name == '': self.set_name()
        self.isrunning = check_submitted_job(self.name)

    def add_registration_data(self, cluster: str=set_cluster(), user: str=set_user()) -> None:
        self.registration_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        self.registration_cluster = cluster 
        self.registration_user = user 

    ## Shouldn't  be necessary in the future
    def set_name(self) -> None:
        if hasattr(self,"output_path"):
            if not hasattr(self,"name"): self.name = self.output_path.rstrip().lstrip().split("/")[-1].split(".")[0]
            if self.name == '': self.name = self.output_path.rstrip().lstrip().split("/")[-1].split(".")[0]

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
        if jb.code == code and hasattr(jb,"isfinished") and hasattr(jb,"run_number"):
            if jb.isfinished and jb.run_number > run_number: run_number = jb.run_number
    return int(run_number)
