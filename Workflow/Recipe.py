import sys
import copy
from copy import deepcopy
import os
import numpy as np
from datetime import datetime

from Scope.Adapted_from_cell2mol import labels2formula

from Scope.Workflow import Job
from Scope.Workflow.Job import *

####################
###### RECIPE ######
####################
class recipe(object):
    def __init__(self, subject: object, _branch: object, debug: int=0) -> None:
        self.type             = "recipe"
        self._branch          = _branch
        self.path             = _branch.path
        self.keyword          = _branch.keyword
        self.subject          = subject
        self.jobs             = []
        self.isregistered     = False
        self.isgood           = False
        self.isfinished       = False
        self.results          = dict()

    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  self.results[result.key] = result
        #if result.type == "collection":
        #    if overwrite or result.key not in self.results.keys():  self.results[result.key] = result.datas
        #elif result.type == "data":
        #    if overwrite or result.key not in self.results.keys():  self.results[result.key] = result

    def add_job(self, job_data, debug: int=0):
        new_job               = job(job_data, _recipe=self)
        if debug > 1: print("Job created")
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

    def find_job(self, job_data: object, debug: int=0):
        found_job = False
        if debug > 1: print(f"Searching Job with keyword: '{job_data.keyword}' and hierarchy '{job_data.hierarchy}'")
        for idx, jb in enumerate(self.jobs):
            if jb.keyword == job_data.keyword and jb.hierarchy== job_data.hierarchy and not found_job:
                this_job = jb
                found_job = True
                if debug > 1: print(f"Job found")
        if found_job: return found_job, this_job
        else: return found_job, None

####################
### Registration ###
####################
    def register(self, debug: int=0):
        if debug > 1: print("Registering Recipe:", self.keyword)
        allgood     = True
        allfinished = True
        if len(self.jobs) > 0:
            for idx, job in enumerate(self.jobs):
                if not job.isregistered:                 job.register(debug=debug)
                if not job.isgood:                       allgood     = False
                if not job.isfinished:                   allfinished = False
        else: 
            allgood     = False
            allfinished = False
        if allgood:                 self.isgood       = True
        if allfinished:             self.isfinished   = True
        self.isregistered = True
        #if allgood and allfinished: self.isregistered = True
        if debug > 1: print("Registered Recipe:", self.keyword, "[REG, GOOD, FIN]", self.isregistered, self.isgood, self.isfinished)

#############
### Other ###
#############
    def get_info(self):
    #def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   >>> >>> RECIPE                                  \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Crystal               = {self.subject.refcode}\n'
        if hasattr(self.subject,"spin"): to_print += f' Spin                  = {self.subject.spin}\n'
        if hasattr(self.subject,"phase"): to_print += f' Phase                 = {self.subject.phase}\n'
        to_print += f' Type of Object        = {self.subject.type}\n'
        to_print += f'---------------------------------------------------\n'
        #to_print += f' self.path             = {self.path}\n'
        #to_print += f' self.keyword          = {self.keyword}\n'
        to_print += f' Num Jobs              = {len(self.jobs)}\n'
        if len(self.jobs) > 0: 
            #self.jobs.sort(key=lambda x: (x.hierarchy, x.run_number))
            self.jobs.sort(key=lambda x: x.hierarchy)
            to_print += f' Last Job Keyword      = {self.jobs[-1].keyword}\n'
            to_print += f' Last Job Hierarchy    = {self.jobs[-1].hierarchy}\n'
        #    if self.jobs[-1].run_number > 1: to_print += f' Last Job Run Number   = {self.jobs[-1].run_number}\n'
        to_print += '----------------------------------------------------\n'
        print(to_print)
