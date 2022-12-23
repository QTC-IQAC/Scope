#!/usr/bin/env python3
import sys
import os
import subprocess
import numpy as np

from Scope.Scope_Classes import job, recipe
from Scope.Other import where_in_array
from Scope.Parse_G16_outputs import G16_time_to_sec
from Scope.Parse_General import read_lines_file, search_string 

def set_queues(queues: str='all'):
    list_q = []
    if queues == 'all':
        for i in range(1,11): 
            if i == 3: pass
            elif i < 10: list_q.append(str("iqtc0"+str(i)))
            else: list_q.append(str("iqtc"+str(i)))
    else:
        digested = queues.rstrip(" ").lstrip(" ").split(',')
        for i in range(1,11): 
            if i == 3: pass
            elif i != 3 and str(i) in digested: 
                if i < 10: list_q.append(str("iqtc0"+str(i)))
                else: list_q.append(str("iqtc"+str(i)))
    return list_q     

def check_usage(queues: str='all'):
    list_q = set_queues(queues)
    cpus = 0
    jobs = 0
    raw = subprocess.check_output(['bash','-c', "qstat"])
    dec = raw.decode("utf-8")
    text = dec.rstrip().split("\n")
    for q in list_q:
        try:
            for line in text:
                if q in line:                 
                    blocks = line.split()
                    cpus += int(blocks[8])
                    jobs += 1 
        except Exception as exc:
            print("Exception checking usage:", exc)
    return cpus, jobs

def check_queue_availability(queues: str='all'):
    list_q = set_queues(queues)
    list_q_worked = []
    list_empty_cpus = []
    list_total_empty = []
    for q in list_q:
        try:
            raw = subprocess.check_output(['bash','-c', "qstat -f | grep "+q ]) 
            dec = raw.decode("utf-8") 
            text = dec.rstrip().split("\n")
            empty_cpus = []
            for node in text:
                blocks = node.split()
                queue = blocks[0].split('@')[0]
                queue_check = blocks[0].split('@')[0][0:6]
                node = blocks[0].split('@')[1]
                if queue_check == q:
                    reserved, used, total = blocks[2].split('/')
                    if queue == q+".q": empty_cpus.append(int(total)-int(reserved)-int(used)) 
                else: print("queue_check not passed for", queue_check, q)
      
            list_q_worked.append(q)
            list_empty_cpus.append(empty_cpus)
            list_total_empty.append(np.sum(empty_cpus))
        except: print(q, "does not exist")
    return list_q_worked, list_empty_cpus, list_total_empty

### DOESNT WORK BECAUSE "NAME" part of qstat is too narrow
#def check_running_job(name: str, queues: str='all', debug: int=0):
#    list_q = set_queues(queues)
#    if debug >= 1: print(list_q)
#    for q in list_q:
#        raw = subprocess.check_output(['bash','-c', "qstat -q "+q+".q" ]) 
#        dec = raw.decode("utf-8") 
#        if debug >= 1: print("for q=", q, "dec is:")
#        if debug >= 1: print(dec)
#        text = dec.rstrip().split("\n")[2::] 
#
#        list_jobID = []
#        list_jobnames = []
#        list_user = []
#        list_status = []
#        list_nodes = []
#        list_procs = []
#        for job in text:
#            blocks = job.split()
#            queue_check = blocks[7].split('@')[0][0:6]
#            list_jobID.append(int(blocks[0]))
#            list_jobnames.append(blocks[2])
#            list_user.append(blocks[3])
#            list_status.append(blocks[4])
#            list_nodes.append(blocks[7].split("@")[1])
#            list_procs.append(int(blocks[8]))
#        if name in list_jobnames: isrunning = True
#        else: isrunning = False
#    return isrunning

def check_submitted_job_xml(name: str, queues: str='all', debug: int=0):
    list_q = set_queues(queues)
    if debug >= 1: print("QUEUES:", list_q)
    for q in list_q:
        raw = subprocess.check_output(['bash','-c', "qstat -xml -q "+q+".q" ]) 
        dec = raw.decode("utf-8") 
        flat = dec.replace("\n", "")
        code = str("<JB_name>"+name+"</JB_name>") 
        if debug >= 1: print("code:", code)
        if debug >= 1: print("flat:", flat)
        if code in flat: 
            issubmitted = True 
            return issubmitted
        else: issubmitted = False
    ## Checks Pending Jobs in all queues 
    if not issubmitted:
       raw = subprocess.check_output(['bash','-c', "qstat -xml" ]) 
       dec = raw.decode("utf-8") 
       flat = dec.replace("\n", "")
       code = str("<JB_name>"+name+"</JB_name>") 
       if code in flat: 
           issubmitted = True 
           return issubmitted
    return issubmitted

def set_best_queue(queues):
    list_q, list_empty_cpus, list_total_empty = check_queue_availability(queues)

    maxempty = np.argmax(list_total_empty)
    found = False
    while not found:
        for idx, cp in enumerate(list_empty_cpus[maxempty]):
            if cp >= 8: 
                found = True
                return str(list_q[maxempty])
        list_total_empty[maxempty] = 0
        maxempty = np.argmax(list_total_empty)
        if np.sum(list_total_empty) == 0: break

def check_fairsharing(user):
    raw = subprocess.check_output(['bash','-c', "check_ocupacio" ])
    dec = raw.decode("utf-8") 
    lines = dec.split("\n")
    for l in lines:
        if "The group" in l and "has a limit for each user" in l:
            limit = l.split()[-1]
        elif user in l:
            usage = l.split()[2].split("/")[0]
    if usage < limit: return True
    else: return False

def get_computer_time(max_jobs: int=90, max_procs: int=300, queues: str='All', debug: int=0):
    sent_procs, sent_jobs = check_usage(queues)
    if sent_jobs <= max_jobs and sent_procs <= max_procs: cancontinue = True
    else: cancontinue = False

    if cancontinue:
        #### Finds best number of processors and queue
        askqueue = set_best_queue('8,9,10')
        if askqueue == 'iqtc08': askprocs = 7
        else: askprocs = 8
    else:
        askqueue = ''
        askprocs = 0
    return cancontinue, askqueue, askprocs

def check_recipe_requisites(recipe: object, code: str, requisites: list=[], constrains: list=[]) -> bool:

    requisites_fulfilled = np.zeros((len(requisites)))  ## To be able to continue, all must be 1
    constrains_fulfilled = np.zeros((len(constrains)))  ## To be able to continue, all must be 0
    for idx, job in enumerate(recipe.jobs):
        if hasattr(job,"isfinished") and hasattr(job,"isregistered"):
            if job.code in requisites and job.isfinished and job.isgood: requisites_fulfilled[where_in_array(requisites,job.code)[0]] = 1
            elif job.code in constrains and job.isfinished and job.isgood: constrains_fulfilled[where_in_array(constrains,job.code)[0]] = 1
            elif job.code == code and 'self' in constrains and job.isfinished and job.isgood: constrains_fulfilled[where_in_array(constrains,'self')[0]] = 1

    if all(a == 1 for a in requisites_fulfilled) or len(requisites) == 0: 
        can_continue = True 
        if all(b == 0 for b in constrains_fulfilled) or len(constrains) == 0: can_continue = True
        else: can_continue = False
    else: can_continue = False

    return can_continue

def find_same_job(recipe: object, code: str):
    found_job = False
    for idx, jb in enumerate(recipe.jobs):
        if jb.code == code and not found_job:
            this_job = jb
            found_job = True
    if found_job: return found_job, this_job
    else: return found_job, None


def register_job(this_job, lines: str, software='G16', print_output: bool=True):
    isfinished = False
    isgood = False
    if os.path.isfile(this_job.output_path): isfinished = True
    if isfinished:

        ###################
        ### Gaussian 16 ###
        ###################
        if software == 'G16':
            line_num, found_good = search_string("Normal termination", lines, typ='all')
            line_num, found_time = search_string("Elapsed time:", lines, typ='last')
            if found_time:
                elapsed_time_list = lines[line_num].split()[2:]
                elapsed_time = G16_time_to_sec(elapsed_time_list)
            else: elapsed_time = float(0)
        else:
            found_good = False
            print(f"Registry of this software is not implemented. See register_jobs in Scope Classes")
        ###################

        if found_good: isgood = True
        else: isgood = False

    this_job.isregistered = True
    this_job.isfinished = isfinished
    this_job.isgood = isgood
    this_job.elapsed_time = elapsed_time

    if print_output:
        print(f"    REGISTERED for Job: {this_job.code}. LOG:{this_job.output_path}")
        if this_job.issubmitted: print(f"     -> Submitted")
        if this_job.isfinished:  print(f"     -> Finished")
        if this_job.isgood:      print(f"     -> Good")
        else:               print(f"     -> BAD !! ")
        if this_job.isfinished:  print(f"     -> Elapsed Time: {this_job.elapsed_time}")
