#!/usr/bin/env python3
import sys
import os
import subprocess
import numpy as np

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
