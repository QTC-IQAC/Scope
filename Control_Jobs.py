#!/usr/bin/env python3
import sys
import os
import pwd
import subprocess
import numpy as np

from Scope.Parse_General import read_lines_file, search_string 

def set_user():
    return pwd.getpwuid( os.getuid() ).pw_name

def set_cluster():
    return os.uname()[1]

def send_command(commandtype: str, filename: str=None, cluster: str=set_cluster(), user: str=set_user(), queue: str='', debug: int=0):
    if 'portal' in cluster:
        if commandtype == "qstat": raw = subprocess.check_output(['bash','-c', "qstat"])
        elif commandtype == "queue_stat": raw = subprocess.check_output(['bash','-c', 'qstat -f | grep'+queue])
        elif commandtype == "check_job":  raw = subprocess.check_output(['bash','-c', "qstat -xml"]) #-q "+queue+".q" ]) 
        elif commandtype == "submit":     subprocess.run(['bash','-c', 'qsub '+filename]) 
    elif 'login' in cluster or 'csuc' in cluster:
        if commandtype == "qstat":
            tmp = 'squeue -o "%.9P %.50j %.12u %.2t %.12M %.5C %.3D %R" | grep '+str(user)
            try:    raw = subprocess.check_output(['bash','-c', tmp])
            except: raw = ""
        elif commandtype == "queue_stat":
            try:    raw = subprocess.check_output(['bash','-c', 'sinfo | grep '+queue+' grep idle'])
            except: raw = ""
        elif commandtype == "check_job":
            tmp = 'squeue -o "%.60j %.12u"'
            try:    raw = subprocess.check_output(['bash','-c', tmp ]) 
            except: raw = ""
        elif commandtype == "submit":     subprocess.run(['bash','-c', 'sbatch '+filename]) 
    else: print("Error in send_command function. Cluster not recognizsed")
    if commandtype == "qstat" or commandtype == "queue_stat" or commandtype == "check_job": return raw


def set_queues(cluster: str=set_cluster(), queues: str='all'):
    list_q = []
    if 'login' in cluster or 'csuc' in cluster:
        list_q.append("std")
    elif 'portal' in cluster:
        if queues == 'all':
                for i in range(1,11): 
                    if i == 3 or i == 1: pass
                    elif i < 10: list_q.append(str("iqtc0"+str(i)))
                    else: list_q.append(str("iqtc"+str(i)))
        else:
            digested = queues.rstrip(" ").lstrip(" ").split(',')
            for i in range(1,11): 
                if i == 3 or i == 1: pass
                elif i != 3 and str(i) in digested: 
                    if i < 10: list_q.append(str("iqtc0"+str(i)))
                    else: list_q.append(str("iqtc"+str(i)))
    return list_q     

def check_usage(cluster: str=set_cluster(), user: str=set_user(), queues: str='all', debug: int=0):
    cpus = 0
    jobs = 0
    #raw = subprocess.check_output(['bash','-c', "qstat"])
    raw = send_command("qstat", cluster=cluster,user=user,debug=debug)

    if 'login' in cluster or 'csuc' in cluster:
        try: 
            #print("Raw:", raw)
            #text = raw.rstrip().split("\n")
            dec = raw.decode("utf-8")
            text = dec.rstrip().split("\n")
            #print("text:", text)
            for line in text:
                blocks = line.split()
                cpus += int(blocks[5])
                jobs += 1 
        except: 
            cpus = int(0)
            jobs = int(0)

    elif 'portal' in cluster:
        dec = raw.decode("utf-8")
        text = dec.rstrip().split("\n")
        list_q = set_queues(cluster, queues)
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

def check_queue_availability(cluster: str=set_cluster(), queues: str='all', debug: int=0):
    list_q = set_queues(cluster, queues)
    list_q_worked = []
    list_empty_cpus = []
    list_total_empty = []
    for q in list_q:
        try:
            raw = send_command("queue_stat",cluster,queue=q, debug=debug) 
            dec = raw.decode("utf-8") 
            empty_cpus = []

            #if 'login' in cluster or 'csuc' in cluster:
            #    text = dec.rstrip().split("\n")
            #    if len(text) > 0: list_total_empty = int(1) 

            if 'portal' in cluster:
                text = dec.rstrip().split("\n")
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

##########################################
### Keeping the old version, pre CSUC: ###
##########################################

#def check_submitted_job_xml(name: str, queues: str='all', debug: int=0):
#    issubmitted = False
#    list_q = set_queues(queues)
#    if debug >= 1: print("QUEUES:", list_q)
#    for q in list_q:
#        raw = subprocess.check_output(['bash','-c', "qstat -xml -q "+q+".q" ])
#        dec = raw.decode("utf-8")
#        flat = dec.replace("\n", "")
#        code = str("<JB_name>"+name+"</JB_name>")
#        if debug >= 1: print("code:", code)
#        if debug >= 1: print("flat:", flat)
#        if code in flat:
#            issubmitted = True
#            return issubmitted
#        else: issubmitted = False
#    ## Checks Pending Jobs in all queues
#    if not issubmitted:
#       raw = subprocess.check_output(['bash','-c', "qstat -xml" ])
#       dec = raw.decode("utf-8")
#       flat = dec.replace("\n", "")
#       code = str("<JB_name>"+name+"</JB_name>")
#       if code in flat:
#           issubmitted = True
#           return issubmitted
#    return issubmitted

def check_submitted_job(name: str, cluster: str=set_cluster(), debug: int=0):
    issubmitted = False

    raw = send_command("check_job", cluster, debug=debug)
    dec = raw.decode("utf-8") 
    flat = dec.replace("\n", "")

    if 'login' in cluster or 'csuc' in cluster:
        if name in flat: issubmitted = True 
        else:            issubmitted = False

    elif 'portal' in cluster:
        code = str("<JB_name>"+name+"</JB_name>") 
        if code in flat: issubmitted = True 
        else:            issubmitted = False

    return issubmitted

def get_queue_and_procs(resources: str="light", cluster: str=set_cluster(), debug: int=0):
        #### Finds best number of processors and queue
    if "portal" in cluster:
        if resources.lower() == "light": mult = 1
        elif resources.lower() == "medium": mult = 2
        elif resources.lower() == "heavy": mult = 4
        else: mult = 1

        askqueue = set_best_queue('8,9,10')
        if askqueue == 'iqtc08': askprocs = 7*mult
        else: askprocs = 8*mult

    elif "login" in cluster or "csuc" in cluster:

        askqueue = "std"
        if resources.lower() == "light": mult = 1
        elif resources.lower() == "medium": mult = 2
        elif resources.lower() == "heavy": mult = 4

        askprocs = 8*mult
    else:
        askqueue = ''
        askprocs = 0

    return askqueue, askprocs
#####################################
##### Portal Specific Functions #####
#####################################
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

def check_fairsharing(user: str=set_user()):
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

