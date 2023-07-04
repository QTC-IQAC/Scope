#!/usr/bin/env python3
import sys
import os
import pwd
import subprocess
import numpy as np

from Scope.Parse_General import read_lines_file, search_string 
#from Scope.Classes_Job import check_job_requisites, find_job

def set_user():
    return pwd.getpwuid( os.getuid() ).pw_name

def set_cluster():
    return os.uname()[1]

########################
def set_paths(cluster: str=set_cluster(), user: str=set_user(), debug: int=0):
    if 'login' in cluster or 'csuc' in cluster:
        cell2mol_path = "/scratch/svela/SCOPE/Database_SCO/4-Merged"
        calcs_path    = "/home/svela/SCOPE/Database_SCO/5-Complexes_Iso"
    elif 'portal' in cluster or 'node' in cluster or 'visual' in cluster:
        if 'g2vela' in user:
            cell2mol_path = "/home/g4vela/SCOPE/Database_SCO/4-Merged"
            calcs_path    = "/home/g2vela/SCOPE/Database_SCO/5-Complexes_Iso"
        elif 'g4vela' in user:
            cell2mol_path = "/home/g4vela/SCOPE/Database_SCO/4-Merged"
            calcs_path    = "/home/g4vela/SCOPE/Database_SCO/5-Complexes_Iso"
    elif 'lemma' in cluster:
        cell2mol_path = "/Users/sergivela/Documents/SCOPE/Database_SCO/4-Merged"
        calcs_path    = "/Users/sergivela/Documents/SCOPE/Database_SCO/5-Complexes_Iso"
    elif 'uam' in cluster:
        cell2mol_path = "/home/proyectos/ub100/SCOPE/4-Merged"
        calcs_path    = "/home/proyectos/ub100/SCOPE/5-Complexes_Iso"
    else: 
        print(f"Cluster {cluster} not recognized")
    if cell2mol_path[-1] != '/': cell2mol_path += '/'
    if calcs_path[-1] != '/': calcs_path += '/'
    return cell2mol_path, calcs_path

########################
def set_PP_Library(cluster: str=set_cluster(), user: str=set_user(), debug: int=0):
    if 'login' in cluster or 'csuc' in cluster:
        PP_Library= "/home/svela/Programes/PP_Library"
    elif 'portal' in cluster:
        PP_Library= "/home/g4vela/Programes/PP_Library"
    elif 'uam' in cluster:
        PP_Library= "/lustre/home/ub100435/Programes/PP_Library"
    else: 
        print("Cluster not recognized")
    return PP_Library

########################
def send_command(commandtype: str, filename: str=None, cluster: str=set_cluster(), user: str=set_user(), queue: str='', debug: int=0):
    if debug > 0: print("SEND_COMMAND: evaluating", commandtype, "for", cluster)
    if 'portal' in cluster:
        if commandtype == "qstat": raw = subprocess.check_output(['bash','-c', "qstat"])
        elif commandtype == "queue_stat": 
            if user[2].isdigit(): group = user[0:3]
            else:               group = user[0:2]
            tmp = str(f"qstat -f | grep {queue} | grep {group}")
            if debug > 0: print("SEND_COMMAND: command=", tmp)
            raw = subprocess.check_output(['bash','-c', tmp ])
        elif commandtype == "check_job":  raw = subprocess.check_output(['bash','-c', "qstat -xml"]) #-q "+queue+".q" ]) 
        elif commandtype == "submit":     subprocess.run(['bash','-c', 'qsub '+filename]) 
    elif 'login' in cluster or 'csuc' in cluster or 'uam' in cluster:
        if commandtype == "qstat":
            tmp = 'squeue -o "%.9P %.50j %.12u %.2t %.12M %.5C %.3D %R" | grep '+str(user)
            try:  
                raw = subprocess.check_output(['bash','-c', tmp ])
            except: 
                raw = subprocess.check_output(['bash','-c', 'echo']) 
        elif commandtype == "queue_stat":
            tmp = 'sinfo | grep '+queue+' | grep idle' 
            raw = subprocess.check_output(['bash','-c', tmp ])
        elif commandtype == "check_job":
            tmp = 'squeue -o "%.60j %.12u"'
            raw = subprocess.check_output(['bash','-c', tmp ]) 
        elif commandtype == "submit":     subprocess.run(['bash','-c', 'sbatch '+filename]) 
    elif 'lemma' in cluster:
        if commandtype == "qstat": raw = '' 
        elif commandtype == "queue_stat": raw = ''
        elif commandtype == "check_job": raw = '' 
        elif commandtype == "submit": pass
    else: print("Error in send_command function. Cluster not recognised"); return ''
    if commandtype == "qstat" or commandtype == "queue_stat" or commandtype == "check_job": 
        return raw

########################
def set_queues(cluster: str=set_cluster(), queues: str='all', debug: int=0):
    list_q = []
    list_of_exceptions = [1, 3, 5, 7]
    if 'login' in cluster or 'csuc' in cluster:
        list_q.append("std")
    elif 'uam' in cluster:
        list_q.append("class_a")
    elif 'portal' in cluster:
        if queues == 'all':
                for i in range(1,11): 
                    if i in list_of_exceptions: pass
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

########################
def check_usage(cluster: str=set_cluster(), user: str=set_user(), queues: str='all', debug: int=0):
    cpus = 0
    jobs = 0
    #raw = subprocess.check_output(['bash','-c', "qstat"])
    raw = send_command("qstat", cluster=cluster,user=user,debug=debug)

    if 'login' in cluster or 'csuc' in cluster or 'uam' in cluster:
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
        list_q = set_queues(cluster=cluster, queues=queues, debug=debug)
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

#####################
## PORTAL SPECIFIC ##
#####################
def check_queue_availability(queues: str='all', cluster: str=set_cluster(), debug: int=0):
    list_q = set_queues(cluster=cluster, queues=queues, debug=debug)
    list_q_worked = []
    list_empty_cpus = []
    list_total_empty = []
    list_ratio = []
    list_score = []
    for q in list_q:
        if   q == 'iqtc02': value = int(1)
        elif q == 'iqtc04': value = int(2)
        elif q == 'iqtc06': value = int(4)
        elif q == 'iqtc08': value = int(8)
        elif q == 'iqtc09': value = int(8)
        elif q == 'iqtc10': value = int(8)
        else:               value = int(0)
        try:
            if debug > 0: print("CHECK_QUEUE_AVAIL: sending queue_stat for q=",q)
            raw = send_command("queue_stat",cluster,queue=q, debug=debug) 
            dec = raw.decode("utf-8") 
            empty_cpus = []
            total_nodes = int(0)

            if 'portal' in cluster:
                text = dec.rstrip().split("\n")
                for node in text:
                    blocks = node.split()
                    queue = blocks[0].split('@')[0]
                    queue_check = blocks[0].split('@')[0][0:6]
                    node = blocks[0].split('@')[1]
                    if debug > 1: print("CHECK_QUEUE_AVAIL: queue=", queue)
                    if debug > 1: print("CHECK_QUEUE_AVAIL: queue_check=", queue_check)
                    if debug > 1: print("CHECK_QUEUE_AVAIL: node=", node)
                    if queue_check == q:
                        reserved, used, total = blocks[2].split('/')
                        if debug > 1: print("    reserved=", reserved)
                        if debug > 1: print("    used=", used)
                        if debug > 1: print("    total=", total)
                        total_nodes += int(total)
                        if queue == q+".q": empty_cpus.append(int(total)-int(reserved)-int(used)) 
                    else: print("queue_check not passed for", queue_check, q)
                list_q_worked.append(q)
                list_empty_cpus.append(empty_cpus)
                list_total_empty.append(np.sum(empty_cpus))
                list_ratio.append(float(np.sum(empty_cpus))/float(total_nodes))
                list_score.append(float(np.sum(empty_cpus))*float(value)/float(total_nodes))

        except: print(f"CHECK_QUEUE_AVAIL: error evaluating queue {q}")
    return list_q_worked, list_empty_cpus, list_total_empty, list_ratio, list_score

########################
def check_submitted(name: str, cluster: str=set_cluster(), debug: int=0):
    issubmitted = False

    if debug > 0: print("check_submitted: cluster is:", cluster)
    if 'node' in cluster or 'lemma' in cluster: 
        issubmitted = False
    else: 
        raw = send_command("check_job", cluster, debug=debug)
        dec = raw.decode("utf-8") 
        flat = dec.replace("\n", "")
    
        if 'login' in cluster or 'csuc' in cluster: 
            if name in flat: issubmitted = True 
            else:            issubmitted = False
        elif 'uam' in cluster:
            if name in flat: issubmitted = True 
            else:            issubmitted = False
        elif 'portal' in cluster:
            code = str("<JB_name>"+name+"</JB_name>") 
            if code in flat: issubmitted = True 
            else:            issubmitted = False
        else: print("    Check_Submitted_Job: Cluster not recognized")
    
    return issubmitted

########################
def get_queue_and_procs(environment: object, debug: int=0):
    #### Finds best number of processors and queue

    # Rename to simplify call
    resources = environment.resources
    min_procs = environment.min_procs
    method    = environment.method
    cluster   = environment.cluster
    queues    = environment.queues

    if "portal" in cluster:
        ## Multiplier
        if   resources == "light":  mult = 1
        elif resources == "medium": mult = 2
        elif resources == "heavy":  mult = 4
        else: mult = 1

        askqueue, maxprocs = set_best_queue(min_procs=min_procs, method=method, queues=queues, debug=debug)
        if askqueue == 'iqtc08': askprocs = 7*mult
        elif askqueue == 'iqtc04': askprocs = 6*mult
        elif askqueue == 'iqtc02': askprocs = 4*mult
        else: askprocs = 8*mult

    elif "login" in cluster or "csuc" in cluster: 
        askqueue = "std"
        if   resources == "light":  mult = 1
        elif resources == "medium": mult = 2
        elif resources == "heavy":  mult = 4
        askprocs = 8*mult

    elif 'uam' in cluster:
        askqueue = "class_a"
        if   resources == "light":  mult = 1
        elif resources == "medium": mult = 2
        elif resources == "heavy":  mult = 4
        askprocs = 8*mult

    else:
        askqueue = ''
        askprocs = 0

    return askqueue, askprocs

#####################################
##### Portal Specific Functions #####
#####################################
def set_best_queue(min_procs: int=6, method: str='score', queues: str='all', debug: int=0):
    list_q, list_empty_cpus, list_total_empty, list_ratio, list_score = check_queue_availability(queues=queues, debug=debug)

    if debug > 0: print("SET_BEST_QUEUE: results of check_queue_availability:") 
    if debug > 0: print(list_q)
    if debug > 0: print(list_total_empty)
    if debug > 0: print(list_ratio)
    if debug > 0: print(list_score)

    queue_with_most_empty = np.argmax(list_total_empty)
    queue_with_best_ratio = np.argmax(list_ratio)
    queue_with_best_score = np.argmax(list_score)
    if debug > 0: print("SET_BEST_QUEUE: queue_with_most_empty:", list_q[queue_with_most_empty], np.max(list_total_empty))
    if debug > 0: print("SET_BEST_QUEUE: queue_with_best_ratio:", list_q[queue_with_best_ratio], np.max(list_ratio))
    if debug > 0: print("SET_BEST_QUEUE: queue_with_best_score:", list_q[queue_with_best_score], np.max(list_score))

    if method == 'total':
        found = False
        while not found:
            max_empty = np.max(list_empty_cpus[queue_with_most_empty]) 
            if debug > 0: print("SET_BEST_QUEUE: max_empty:", max_empty) 
            for idx, cp in enumerate(list_empty_cpus[queue_with_most_empty]):
                if cp >= min_procs: 
                    found = True
                    return str(list_q[queue_with_most_empty]), max_empty
            list_total_empty[queue_with_most_empty] = 0
            queue_with_most_empty = np.argmax(list_total_empty)
            if np.sum(list_total_empty) == 0: break

    elif method == 'ratio':
        found = False
        while not found:
            max_ratio = np.max(list_ratio[queue_with_best_ratio]) 
            num_empty = np.max(list_empty_cpus[queue_with_best_ratio]) 
            if debug > 0: print("SET_BEST_QUEUE: max_ratio:", max_ratio) 
            for idx, cp in enumerate(list_empty_cpus[queue_with_best_ratio]):
                if cp >= min_procs: 
                    found = True
                    return str(list_q[queue_with_best_ratio]), num_empty
            list_ratio[queue_with_best_ratio] = 0.0
            queue_with_best_ratio = np.argmax(list_ratio)
            if np.sum(list_ratio) == 0: break

    elif method == 'score':
        found = False
        while not found:
            max_score = np.max(list_score[queue_with_best_score])
            num_empty = np.max(list_empty_cpus[queue_with_best_score])
            if debug > 0: print("SET_BEST_QUEUE: max_score:", max_score)
            for idx, cp in enumerate(list_empty_cpus[queue_with_best_score]):
                if cp >= min_procs:
                    found = True
                    return str(list_q[queue_with_best_score]), num_empty
            list_ratio[queue_with_best_score] = 0.0
            queue_with_best_score = np.argmax(list_score)
            if np.sum(list_score) == 0: break

    else: print("SET BEST QUEUE: wrong method specified"); return None, 0   

