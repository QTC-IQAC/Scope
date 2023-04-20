#!/usr/bin/env python3
import sys
import os
import pwd

from Test_V3.Classes_Input import *
from Test_V3.Classes_SCO import sco_system, crystal
from Test_V3.Workflow import Recipe, Job, Computation 
from Test_V3.Workflow.Recipe import *
from Test_V3.Workflow.Job import *
from Test_V3.Workflow.Computation import *
from Test_V3.Environment import check_usage, get_queue_and_procs, send_command, set_cluster, set_user
from Test_V3.Read_Write import load_binary, save_binary 

######################
def execute_job(sys_path: str, job_path: str, debug: int=0):

    #### 0-Verifies Files
    files = False
    if os.path.isfile(sys_path) and os.path.isfile(job_path): files = True 
    if debug > 1: print("Execute_JOB, step 0: sys_path=", sys_path)
    if debug > 1: print("Execute_JOB, step 0: job_path=", job_path)
    if debug > 1: print("Execute_JOB, step 0: files=", files)
    if not files: return None

    #### 1-Reads Input Data
    resources = set_resources_data(job_path, section="&resources", debug=0)
    options   = set_options_data(job_path, section="&options"    , debug=0)
    job_data  = set_job_data(job_path, section="&job_data"       , debug=0)
    qc_data   = set_qc_data(job_path, section="&qc_data"         , debug=0)

    #### 2-Fixes Some Data Depending on Cluster
    cluster = set_cluster()
    if 'lemma' in cluster or 'node' in cluster:
        options.dct['want_submit'] = False         ## I think I tried options.want_submit = False and didn't work
        options.dct['overwrite_inputs'] = False    ## " " same above
        options.dct['overwrite_logs'] = False      ## " " same above
    
    ## 3.1-Loads the System
    if debug > 1: print(" ")
    sys = load_binary(sys_path)
    updated = False
    if debug > 1: print("Execute_JOB, step 3: system loaded")
    
    ## 3.2-Sets instructions for both HS and LS
    if job_data.recipe_code.lower()   == "solid":    object_list = list([sys.HS_ref_crys.cell,sys.LS_ref_crys.cell])
    elif job_data.recipe_code.lower() == "isolated": object_list = list([sys.HSref,sys.LSref])
    for idx, gmol in enumerate(object_list):

        ############## 
        ### RECIPE ###
        ##############
        ## 4-Finds the recipe. If it does not exist, it is created.
        exists, recipe = sys.get_recipe(gmol, job_data.recipe_code, debug=debug)
        if exists and debug > 1: print("Execute_JOB, step 4: recipe found with keyword:", recipe.keyword)
        if not exists: 
            updated = True
            if debug > 1: print("Execute_JOB, step 4: recipe created with keyword:", recipe.keyword)
            if debug > 1: print("Execute_JOB, step 4: recipe folder created at:", recipe.path)

        ###########
        ### JOB ###
        ###########
        ## 5-Finds the job. If it does not exist, it is NOT created.
        if debug > 0: recipe.get_info()
        exists, this_job = recipe.find_job(job_data, debug=debug)
    
        ## 5.1 If necessary, creates the job
        if not exists: 
            this_job = recipe.add_job(job_data); updated = True
            if debug > 1: print("Execute_JOB, step 5.1: job does not exist, creating")
           
        ## 6-Checks that all requisites and constrains of the job are fulfilled, and sets run_number
        cancontinue = this_job.check_requisites(debug=debug)
        if not cancontinue and debug > 1: print(f"    Requisites not met, or job already run!")
        if not cancontinue: continue    # I know if might seem misleading. Here, continue means "skip this one"
        this_job.set_run_number()
        if debug > 1: print("Execute_JOB, step 6: run_number set:", this_job.run_number)

        ####################
        ### COMPUTATIONS ###
        ####################
        ## 7-Sets the computation(s), meaning that it will check if they exist, and if not, it creates them
        this_job.set_computations_from_type(qc_data)
        if debug > 1: print("Execute_JOB, step 7: computations set")
        if debug > 0: this_job.get_info()
    
        for jdx, comp in enumerate(this_job.computations):

            ## 8.1-Checks files 
            if not hasattr(comp,"output_exists"): comp.check_files()
            ## 8.2-Evaluates Submission
            if not comp.output_exists:
                comp.check_submission_status()
                if not comp.isrunning:             comp.run(resources, options, debug=debug); updated = True
            else:
                ## 8.3-If output exists, and is not registered, it does it
                if not comp.isregistered:          comp.register(debug=debug);     updated = True
                #if not comp._job.isregistered:     this_job.register(debug=debug); updated = True
            if debug > 0: comp.get_info() 

        # Updates Job Registry information
        this_job.register(debug=debug)

    if updated: save_binary(sys, sys_path)
    return None


#def run_computation(sys: object, gmol: object, resources: object, options: object, job_data: object, qc_data: object, debug: int=0):
#
#    if options.want_submit: sent_procs, sent_jobs = check_usage()
#    else:                   sent_procs = 0; sent_jobs = 0
#
#    if sent_procs >= resources.max_procs or sent_jobs >= resources.max_jobs: 
#        if debug > 0: print(f"    Over maximum jobs and cores reached")
#        return gmol, False
#
#    updated = False
#    ## 1-Finds the recipe and Evaluates the State of Queues. If recipe does not exist, it creates it.
#    exists, recipe = sys.get_recipe(gmol, job_data.recipe_code, debug=debug)
#    if not exists: updated = True
#
#    ## 2-Checks that all requisites and constrains of the job are fulfilled
#    cancontinue = check_job_requisites(recipe, job_data)
#    if not cancontinue and debug > 0: print(f"    Requisites not met, or job already run!"); return gmol, False
#
#    ## 3-Determines Run Number and sets it
#    run_number = get_run_number(recipe, job_data)
#    setattr(job_data,"run_number",run_number)
#
#    ## 4-Sets_Paths of Files
#    paths = set_path_data(sys, gmol, recipe, job_data)
#
#    ## 5-Sets Job.  It does NOT create new jobs, unlike get_recipe, so it must be done separately
#    exists, this_job = find_job(recipe,job_data)
#    if not exists: this_job = job(job_data, paths); recipe.jobs.append(this_job); updated = True
#
#    ## 6-Evaluates Submission
#    if options.want_submit:
#        ## 6.1-Defines resources
#        askqueue, askprocs = get_queue_and_procs(resources=job_data.resources)
#        ## 6.2-Adds resources to job
#        this_job.add_submission_init(nprocs=askprocs, queue=askqueue)
#        ## 6.3-Creates Files
#        if not this_job.input_exists or options.overwrite_inputs:
#            if this_job.software == 'g16':   gen_G16_input(gmol, qc_data, paths, debug=debug)
#            elif this_job.software == 'qe':  pass # gen_G16_input(gmol, qc_data, paths, debug=debug)
#        if not this_job.subfile_exists or options.overwrite_inputs:
#            if this_job.software == 'g16':   gen_G16_subfile(paths, procs=askprocs, queue=askqueue)
#            elif this_job.software == 'qe':  pass # gen_G16_subfile(paths, procs=askprocs, queue=askqueue)
#
#    ## 7-Evaluates if the output exists
#    this_job.add_submission_end(want_to_read=False)
#
#    ## 8-If it exists, prompts for registration
#    if this_job.output_exists and not this_job.isregistered:
#        print(f"    Output file Found Pending to be REGISTERED")
#        print(f"    {this_job.output_path}")
#        print(f"    ")
#
#    ## 9-Evaluates Submission
#    if options.want_submit:
#        can_submit = True
#        ## 9.1-Evaluates if output exists
#        if this_job.output_exists and not options.overwrite_logs:
#            can_submit = False
#            if debug > 0: print("Output exists and not overwriting logs")
#
#        ## 9.2-Evaluates if output is running
#        if can_submit and not options.ignore_submitted:
#            this_job.check_submission_status()   ### retrieves this_job.isrunning
#            if this_job.isrunning:
#                can_submit = False
#                if debug > 0: print("Job already running")
#
#        ## 9.3-Submits if possible
#        if can_submit:
#            os.chdir(recipe.path)
#            send_command("submit", filename=paths.sub_name)
#            if debug > 0: print(f"Job {this_job.output_path} submitted")
#
#    return gmol, updated
