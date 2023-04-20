#!/usr/bin/env python3
import sys
import os
import pwd
import subprocess
import numpy as np

from Test_V3.Workflow import Recipe, Job, Computation 
from Test_V3.Environment import check_usage, get_queue_and_procs, send_command
from Test_V3.Write_G16_Inputs import *
from Test_V3.Write_QE_Inputs import *
from Test_V3.Write_ORCA_Inputs import *

######################
##### New for V3 #####
######################
def run_computation(sys: object, gmol: object, resources: object, options: object, job_data: object, qc_data: object, debug: int=0):

    if options.want_submit: sent_procs, sent_jobs = check_usage()
    else:                   sent_procs = 0; sent_jobs = 0

    if sent_procs >= resources.max_procs or sent_jobs >= resources.max_jobs: 
        if debug > 0: print(f"    Over maximum jobs and cores reached")
        return gmol, False

    updated = False
    ## 1-Finds the recipe and Evaluates the State of Queues. If recipe does not exist, it creates it.
    exists, recipe = sys.get_recipe(gmol, job_data.recipe_code, debug=debug)
    if not exists: updated = True

    ## 2-Checks that all requisites and constrains of the job are fulfilled
    cancontinue = check_job_requisites(recipe, job_data)
    if not cancontinue and debug > 0: print(f"    Requisites not met, or job already run!"); return gmol, False

    ## 3-Determines Run Number and sets it
    run_number = get_run_number(recipe, job_data)
    setattr(job_data,"run_number",run_number)

    ## 4-Sets_Paths of Files
    paths = set_path_data(sys, gmol, recipe, job_data)

    ## 5-Sets Job.  It does NOT create new jobs, unlike get_recipe, so it must be done separately
    exists, this_job = find_job(recipe,job_data)
    if not exists: this_job = job(job_data, paths); recipe.jobs.append(this_job); updated = True

    ## 6-Evaluates Submission
    if options.want_submit:
        ## 6.1-Defines resources
        askqueue, askprocs = get_queue_and_procs(resources=job_data.resources)
        ## 6.2-Adds resources to job
        this_job.add_submission_init(nprocs=askprocs, queue=askqueue)
        ## 6.3-Creates Files
        if not this_job.input_exists or options.overwrite_inputs:
            if this_job.software == 'g16':   gen_G16_input(gmol, qc_data, paths, debug=debug)
            elif this_job.software == 'qe':  pass # gen_G16_input(gmol, qc_data, paths, debug=debug)
        if not this_job.subfile_exists or options.overwrite_inputs:
            if this_job.software == 'g16':   gen_G16_subfile(paths, procs=askprocs, queue=askqueue)
            elif this_job.software == 'qe':  pass # gen_G16_subfile(paths, procs=askprocs, queue=askqueue)

    ## 7-Evaluates if the output exists
    this_job.add_submission_end(want_to_read=False)

    ## 8-If it exists, prompts for registration
    if this_job.output_exists and not this_job.isregistered:
        print(f"    Output file Found Pending to be REGISTERED")
        print(f"    {this_job.output_path}")
        print(f"    ")

    ## 9-Evaluates Submission
    if options.want_submit:
        can_submit = True
        ## 9.1-Evaluates if output exists
        if this_job.output_exists and not options.overwrite_logs:
            can_submit = False
            if debug > 0: print("Output exists and not overwriting logs")

        ## 9.2-Evaluates if output is running
        if can_submit and not options.ignore_submitted:
            this_job.check_submission_status()   ### retrieves this_job.isrunning
            if this_job.isrunning:
                can_submit = False
                if debug > 0: print("Job already running")

        ## 9.3-Submits if possible
        if can_submit:
            os.chdir(recipe.path)
            send_command("submit", filename=paths.sub_name)
            if debug > 0: print(f"Job {this_job.output_path} submitted")

    return gmol, updated
