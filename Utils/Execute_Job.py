#!/usr/bin/env python3
import sys
import os
import pwd

from Scope.Classes_Input import *
from Scope.Classes_SCO import sco_system, crystal
from Scope.Workflow import Recipe, Job, Computation
from Scope.Workflow.Recipe import *
from Scope.Workflow.Job import *
from Scope.Workflow.Computation import *
#from Scope.Environment import check_usage, get_queue_and_procs, send_command, set_cluster, set_user
from Scope.Read_Write import load_binary, save_binary

######################
def execute_job(sys_path: str, job_path: str, handle_errors: bool=False, debug: int=0):

    print("ENTERED EXECUTE JOB with job_path:", job_path)
    report = ''

    if debug > 1: print("")
    if debug > 1: print("----------- NEW JOB ----------")
    if debug > 1: print("")

    #### 0-Verifies Files
    files = False
    if os.path.isfile(sys_path) and os.path.isfile(job_path): files = True
    if debug > 1: print("Execute_JOB, step 0: sys_path=", sys_path)
    if debug > 1: print("Execute_JOB, step 0: job_path=", job_path)
    if debug > 1: print("Execute_JOB, step 0: files=", files)
    if debug > 1: print("------------------------------")
    if not files: return None

    #### 1-Reads Input Data
    environment = set_environment_data(job_path, section="&environment", debug=0)
    options     = set_options_data(job_path, section="&options"      , debug=0)
    job_data    = set_job_data(job_path, section="&job_data"         , debug=0)
    qc_data     = set_qc_data(job_path, section="&qc_data"           , debug=0)

    #### 2-Fixes Some Data Depending on Cluster
    cluster = set_cluster()
    if 'lemma' in cluster or 'node' in cluster:
        options.dct['want_submit'] = False         ## I think I tried options.want_submit = False and didn't work
        options.dct['overwrite_inputs'] = False    ## " " same above
        options.dct['overwrite_logs']   = False      ## " " same above

        options.want_submit      = False         ## I think I tried options.want_submit = False and didn't work
        options.overwrite_inputs = False    ## " " same above
        options.overwrite_logs   = False      ## " " same above

    ##############
    ### SYSTEM ###
    ##############
    ## 3.1-Loads the System
    if debug > 1: print(" ")
    sys = load_binary(sys_path)
    updated = False
    if debug > 1: print("Execute_JOB, step 3: system loaded")

    ##############
    ### BRANCH ###
    ##############
    exists, this_branch = sys.find_branch(job_data.branch, debug=debug)
    if debug > 1 and not exists: print("Execute_JOB, step 4: creating branch")
    if not exists: this_branch = sys.add_branch(job_data.branch, job_data.target, debug=debug); updated = True
    if debug > 1: print("Execute_JOB, step 4: branch loaded")

    ##############
    ### RECIPE ###
    ##############
    for idx, recipe in enumerate(this_branch.recipes):

        ###########
        ### JOB ###
        ###########
        ## 5-Finds the job. If it does not exist, it is created in 5.1
        exists, this_job = recipe.find_job(job_data=job_data, debug=debug)

        ## 5.1 If job does not exist, creates the job
        if not exists: this_job = recipe.add_job(job_data); updated = True
        ## 5.2 If job exists, it checks for input changes (also in qc_data for the computations)
        else:
            this_job.check_input(job_path=job_path, debug=debug)
        if debug > 1: print("---------------------------------------------------")
        if debug > 1: print("Execute_JOB, step 5: job loaded for spin", this_job._recipe.subject.spin)
        if debug > 1: print("---------------------------------------------------")

        ## 6-Checks that all requisites and constrains of the job are fulfilled
        cancontinue = this_job.check_requisites(debug=debug)
        if not cancontinue:
            if debug > 1:   print("  Execute_JOB, step 6: requisites NOT met or job already run!")
            continue        # I know if might seem misleading. Here, "continue" means "skip this one"
        else:
            if debug > 1:   print("  Execute_JOB, step 6: requisites fulfilled")

        ####################
        ### COMPUTATIONS ###
        ####################
        ## 7-Sets the computation(s), meaning that it will check if they exist, and if not, it creates them
        this_job.set_computations_from_setup(qc_data, debug=debug)
        if debug > 1: print("Execute_JOB, step 7: computations set:")

        for jdx, comp in enumerate(this_job.computations):

            #if comp.has_update and comp.isregistered: continue # Skip jobs with update (i.e. with other related computations with higher run_number)
            if debug > 1: print("-----------------------------------------------------------------------------------")
            if debug > 1: print(f"Execute_JOB, step 7.0: evaluating job, and computation with indices: {recipe.jobs.index(this_job)+1}/{len(recipe.jobs)}, {jdx+1}/{len(this_job.computations)}")
            ## 8.1-Checks files
            if debug > 1: print("Execute_JOB, step 7.1: doing computation with keyword and run_number:", comp.keyword, comp.run_number)
            if debug > 1: print("Execute_JOB, step 7.1: out_file:", comp.out_path)
            if debug > 1: print("Execute_JOB, step 7.1: is_update:", comp.is_update)
            if debug > 1: print("-----------------------------------------------------------------------------------")

            comp.check_files()

            ## 8.2-Evaluates Submission
            if debug > 1: print("Execute_JOB, step 7.2: checking files [inp, out, sub]:", comp.input_exists, comp.output_exists, comp.subfile_exists)
            if not comp.output_exists: # and comp.input_exists:
                if options.want_submit:
                    comp.check_submission_status(debug=debug)
                    if debug > 1: print("Execute_JOB, step 7.3a: checking submission status: isrunning=",comp.isrunning)
                    if debug > 1: print("Execute_JOB, step 7.3a: initial state is", comp.qc_data.istate)
                    if debug > 1: print("Execute_JOB, step 7.3a: is_update:", comp.is_update)
                    if not comp.isrunning:             comp.run(environment, options, debug=debug); updated = True
                else: 
                    if debug > 1: print("Execute_JOB, step 7.3a: want_submit is False")
            elif comp.output_exists and not comp.input_exists:  
                report += f"Investigate {comp.out_path} \n"
                print(f"Investigate {comp.out_path}")
            else:
                ## 8.3-If output exists, and is not registered, it does it
                if not comp.isregistered:
                    worked = comp.register(debug=debug)

                    # If registration fails, either...
                    if not worked: 
                        if handle_errors: # ...takes default action 
                            comp.read_lines()
                            if len(comp.output_lines) > 0: comp.store(debug=debug)  # Creates Copy of output 
                            else:                          this_job.remove_computation(comp_index=comp.index)
                            report += f"Errors handled for {comp.out_path} \n"
                            print(f"Errors handled for {comp.out_path}")
                        else:              # ...or warns the user
                            report += f"Check Registration of {comp.out_path} \n"
                            print(f"Check Registration of {comp.out_path}")
                    updated = True

                ## 8.4-If output exists, is registered, but if is not good, and it must_be_good:
                if comp.isregistered and not comp.isgood and this_job.must_be_good:
                    comp.has_update = True

                    # We make sure that the new_run does not exist:
                    exists, new_comp = this_job.find_computation(keyword=comp.keyword, index=comp.index+1)
                    if exists:        print("Execute_JOB, step 7.3b: Continuation Computation exists")
                    else:
                        if debug > 1: print("Execute_JOB, step 7.3b: Creating new computation to continue job:")

                        ## Creates new computation
                        new_comp = this_job.add_computation(comp.index+1, qc_data, path=comp.path, comp_keyword=comp.keyword, is_update=True, debug=debug)
                        new_comp.qc_data = deepcopy(qc_data)                    ## I need to deepcopy, as it will use a modified version of the qc_data object

                    ## In continuation computations, istate must be the modified, so it continues from the last attempt
                    if  hasattr(comp.qc_data,"fstate"):  
                        new_comp.qc_data._mod_attr("istate",comp.qc_data.fstate)         
                        print("Execute_JOB, step 7.3b: istate of new computation is modified to", new_comp.qc_data.istate)
                    elif hasattr(this_job,"fstate"): 
                        new_comp.qc_data._add_attr("istate",this_job.fstate)         
                        print("Execute_JOB, step 7.3b: istate of new computation is modified to", new_comp.qc_data.istate)
                    else: 
                        print("Execute_JOB, step 7.3b: Could not find valid 'istate' for continuation computation")

        # Updates Job Registry information
        this_job.register(debug=0)
    #recipe.register(debug=0)
    this_branch.register(debug=0)

    if updated: print("saving binary"); save_binary(sys, sys_path)
    return report


