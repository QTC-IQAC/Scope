#!/usr/bin/env python3
import sys
import os
import pwd

from Scope.Classes_Input import *
from Scope.Classes_State import *
from Scope.Classes_SCO import sco_system, crystal
#from Scope.Environment import check_usage, get_queue_and_procs, send_command, set_cluster, set_user
from Scope.Read_Write import load_binary, save_binary

######################
def execute_job(sys_path: str, job_path: str, global_env: object, handle_errors: bool=False, calc_folder: str=None, debug: int=0):
    ### !!!!
    ### It is better that the environment is already loaded, so it is not necessary to load it for every system
    ### !!!!

    if calc_folder is None: calc_folder = sys_path ## Temporary Measure for Unique Ligands

    print("ENTERED EXECUTE JOB with job_path:", job_path)
    report = ''

    if debug > 1: print("")
    if debug > 1: print("----------- NEW JOB ----------")
    #if debug > 1: print("")

    #### 0-Verifies Files
    files = False
    if os.path.isfile(sys_path) and os.path.isfile(job_path): files = True
    if debug > 1: print("EXECUTE_JOB, step 0: sys_path=", sys_path)
    if debug > 1: print("EXECUTE_JOB, step 0: job_path=", job_path)
    if debug > 1: print("EXECUTE_JOB, step 0: files=", files)
    if debug > 1: print("------------------------------")
    if not files: return None

    #### 1-Reads Input Data
    user_environment = set_environment_data(job_path, section="&environment", debug=0)
    options          = set_options_data(job_path, section="&options"        , debug=0)
    job_data         = set_job_data(job_path, section="&job_data"           , debug=0)
    qc_data          = set_qc_data(job_path, section="&qc_data"             , debug=0)

    #### 2a-Enrich Environment with User Choices:
    global_env.read_local_environment(job_path, debug=0)

    #### 2b-Forces some options in case the environment is not that of a computation cluster
    if global_env.management_type == 'None':
        options._mod_attr('want_submit',False)      
        options._mod_attr('overwrite_inputs',False)  
        options._mod_attr('overwrite_logs',False)     

    ##############
    ### SYSTEM ###
    ##############
    ## 3.1-Loads the System
    if debug > 1: print(" ")
    sys = load_binary(sys_path)
    updated = False
    if debug > 1: print(f"EXECUTE_JOB, step 3a: system in {sys_path} loaded")

    ## 3.2-Changes paths if necessary
    #if not os.path.isdir(sys.sys_path):
    #    if debug > 1: print(f"EXECUTE_JOB, step 3b: {sys.sys_path} does not exist")
    if global_env.check_paths(debug=1): 
        if debug > 1: print(f"EXECUTE_JOB, step 3b: global environment found with correct paths. RESETTING")
        updated = sys.reset_paths(global_env, debug=0)
        if updated and debug > 1: print(f"EXECUTE_JOB, step 3b: system paths reset")

    ##############
    ### BRANCH ###
    ##############
    if sys.type == "sco_system":
        exists, this_branch = sys.find_branch(job_data.branch, debug=0)
        if not exists: this_branch = sys.add_branch(job_data.branch, job_data.target, debug=debug); updated = True
    elif sys.type == "Ligand":
        assert job_data.target == 'self'
        from Scope.Gmol_ops import find_branch_gmol, add_branch_gmol
        exists, this_branch = find_branch_gmol(sys, job_data.branch, debug=0)
        if not exists: this_branch = add_branch_gmol(sys, job_data.branch, calc_folder, debug=debug); updated = True
    if debug > 1: print("EXECUTE_JOB, step 4: branch loaded")

    ##############
    ### RECIPE ###
    ##############
    for idx, recipe in enumerate(this_branch.recipes):

        ###########
        ### JOB ###
        ###########
        ## 5-Finds the job. If it does not exist, it is created in 5.1
        exists, this_job = recipe.find_job(job_data=job_data, debug=0)

        ## 5.1 If job does not exist, creates the job
        if not exists: this_job = recipe.add_job(job_data); updated = True
        ## 5.2 If job exists, it checks for input changes (also in qc_data for the computations)
        else:
            this_job.check_input(job_path=job_path, debug=debug)

        if debug > 1: print("---------------------------------------------------")
        if debug > 1: print("EXECUTE_JOB, step 5: job loaded")
        if debug > 1: print("---------------------------------------------------")

        ## 6-Checks that all requisites and constrains of the job are fulfilled
        cancontinue = this_job.check_requisites(debug=0)
        if not cancontinue:
            if debug > 1:   print("EXECUTE_JOB, step 6: requisites NOT met or job already run!")
            continue        # I know if might seem misleading. Here, "continue" means "skip this one"
        else:
            if debug > 1:   print("EXECUTE_JOB, step 6: requisites fulfilled")

        ####################
        ### COMPUTATIONS ###
        ####################
        ## 7-Sets the computation(s), meaning that it will check if they exist, and if not, it creates them
        this_job.set_computations_from_setup(qc_data, debug=debug)
        if debug > 1: print("EXECUTE_JOB, step 7: computations set:")

        for jdx, comp in enumerate(this_job.computations):

            #if comp.has_update and comp.isregistered: continue # Skip jobs with update (i.e. with other related computations with higher run_number)
            if debug > 1: print("")
            if debug > 1: print("########################################################################")
            if debug > 1: print(f"    {sys.refcode} -> {this_job._recipe.subject.spin} -> {this_job.keyword} -> {comp.run_number}")
            if debug > 1: print("########################################################################")
            if debug > 1: print(f"EXECUTE_JOB, step 7.0: evaluating job, and computation with indices: {recipe.jobs.index(this_job)+1}/{len(recipe.jobs)}, {jdx+1}/{len(this_job.computations)}")

            ## 8.1-Checks files
            comp.check_files()
            if debug > 1: print("EXECUTE_JOB, step 7.1: doing computation with keyword and run_number:", comp.keyword, comp.run_number)
            if debug > 1: print("EXECUTE_JOB, step 7.1: out_file:", comp.out_path)
            if debug > 1: print("EXECUTE_JOB, step 7.1: is_update:", comp.is_update)
            if debug > 1: print("EXECUTE_JOB, step 7.1: checking files [inp, out, sub]:", comp.input_exists, comp.output_exists, comp.subfile_exists)
            #if debug > 1: print("-----------------------------------------------------------------------------------")

            ## 8.2-Evaluates Submission
            if not comp.output_exists: # and comp.input_exists:
                if options.want_submit:
                    comp.check_submission_status(global_env, debug=debug)
                    if debug > 1: print("EXECUTE_JOB, step 7.3a: checking submission status: isrunning=",comp.isrunning)
                    if debug > 1: print("EXECUTE_JOB, step 7.3a: initial state is", comp.qc_data.istate)
                    if debug > 1: print("EXECUTE_JOB, step 7.3a: is_update:", comp.is_update)
                    if not comp.isrunning:        comp.run(global_env, options, debug=debug); updated = True
                else: 
                    if debug > 1: print("EXECUTE_JOB, step 7.3a: want_submit is False")
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
                    updated = True

                ## 8.4 sets continuation computations
                if comp.status == 'scf_convergence':
                #if comp.isregistered and comp.status == 'scf_convergence':
                    new_comp = this_job.set_continuation_computation(comp, "scf", debug=debug)
                    #if new_comp.run_number >= 10: report += f"Investigate {new_comp.out_path} \n"
                elif not comp.isgood and this_job.must_be_good:
                #elif comp.isregistered and not comp.isgood and this_job.must_be_good:
                    new_comp = this_job.set_continuation_computation(comp, "opt", debug=debug)
                    #if new_comp.run_number >= 10: report += f"Investigate {new_comp.out_path} \n"


                print("HERE WITH:", comp.out_path, comp.isregistered, comp._job.keyword, hasattr(comp.qc_data,"fstate") or hasattr(comp._job,"fstate"))
                if comp.isregistered and "freq" in comp._job.keyword and (hasattr(comp.qc_data,"fstate") or hasattr(comp._job,"fstate")):
                    if hasattr(comp.qc_data,"fstate"): fstate = comp.qc_data.fstate
                    else:                              fstate = comp._job.fstate
                    exists, state = find_state(comp._job._recipe.subject, fstate)
                    print("State",fstate,"exist=", exists)
                    if exists and hasattr(state,"VNMs"):
                        if not hasattr(state.VNMs,"xs"): 
                            print("RE-REGISTERING:", comp.out_path)
                            worked = comp.register(debug=debug)
#                            for v in state.VNMs: 
#                                print(v.index, v.xs[0])
                    else:
                        print("State",fstate,"does not exist")

        # Updates Job Registry information
        this_job.register(debug=0)
    #recipe.register(debug=0)
    this_branch.register(debug=0)

    if updated: print("saving binary"); save_binary(sys, sys_path)
    return report


