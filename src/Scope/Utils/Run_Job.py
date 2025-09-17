#!/usr/bin/env python3
import os
import sys
import pwd
from Scope.Classes_Environment  import environment, read_job_specs
from Scope.Classes_Input        import *
from Scope.Classes_State        import *
from Scope.Read_Write           import load_binary

######################
def run_job(sys_path: str, job_path: str, global_env_path: str, handle_errors: bool=False, debug: int=0):
    """
    Runs a SCOPE Task (defined in the job_path file) on a SCOPE system (in the sys_path). 
    The Configuration of the Computer is read from the GLOBAL_ENVIRONMENT, which must be configured before and given as a binary.
    This function performs the following steps:
    1. Verifies the existence of the system and job files.
    2. Reads input data from the job file, including environment, options, job data, and QC_data.
    3. Updates the global environment with user-specific choices.
    4. Adjusts options if no queue management is detected.
    5. Loads the system object from the binary file and updates paths if necessary. Handy if registration and execution of tasks are performed in different computers
    6. Finds or creates the required branch and recipe in the system.
    7. Finds or creates the job within the recipe, and checks for input changes.
    8. Validates job requisites and continues only if they are fulfilled.
    9. Sets up computations for the job, checks file existence, and handles submission.
    10. Registers computations, handles errors, and manages continuation or repetitive computations as needed.
    11. Updates registry information and saves the system binary if changes occurred.
    ----------
    Parameters
    ----------
    sys_path : str
        Path to the system binary file.
    job_path : str
        Path to the job configuration file.
    global_env : object
        Global environment object, expected to be pre-loaded.
    handle_errors : bool, optional
        If True, handles errors encountered during registration (default is False).
    debug : int, optional
        Debug level for verbose output (default is 0).
    Returns
    -------
    report : str or None
        A report string summarizing any issues or actions taken during run,
        or None if the required files do not exist.
    """

    report = ''

    if debug > 1: print("")
    if debug > 1: print("----------- NEW JOB ----------")

    #### 0-Verifies Files
    files = False
    if os.path.isfile(sys_path) and os.path.isfile(job_path) and os.path.isfile(global_env_path): files = True
    if debug > 1: print("RUN_JOB, step 0: sys_path=", sys_path)
    if debug > 1: print("RUN_JOB, step 0: job_path=", job_path)
    if debug > 1: print("RUN_JOB, step 0: files=", files)
    if debug > 1: print("------------------------------")
    if not files: return None

    #### 1a-Reads Input Data
    user_environment = set_environment_data(job_path, section="&environment", debug=0)  ## For completeness, but env data is read in global_env.read_job_specs below 
    options          = set_options_data(job_path, section="&options"        , debug=0)
    job_data         = set_job_data(job_path, section="&job_data"           , debug=0)
    qc_data          = set_qc_data(job_path, section="&qc_data"             , debug=0)

    #### 2a-Checks status of system and branch. This is to avoid loading the system if calculations for that branch are already finished
    status = get_status(sys_path, job_data.branch, debug=debug)
    if status != 'active':
        if debug > 0: 
            print(f"RUN_JOB, step 2a: status of {sys_path} for branch {job_data.branch} is {status}. Skipping")
            if   status == "terminated": print(f"RUN_JOB, step 2a: I you wish to activate it, remove file 'TERMINATED' in {sys_path}")
            elif status == "finished":   print(f"RUN_JOB, step 2a: I you wish to activate it, remove file '{job_data.branch.name}_finished' in {sys_path}")
        return None

    #### 2b-Load Environment and Enriches with User Choices:
    global_env      = load_binary(global_env_path)
    global_env.read_job_specs(job_path, debug=0)

    #### 2c-Forces some options in case the environment is not that of a computation cluster
    if global_env.management_type == 'None':
        print("WARNING!!! No Queue Management has been detected in cluster, disabling submission")
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
    if debug > 1: print(f"RUN_JOB, step 3a: system in {sys_path} loaded")

    ## 3.2-Changes paths if necessary
    if global_env.check_paths(debug=1): 
        if debug > 1: print(f"RUN_JOB, step 3b: global environment found with correct paths")
        try: 
            updated = sys.reset_paths(global_env, debug=0)
            if updated and debug > 1: print(f"RUN_JOB, step 3b: system paths reset")
        except Exception as exc: 
            pass

    ##############
    ### BRANCH ###
    ##############
    exists, this_branch        = sys.find_branch(job_data.branch, debug=0)
    if not exists: this_branch = sys.add_branch(job_data.branch, debug=debug); updated = True

    ##############
    ### RECIPE ###
    ##############
    for rec in job_data.recipe if isinstance(job_data.recipe, list) else list([job_data.recipe]):   ### Works when job_data.recipe is a str or a list

        exists, recipe             = this_branch.find_recipe(rec)
        if not exists: recipe      = this_branch.add_recipe(rec); updated = True

        ###########
        ### JOB ###
        ###########
        ## 5.1 Finds or creates the job.
        exists, this_job = recipe.find_job(job_data=job_data, debug=0)
        if not exists: this_job = recipe.add_job(job_data); updated = True

        ## 5.2 If job exists, it checks for input changes (also in qc_data for the computations)
        if exists: this_job.check_job_data(job_path=job_path, debug=debug)

        if debug > 1: print("---------------------------------------------------")
        if debug > 1: print(f"RUN_JOB, step 5: job {this_job.name} loaded")
        if debug > 1: print("---------------------------------------------------")

        ## 6-Checks that all requisites and constrains of the job are fulfilled
        cancontinue = this_job.check_requisites(debug=debug)
        if not cancontinue:
            if debug > 1:   
                print("RUN_JOB, step 6: requisites NOT met or job already run. Printing job")
                print(this_job)
            continue        # I know if might seem misleading. Here, "continue" means "skip this one"
        else:
            if debug > 1:   print("RUN_JOB, step 6: requisites fulfilled")

        ####################
        ### COMPUTATIONS ###
        ####################
        ## 7-Sets the computation(s), meaning that it will check if they exist, and if not, it creates them
        this_job.set_computations_from_setup(qc_data, debug=debug)
        if debug > 1: print("RUN_JOB, step 7: computations set:")

        for jdx, comp in enumerate(this_job.computations):

            #if comp.has_update and comp.isregistered: continue # Skip jobs with update (i.e. with other related computations with higher run_number)
            if debug > 1: print("")
            if debug > 1: print("########################################################################")
            if debug > 1: print(f"    {sys.name} -> {this_job._recipe.name} -> {this_job.name} -> {comp.step} -> {comp.run_number}")
            if debug > 1: print("########################################################################")
            if debug > 1: print(f"RUN_JOB, step 7.0: evaluating job, and computation with indices: {recipe.jobs.index(this_job)+1}/{len(recipe.jobs)}, {jdx+1}/{len(this_job.computations)}")

            ## 7.0-Checks files and updates
            qc_has_updated = comp.check_qc_data(job_path=job_path, debug=debug)  ## Checks wether the user has updated the qc_data
            comp.check_updates()                                                 ## Checks for not-registered update computations 
            comp.check_files()

            if debug > 1: print("RUN_JOB, step 7.1: doing computation with keyword and run_number:", comp.keyword, comp.run_number)
            if debug > 1: print("RUN_JOB, step 7.1: out_file:", comp.out_path)
            if debug > 1: print("RUN_JOB, step 7.1: is_update:", comp.is_update)
            if debug > 1: print("RUN_JOB, step 7.1: has_update:", comp.has_update)
            if debug > 1: print("RUN_JOB, step 7.1: checking files [inp, out, sub]:", comp.input_exists, comp.output_exists, comp.subfile_exists)
            if debug > 1: print("RUN_JOB, step 7.1: qc has been updated", qc_has_updated)
            #if debug > 1: print("-----------------------------------------------------------------------------------")

            ## 7.2a-Evaluates Submission
            if not comp.output_exists:
                if options.want_submit and not comp.has_update:
                    comp.check_submission_status(global_env, debug=debug)
                    if debug > 1: print("RUN_JOB, step 7.2a: checking submission status: isrunning=",comp.isrunning)
                    if debug > 1: print("RUN_JOB, step 7.2a: ignore_submitted=", options.ignore_submitted)
                    if debug > 1: print("RUN_JOB, step 7.2a: initial state is", comp.qc_data.istate)
                    if not comp.isrunning or options.ignore_submitted:       ## options.ignore_submitted will also be checked in comp.run 
                        comp.check_qc_data(job_path=job_path, debug=debug)
                        comp.run(global_env, options, debug=debug); updated = True
                else: 
                    if debug > 1: print("RUN_JOB, step 7.2a: want_submit is False or comp.has_update")

            ## 7.2b-Warns if output exists but not input
            elif comp.output_exists and not comp.input_exists:  
                report += f"Investigate {comp.out_path} \n"
                print(f"Investigate {comp.out_path}")

            ## Step 8: registration
            elif comp.output_exists and comp.input_exists:
                ## 8.1-If output exists, and is not registered, it does it
                ## This means that the output is read, parsed. The parsed data depends on the type of computation
                if not comp.isregistered:
                    if debug > 1: print(f"RUN_JOB, step 8: registration")
                    worked = comp.register(debug=debug) 

                    if debug > 1: print(f"RUN_JOB, step 8.1: registration {worked=}")
                    if debug > 1: print(f"RUN_JOB, step 8.1: {comp.has_update=}")
                    if debug > 1: print(f"RUN_JOB, step 8.1: {options.overwrite_inputs=}")

                    ## 8.2 sets continuation computations. These are added to JOB object 
                    if comp.status == 'no_scf_convergence': 
                        if debug > 1: print(f"RUN_JOB, step 8.2: setting continuation computation with typ=scf")  
                        new_comp = this_job.set_continuation_computation(comp, "scf", debug=debug)
                        if new_comp.run_number >= 10: report += f"Investigate {new_comp.out_path} \n"
                    elif not comp.isgood and this_job.must_be_good:
                        if debug > 1: print(f"RUN_JOB, step 8.2: setting continuation computation with typ=opt")  
                        new_comp = this_job.set_continuation_computation(comp, "opt", debug=debug)
                        if new_comp.run_number >= 10: report += f"Investigate {new_comp.out_path} \n"
                    
                    ## 8.3 Cases meant to be repetitive. Next step added to JOB object 
                    if this_job.setup == "rep_opt" and comp.isgood:
                        from Scope.Other import check_convergence
 
                        ## 8.3.1 Collects energies from state in this step
                        if hasattr(comp.qc_data,"fstate"): fstate = comp.qc_data.fstate
                        else:                              fstate = comp._job.fstate
                        exists, state = find_state(comp._job._recipe.source, fstate)
                        if comp.step > len(this_job.energies): 
                            this_job.energies = np.append(this_job.energies,int(0))
                        this_job.energies[int(comp.step)-1] = state.results['energy'].value

                        ## 8.3.2 Checks the energy convergence
                        this_job.isconverged = False
                        if comp.step > 1: this_job.isconverged = check_convergence(this_job.energies, comp.step-1, this_job.job_data.energy_thres)

                        ## 8.3.3 Continutes if not converged and below max_steps
                        if this_job.isconverged: 
                            print(f"RUN_JOB, step 8.3: repetitive opt reached convergence")  
                        elif not this_job.isconverged and comp.step <= this_job.job_data.max_steps:
                            new_comp = this_job.set_continuation_computation(comp, "rep_opt", debug=debug)     
                        else:
                            print(f"RUN_JOB, step 8.3: maximum steps reached without convergence")  

                    # Checks for common stupid errors and handles files
                    if not worked:
                        if handle_errors:          # ...takes default action 
                            comp.read_lines()
                            if len(comp.output_lines) > 0: comp.store(debug=debug)  # Creates Copy of output 
                            else:                          this_job.remove_computation(comp_index=comp.index)
                            report += f"Errors handled for {comp.out_path}. Please Re-Submit \n"
                            print(f"RUN_JOB: Errors handled for {comp.out_path}")
                        else:
                            report += f"Error registering {comp.out_path} . Please Re-Submit \n"
                            print(f"RUN_JOB: Error registering {comp.out_path}")
                    updated = True

                ### 
                ### Re-reads eigenvectors of a frequency computation 
                ### Not sure I need this
                ### 
                #if comp.isregistered and "freq" in comp._job.keyword and (hasattr(comp.qc_data,"fstate") or hasattr(comp._job,"fstate")):
                #    if hasattr(comp.qc_data,"fstate"): fstate = comp.qc_data.fstate
                #    else:                              fstate = comp._job.fstate
                #    exists, state = find_state(comp._job._recipe.source, fstate)
                #    print("State",fstate,"exist=", exists)
                #    if exists and hasattr(state,"VNMs"):
                #        if not hasattr(state.VNMs,"xs"): 
                #            print("RE-REGISTERING:", comp.out_path)
                #            worked = comp.register(debug=debug)
                #    else:
                #        print("State",fstate,"does not exist")

        # Updates Job Registry information
        this_job.register(debug=0)
    #recipe.register(debug=0)
    this_branch.register(debug=0)

    if updated: print("Saving System"); sys.save()
    return report

def get_status(sys_path: str, branch_name, debug: int=0):
    ##########################################################
    ## Function to determine whether the computations of    ##
    ## a Branch for a given system have already finished    ##
    ## without having to load the system file, to save time ##
    ##########################################################
    if sys_path[-1] != '/': sys_path += '/'

    ## If system path does not exist
    if not os.path.isdir(f"{sys_path}"):
        print(f"GET_STATUS: System path {sys_path} does not exist")
        return "absent"

    ## Otherwise.. it looks for files with standardized names (see setup functions in Scope.Workflow.Branch)
    if os.path.isfile(f"{sys_path}TERMINATED"):                return "terminated"
    if os.path.isfile(f"{sys_path}{branch_name}_FINISHED"):    return "finished"  
    return "active"
