#!/usr/bin/env python3
import os
from scope.classes_environment  import *
from scope.classes_input        import *
from scope.classes_state        import *
from scope.classes_system       import *
from scope.read_write           import load_binary

######################
def run_task(sys_path: str, inp_paths: list, global_env: str | object, handle_errors: bool=False, debug: int=0):
    """
    Runs a SCOPE Input (defined in the inp_paths file) on a SCOPE system (in the sys_path). 
    The Configuration of the Computer is read from the GLOBAL_ENVIRONMENT, which must be configured before and given as a binary.
    This function performs the following steps:
    - Verifies the existence of the system and input files
    - Reads input data from the input file, including environment, options, job data, and QC_data.
    - Loads the system object from the binary file and updates paths if necessary. Handy if registration and execution of tasks are performed in different computers
    - Finds or creates the required branch and workflow in the system.
    - Finds or creates the job, and checks for input changes.
    - Validates job requisites and continues only if they are fulfilled.
    - Sets up computations for the job, checks file existence, and handles submission.
    - Registers computations, handles errors, and manages continuation or repetitive computations as needed.
    - Updates registry information and saves the system binary if changes occurred.
    ----------
    Parameters
    ----------
    sys_path : str
        Path to the system binary file
    inp_paths : list
        List of Paths to the job files. 
    global_env : str or object
        Global environment as a path or directly the object.
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

    ##############
    ### HEADER ###
    ##############
    print(f"##################")
    print(f"- Starting SCOPE -")
    print(f"##################")
    print(f"")

    #############################
    ### STEP 0: Loading Files ###
    #############################
    ## 0.1 Verifies that all files exist
    if not verify_files(sys_path, global_env, inp_paths, debug=debug): raise ValueError("RUN_TASK: File verification failed")

    ## 0.2-Checks status of branch in system. This is to avoid loading the system if calculations for that branch are already finished
    if not verify_status(sys_path, inp_paths, debug=debug): 
        if debug > 0: print("RUN_TASK: All branches are either terminated or finished. No need to load system")
        return report 

    ## 0.3-Loads the System
    sys = load_binary(sys_path)
    if debug > 0: print(f"RUN_TASK, step 0.3: System {sys.name} loaded from {sys_path}") 
    updated = False

    ## 0.4-Loads the Global Environment, if not provided already as an object 
    if isinstance(global_env, str): 
        global_env = load_binary(global_env) 
        if debug > 0: print(f"RUN_TASK, step 0.4: Environment Loaded")

    ### 0.5 Verifies that all relevant paths inside system and environment.
    if not verify_paths(global_env, sys, debug=debug): raise ValueError("RUN_TASK: Path verification failed")

    ###################
    #### MAIN LOOP ####
    ###################
    for inp_idx, inp_path in enumerate(inp_paths):

        if debug > 0: print("")
        if debug > 0: print(f"#####################")
        if debug > 0: print(f"- Starting JOB {inp_idx+1}/{len(inp_paths)} -")
        if debug > 0: print(f"#####################")
        
        print("ABSPATH:", os.path.abspath('.'))

        #### Step 1.1-Reads Input Data
        user_environment = set_environment_data(inp_path, section="&environment", debug=0)  ## For completeness, but env data is read in global_env.read_job_specs below 
        options          = set_options_data(inp_path, section="&options"        , debug=0)
        job_data         = set_job_data(inp_path, section="&job_data"           , debug=0)
        qc_data          = set_qc_data(inp_path, section="&qc_data"             , debug=0)

        #### Step 1.2-Loads Environment and Enriches with User Choices
        global_env.read_job_specs(inp_path, debug=0)

        #### Step 1.3-Forces some options in case the environment is not that of a computation cluster
        if global_env.scheduler == 'local':
            print("WARNING!!! No Queue Management has been detected in cluster, disabling submission")
            options._mod_attr('want_submit',False)      
            options._mod_attr('overwrite_inputs',False)  
            options._mod_attr('overwrite_logs',False)     

        ########################
        ### STEP 1.4: BRANCH ###
        ########################
        exists, this_branch        = sys.find_branch(job_data.branch, debug=0)
        if not exists: this_branch = sys.add_branch(job_data.branch, debug=debug); updated = True

        ##########################
        ### STEP 1.5: WORKFLOW ###
        ##########################
        if debug > 0:   print(f"RUN_TASK, step 1.5: Evaluating list of Workflows in input: {job_data.workflow=}")
        job_data.workflow = interpret_workflow_from_user(sys, job_data, debug=debug)
  
        for wrk in job_data.workflow if isinstance(job_data.workflow, list) else list([job_data.workflow]):   ### Works when job_data.workflow is a str or a list

            exists, this_workflow             = this_branch.find_workflow(wrk)
            if not exists: this_workflow      = this_branch.add_workflow(wrk); updated = True
            if this_workflow is None: 
                raise ValueError(f"RUN_TASK, step 2.0: WORKFLOW could not be created. Make sure a Source with name {wrk} exist in System")
            if debug == 1:   print(f"RUN_TASK, step 2.0: Evaluating WORKFLOW with name={this_workflow.name}")
            elif debug > 1:  print(f"RUN_TASK, step 2.0: Evaluating WORKFLOW \n {this_workflow}")

            ###################
            ### STEP 2: JOB ###
            ###################
            ## 2.1 Finds or creates the job.
            print(f"RUN_TASK, step 2.1: Evaluating JOB with keyword={job_data.keyword}")
            exists, this_job = this_workflow.find_job(job_data=job_data, debug=debug)
            if not exists: this_job = this_workflow.add_job(job_data); updated = True

            print(f"RUN_TASK, step 2.2: Checking its job_data")
            ## 2.2 If job exists, it checks for input changes (also in qc_data for the computations)
            print("ABSPATH:", os.path.abspath('.'))
            if exists: this_job.check_job_data(inp_path=inp_path, debug=debug)

            ## 2.3-Checks that all requisites and constrains of the job are fulfilled
            cancontinue = this_job.check_requisites(debug=debug)
            if not cancontinue:
                if debug > 0:   
                    print(f"RUN_TASK, step 2.3: requisites NOT met or job already run:")
                    print(f"RUN_TASK, step 2.3: requisites: {this_job.requisites}")
                    print(f"RUN_TASK, step 2.3: constrains: {this_job.constrains}")
                continue        # I know if might seem misleading. Here, "continue" means "skip this one"
            print("RUN_TASK, step 2.3: Requisites passed")

            ############################
            ### STEP 3: COMPUTATIONS ###
            ############################
            ## 3-Sets the computation(s), meaning that it will check if they exist, and if not, it creates them
            if debug > 0: print(f"RUN_TASK, step 2.4: setting job computations with job_setup={this_job.job_setup}")
            this_job.set_computations_from_setup(qc_data, debug=debug)
            for jdx, comp in enumerate(this_job.computations):

                #if comp.has_update and comp.isregistered: continue # Skip jobs with update (i.e. with other related computations with higher run_number)
                if debug > 0: print(f"--------------------------------------------------------------------------")
                if debug > 0: print(f" {sys.name} -> {this_branch.name} -> {this_workflow.name} -> {this_job.keyword} -> Step: {comp.step} -> Run: {comp.run_number}")
                if debug > 0: print(f"--------------------------------------------------------------------------")
                if debug > 0: print(f"RUN_TASK, step 3.0: evaluating job, and computation with indices: {this_workflow.jobs.index(this_job)+1}/{len(this_workflow.jobs)}, {jdx+1}/{len(this_job.computations)}")

                ## Step 3.1-Checks files and updates
                qc_has_updated = comp.check_qc_data(inp_path=inp_path, debug=debug)  ## Checks wether the user has updated the qc_data
                comp.check_updates()                                                 ## Checks for not-registered update computations 
                comp.check_files()

                if debug > 0: print("RUN_TASK, step 3.1: doing computation with keyword and run_number:", comp.keyword, comp.run_number)
                if debug > 0: print("RUN_TASK, step 3.1: out_file:", comp.out_path)
                if debug > 0: print("RUN_TASK, step 3.1: is_update:", comp.is_update)
                if debug > 0: print("RUN_TASK, step 3.1: has_update:", comp.has_update)
                if debug > 0: print("RUN_TASK, step 3.1: checking files [inp, out, sub]:", comp.input_exists, comp.output_exists, comp.subfile_exists)
                if debug > 0: print("RUN_TASK, step 3.1: qc has been updated", qc_has_updated)
                #if debug > 1: print("-----------------------------------------------------------------------------------")

                ## Step 3.2a-Evaluates Submission
                if not comp.output_exists:
                    if options.want_submit and not comp.has_update:
                        comp.check_submission_status(global_env, debug=debug)
                        if debug > 0: print("RUN_TASK, step 3.2a: checking submission status: isrunning=",comp.isrunning)
                        if debug > 0: print("RUN_TASK, step 3.2a: ignore_submitted=", options.ignore_submitted)
                        if debug > 0: print("RUN_TASK, step 3.2a: initial state is", comp.qc_data.istate)
                        if not comp.isrunning or options.ignore_submitted:       ## options.ignore_submitted will also be checked in comp.run 
                            comp.check_qc_data(inp_path=inp_path, debug=debug)
                            comp.run(global_env, options, debug=debug); updated = True
                    else: 
                        if debug > 0: print("RUN_TASK, step 3.2a: want_submit is False or comp.has_update")

                ## Step 3.2b-Warns if output exists but not input
                elif comp.output_exists and not comp.input_exists:  
                    report += f"Investigate {comp.out_path} \n"
                    print(f"Investigate {comp.out_path}")

                ############################
                ### STEP 4: REGISTRATION ###
                ############################
                elif comp.output_exists and comp.input_exists:
                    ## 4.1-If output exists, and is not registered, it does it
                    ## This means that the output is read, parsed. The parsed data depends on the type of computation
                    if comp.isregistered:
                        if debug > 0: print(f"RUN_TASK, step 4: job already registered")
                    else:
                        if debug > 0: print(f"RUN_TASK, step 4: registration")
                        worked = comp.register(debug=debug) 

                        if debug > 0: print(f"RUN_TASK, step 4.1: registration {worked=}")
                        if debug > 0: print(f"RUN_TASK, step 4.1: {comp.has_update=}")
                        if debug > 0: print(f"RUN_TASK, step 4.1: {options.overwrite_inputs=}")

                        ## 4.2 sets continuation computations. These are added to JOB object 
                        if comp.status == 'no_scf_convergence': 
                            if debug > 0: print(f"RUN_TASK, step 4.2: setting continuation computation with typ=scf")  
                            new_comp = this_job.set_continuation_computation(comp, "scf", debug=debug)
                            if new_comp.run_number >= 10: report += f"Investigate {new_comp.out_path} \n"
                        elif not comp.isgood and this_job.must_be_good:
                            if debug > 0: print(f"RUN_TASK, step 4.2: setting continuation computation with typ=opt")  
                            new_comp = this_job.set_continuation_computation(comp, "opt", debug=debug)
                            if new_comp.run_number >= 10: report += f"Investigate {new_comp.out_path} \n"

                        ## 4.3 Cases meant to be repetitive. Next step added to JOB object 
                        if this_job.job_setup == "rep_opt" and comp.isgood:
                            from scope.other import check_convergence
 
                            ## 4.3.1 Collects energies from state in this step
                            if hasattr(comp.qc_data,"fstate"): fstate_name = comp.qc_data.fstate
                            else:                              fstate_name = comp._job.fstate
                            exists, fstate = comp.source.find_state(fstate_name, debug=debug)         ## This one MUST exist
                            assert exists
                            if comp.step > len(this_job.energies): 
                                this_job.energies = np.append(this_job.energies,int(0))
                            this_job.energies[int(comp.step)-1] = fstate.results['energy'].value

                            ## 4.3.2 Checks the energy convergence
                            this_job.isconverged = False
                            if comp.step > 0: this_job.isconverged = check_convergence(this_job.energies, comp.step-1, this_job.job_data.energy_thres)

                            ## 4.3.3 Continutes if not converged and below max_steps
                            if this_job.isconverged: 
                                print(f"RUN_TASK, step 4.3: repetitive opt reached convergence")  
                            elif not this_job.isconverged and comp.step <= this_job.job_data.max_steps:
                                new_comp = this_job.set_continuation_computation(comp, "rep_opt", debug=debug)     
                            else:
                                print(f"RUN_TASK, step 4.3: maximum steps reached without convergence")  

                        # Checks for common stupid errors and handles files
                        if debug > 0: print(f"RUN_TASK, step 4.4: registration {worked=}")  
                        if not worked:
                            if handle_errors:          # ...takes default action 
                                comp.read_lines()
                                if len(comp.output_lines) > 0: comp.store(debug=debug)  # Creates Copy of output 
                                else:                          this_job.remove_computation(comp_index=comp.index)
                                report += f"Errors handled for {comp.out_path}. Please Re-Submit \n"
                                print(f"RUN_TASK: Errors handled for {comp.out_path}")
                            else:
                                report += f"Error registering {comp.out_path} . Please Re-Submit \n"
                                print(f"RUN_TASK: Error registering {comp.out_path}")
                        updated = True
                        if debug > 0: print("")

            # Updates Job Registry information. Necessary to update class.isregistered attribute, but not very important
            this_job.register(debug=0)
        this_workflow.register(debug=0)
        this_branch.register(debug=0)

        if updated: 
            if debug > 0: print(f"Saving System {sys.name} to {sys.system_file}")
            sys.save()
    return report

#######################
## Related Functions ##
#######################
def get_status(sys_folder: str, branch_name: str, debug: int=0):
    ##########################################################
    ## Function to determine whether the computations of    ##
    ## a Branch for a given system have already finished    ##
    ## without having to load the system file, to save time ##
    ##########################################################

    ## If system file does not exist
    if not os.path.isdir(f"{sys_folder}"):
        print(f"RUN_TASK.GET_STATUS: System path {sys_folder} does not exist")
        return "absent"

    ## Otherwise.. it looks for files with standardized names (see setup functions in scope.classes_workflow)
    if os.path.isfile(f"{sys_folder}TERMINATED"):                return "terminated"
    if os.path.isfile(f"{sys_folder}{branch_name}_FINISHED"):    return "finished"  
    return "active"

######
def verify_files(sys_path: str, global_env: str | object, inp_paths: list, debug=0) -> bool:
    sys_file  = os.path.isfile(sys_path)
    inp_files = all(os.path.isfile(ip) for ip in inp_paths)
    if isinstance(global_env, str):      env_file = os.path.isfile(global_env)
    elif isinstance(global_env, object): env_file = os.path.isfile(global_env.filepath)
    if sys_file and inp_files and env_file: return True
    else:
        if debug > 0: 
            if not sys_file: print(f"RUN_TASK.VERIFY_FILES: System file {sys_path} does not exist")
            if not inp_files: 
                for ip in inp_paths:
                    if not os.path.isfile(ip): print(f"RUN_TASK.VERIFY_FILES: Job file {ip} does not exist")
            if not env_file:
                if isinstance(global_env, str):      print(f"RUN_TASK.VERIFY_FILES: Environment file {global_env} does not exist")
                elif isinstance(global_env, object): print(f"RUN_TASK.VERIFY_FILES: Environment file {global_env.filepath} does not exist")
    return False

######
def verify_paths(global_env: object, sys: object, debug: int=0) -> bool: 
    ## Check Env Paths
    if global_env.check_paths(debug=debug): 
        if debug > 0: print(f"RUN_TASK.VERIFY_PATHS: Environment Paths were Located")
    else:
        raise    ValueError(f"RUN_TASK.VERIFY_PATHS: Environment Paths were NOT Located. Please check paths in environment at {global_env.filepath}")

    ## Check Sys Paths
    if sys.check_paths(debug=debug):
        if debug > 0: print(f"RUN_TASK.VERIFY_PATHS: System Paths were Located")
    else:  
        if debug > 0: print(f"RUN_TASK.VERIFY_PATHS: System Paths were NOT Located. Resetting from Environment")
        update = sys.set_paths_from_environment(global_env, debug=debug)  ## Checks and resets paths if necessary
        if not update: raise ValueError(f"RUN_TASK.VERIFY_PATHS: System Paths could not be RESET. Please check paths in system at {sys.sys_file}")
        else: print(f"RUN_TASK.VERIFY_PATHS: system paths reset")
    return True

######
def verify_status(sys_path: str, inp_paths: list, debug: int=0) -> bool:

    sys_folder = os.path.dirname(sys_path)
    if debug > 0: print(f"RUN_TASK.VERIFY_STATUS: searching status files in {sys_folder}")
    if debug > 0: print(f"RUN_TASK.VERIFY_STATUS: found {os.listdir(sys_folder)} files in folder")

    ## Collects all branches requested in the different job files
    branch_names = []
    for inp_path in inp_paths:
        job_data = set_job_data(inp_path, section="&job_data", debug=0)
        branch_names.append(job_data.branch)
    if debug > 0: print(f"RUN_TASK.VERIFY_STATUS: found {branch_names=}")

    ## Checks Status of all branches
    all_status = []
    for br in branch_names:
        all_status.append(get_status(sys_folder, br, debug=debug))
    if debug > 0: print(f"RUN_TASK.VERIFY_STATUS: found {all_status=}")

    ## If any is active, then proceed (i.e. whether the system must be loaded) is true
    if any(s == 'active' for s in all_status): proceed = True 
    else:                                      proceed = False

    ## Apart from that, informs if any other status is found
    if debug > 0: 
        if any(all_status) == 'terminated':
            tidx = [i for i, status in enumerate(all_status) if status == 'terminated']
            for t in tidx:
                print(f"RUN_TASK.VERIFY_STATUS: Branch {branch_names[t]} is TERMINATED. Job will be skipped.") 
                print(f"RUN_TASK.VERIFY_STATUS: If you wish to activate it, remove file 'TERMINATED' in {sys_folder}/{branch_names[t]}")
        if any(all_status) == 'finished':
            tidx = [i for i, status in enumerate(all_status) if status == 'finished']
            for t in tidx:
                print(f"RUN_TASK.VERIFY_STATUS: Branch {branch_names[t]} is FINISHED. Job will be skipped.") 
                print(f"RUN_TASK.VERIFY_STATUS: If you wish to activate it, remove file '{branch_names[t]}_FINISHED' in {sys_folder}/{branch_names[t]}")
    if debug > 0: print(f"RUN_TASK.VERIFY_STATUS: returning {proceed=}")

    return proceed
