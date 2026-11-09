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
    Run one or more SCOPE input files on a stored system.

    Parameters:
        sys_path (str):                Path to the system binary.
        inp_paths (list):              Input-job file paths.
        global_env (str | object):     Environment path or loaded environment object.
        handle_errors (bool):          Retry startup aborts and apply automatic file/error handling.
        debug (int):                   Verbosity level.

    Returns:
        str: Report string summarizing unresolved issues that require manual attention.

    Stored branches whose computation folder is unavailable are skipped and
    reported once. Their folders and workflow history are not recreated or reset.
    Save after a task only when System.should_save marks pending changes.
    """

    report = ''
    unavailable_branches = set()

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
        
        #### Step 1.1-Reads Input Data
        user_environment, options, job_data, qc_data = set_input_data(inp_path, debug=0)

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
        if not exists: 
            if debug > 0:   print(f"RUN_TASK, step 1.4: Branch not found. Creating it")
            this_branch = sys.add_branch(job_data.branch, debug=debug); sys.should_save = True

        if not os.path.isdir(this_branch.path):
            if this_branch.name not in unavailable_branches:
                message = f"""Skipped {sys.name}/{this_branch.name}: branch folder not found: {this_branch.path}
Check whether it was moved or renamed; restore the folder to resume.
To restart this branch, run in Python:
    from scope.read_write import load_binary
    system = load_binary({sys_path!r})
    system.remove_branch({this_branch.name!r})
    system.save()
Then relaunch the SCOPE tasks.
"""
                print(f"RUN_TASK: {message}", end='')
                report += message
                unavailable_branches.add(this_branch.name)
            continue

        ##########################
        ### STEP 1.5: WORKFLOW ###
        ##########################
        if debug > 0:   print(f"RUN_TASK, step 1.5: Evaluating list of Workflows in input: {job_data.workflow=}")
        job_data.workflow = interpret_workflow_from_user(sys, job_data, debug=debug)
  
        for wrk in job_data.workflow if isinstance(job_data.workflow, list) else list([job_data.workflow]):   ### Works when job_data.workflow is a str or a list

            exists, this_workflow             = this_branch.find_workflow(wrk)
            if not exists: this_workflow      = this_branch.add_workflow(wrk); sys.should_save = True
            if this_workflow is None: 
                raise ValueError(f"RUN_TASK, step 2.0: WORKFLOW could not be created. Make sure a Source with name {wrk} exist in System")
            if debug == 1:   print(f"RUN_TASK, step 2.0: Evaluating WORKFLOW with name={this_workflow.name}")
            elif debug > 1:  print(f"RUN_TASK, step 2.0: Evaluating WORKFLOW below: \n {this_workflow}")

            ###################
            ### STEP 2: JOB ###
            ###################
            ## 2.1 Finds or creates the job.
            print(f"RUN_TASK, step 2.1: Evaluating JOB with name={job_data.job}")
            exists, this_job = this_workflow.find_job(job_data=job_data, debug=debug)
            if not exists: this_job = this_workflow.add_job(job_data); sys.should_save = True

            print(f"RUN_TASK, step 2.2: Checking its job_data")
            ## 2.2 If job exists, it checks for input changes (also in qc_data for the computations)
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
                if debug > 0: print(f" {sys.name} -> {this_branch.name} -> {this_workflow.name} -> {this_job.name} -> Step: {comp.step} -> Run: {comp.run_number}")
                if debug > 0: print(f"--------------------------------------------------------------------------")
                if debug > 0: print(f"RUN_TASK, step 3.0: evaluating job, and computation with indices: {this_workflow.jobs.index(this_job)+1}/{len(this_workflow.jobs)}, {jdx+1}/{len(this_job.computations)}")

                ## Step 3.1-Checks files and updates
                qc_has_updated = comp.check_qc_data(inp_path=inp_path, debug=0)   ## Checks wether the user has updated the qc_data
                comp.check_updates()                                              ## Checks for not-registered update computations 
                comp.check_files()

                # Deleted outputs invalidate cached completion and parser data.
                if not comp.output_exists and (comp.isregistered or getattr(comp,'isfinished',False) or getattr(comp,'isgood',False) or hasattr(comp,'output') or hasattr(comp,'output_lines')):
                    comp.isregistered = False
                    comp.isfinished   = False
                    comp.isgood       = False
                    comp.delete_lines()
                    comp.delete_output()
                    sys.should_save = True

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

                        # A saved State can outlive deleted files from earlier steps.
                        if this_job.job_setup == 'rep_opt' and comp.step > 1:
                            previous_steps = {}
                            for previous in this_job.computations:
                                if previous.keyword == comp.keyword and previous.step < comp.step:
                                    latest = previous_steps.get(previous.step)
                                    if latest is None or previous.run_number > latest.run_number: previous_steps[previous.step] = previous
                            waiting_step = None if comp.step-1 in previous_steps else comp.step-1
                            for previous in previous_steps.values():
                                previous.check_files()
                                if not previous.input_exists or not previous.output_exists or not previous.isregistered or not getattr(previous,'isfinished',False) or previous.has_update:
                                    waiting_step = previous.step if waiting_step is None else min(waiting_step, previous.step)
                            if waiting_step is not None:
                                if debug > 0: print(f"RUN_TASK, step 3.2a: waiting for rep_opt step {waiting_step}, keyword={comp.keyword!r}, to finish and register")
                                continue

                        comp.check_submission_status(global_env, debug=debug)
                        if debug > 0: print("RUN_TASK, step 3.2a: checking submission status: isrunning=",comp.isrunning)
                        if debug > 0: print("RUN_TASK, step 3.2a: ignore_submitted=", options.ignore_submitted)
                        if debug > 0: print("RUN_TASK, step 3.2a: initial state is", comp.qc_data.istate)
                        if not comp.isrunning or options.ignore_submitted:       ## options.ignore_submitted will also be checked in comp.run 
                            # A stored completion flag cannot replace the required geometry.
                            if not comp.verify_state(comp.qc_data.istate):
                                if debug > 0: print(f"RUN_TASK, step 3.2a: initial state '{comp.qc_data.istate}' is missing or lacks geometry")
                                report += f"Retry pending for {comp.out_path}: initial state '{comp.qc_data.istate}' is missing or lacks geometry.\n"
                                continue
                            comp.run(global_env, options, debug=debug); sys.should_save = True
                    else: 
                        if debug > 0: print("RUN_TASK, step 3.2a: want_submit is False or comp.has_update")

                ## Step 3.2b-Warns if output exists but not input
                elif comp.output_exists and not comp.input_exists:  
                    if not comp.has_update:
                        report += f"Missing input for current attempt: {comp.out_path}. Inspect the computation files.\n"
                    if debug > 0: print(f"RUN_TASK: Missing input for {comp.out_path}")

                ############################
                ### STEP 4: REGISTRATION ###
                ############################
                elif comp.output_exists and comp.input_exists:
                    ## 4.1 Parse new outputs; cached registrations still reach continuation decisions in Step 5.
                    worked = True
                    if comp.isregistered:
                        if debug > 0: print(f"RUN_TASK, step 4: job already registered")
                    else:
                        if debug > 0: print(f"RUN_TASK, step 4: registration")
                        worked = comp.register(debug=debug)
                        sys.should_save = True

                        if debug > 0: print(f"RUN_TASK, step 4.1: registration {worked=}")
                        if debug > 0: print(f"RUN_TASK, step 4.1: {comp.has_update=}")
                        if debug > 0: print(f"RUN_TASK, step 4.1: {options.overwrite_inputs=}")

                    if this_job.job_setup == 'rep_opt' and getattr(comp,'energy',None) is not None:
                        # Preserve this attempt's energy even when a later step overwrote fstate.
                        if comp.step > len(this_job.energies): this_job.energies = np.pad(this_job.energies, (0,comp.step-len(this_job.energies)))
                        this_job.energies[comp.step-1] = comp.energy

                    #############################
                    ### STEP 5: CONTINUATIONS ###
                    #############################
                    ## 5.1 Generate or reuse continuations, including when the parent was registered earlier.
                    ## has_update can mean only that a file exists; the native setter must still create its object.
                    ncomputations = len(this_job.computations)
                    continuation_stopped = False
                    retry_running = False
                    if comp.status == 'no_scf_convergence':
                        if debug > 0: print(f"RUN_TASK, step 5.1: setting continuation computation with typ=scf")
                        new_comp = this_job.set_continuation_computation(comp, "scf", debug=debug)
                        if new_comp is None:
                            continuation_stopped = True
                            report += f"Automatic QE SCF recovery exhausted: {comp.out_path} : Review the electronic and SCF settings.\n"
                        elif new_comp.qc_data.software != "qe" and new_comp.run_number >= 10: report += f"Investigate {new_comp.out_path} : Run Number > 10 \n"
                    elif comp.status == "aborted":
                        # Aborted attempts are resubmitted under the same run, never continued.
                        if handle_errors and options.want_submit:
                            comp.check_submission_status(global_env, debug=debug)
                            retry_running = comp.isrunning and not options.ignore_submitted
                    elif not comp.isgood and this_job.must_be_good and this_job.job_setup != 'rep_opt':
                        if comp.type in ['opt','relax','vc-relax']:
                            if debug > 0: print(f"RUN_TASK, step 5.1: setting continuation computation with typ=opt")
                            new_comp = this_job.set_continuation_computation(comp, "opt", debug=debug)
                            if new_comp.run_number >= 10: report += f"Investigate {new_comp.out_path} : Run Number > 10 \n"
                        elif comp.type == 'ts':
                            if comp.qc_data.fctype != "calcall":
                                if debug > 0: print(f"RUN_TASK, step 5.1: setting continuation computation with typ=ts")
                                new_comp = this_job.set_continuation_computation(comp, "ts", debug=debug)
                            else:
                                if debug > 0: print(f"RUN_TASK, step 5.1: continuation computation for TS not set. calcAll already attempted")
                                report += f"Investigate {comp.out_path} : TS could not be found\n"

                    ## 5.2 Continue repetitive optimization from current attempts, retaining cached step energies.
                    if this_job.job_setup == 'rep_opt' and worked and comp.isfinished and comp.status in ['worked','exhausted','stucked'] and not comp.has_update:
                        converged = this_job.check_convergence(comp=comp, debug=debug)

                        ## Reuse an existing next step without calling the setter again for rep_opt.
                        exists, next_comp = this_job.find_computation(keyword=comp.keyword, step=comp.step+1, run_number=1)
                        if converged:
                            print(f"RUN_TASK, step 5.2: repetitive opt reached convergence")
                        elif comp.step < this_job.job_data.max_steps:
                            if not exists:
                                new_comp = this_job.set_continuation_computation(comp, "rep_opt", debug=debug)
                        else:
                            print(f"RUN_TASK, step 5.2: maximum steps reached without convergence")
                            report += f"Repetitive optimization exhausted {this_job.job_data.max_steps} steps without convergence: {comp.out_path}. Check energy, forces and QE optimization status.\n"

                    if len(this_job.computations) > ncomputations: sys.should_save = True

                    ## 5.3 Keep error handling after continuation creation so superseded outputs stay intact.
                    if debug > 0: print(f"RUN_TASK, step 5.3: registration {worked=}")
                    if not worked:
                        if comp.has_update:
                            # Keep failed outputs under their original attempt names.
                            comp.delete_lines()
                            comp.delete_output()
                            if debug > 0: print(f"RUN_TASK: Preserving superseded output: {comp.out_path}")
                        elif retry_running:
                            if debug > 0: print(f"RUN_TASK: Attempt is still queued or running; preserving {comp.out_path}")
                        elif continuation_stopped:
                            if debug > 0: print(f"RUN_TASK: Preserving {comp.out_path} for manual inspection")
                        elif handle_errors:          # ...takes default action
                            comp.read_lines()
                            empty_output = len(comp.output_lines) == 0
                            if not empty_output:
                                comp.store(debug=debug)  # Moves the failed log to _part.
                            else:
                                if os.path.isfile(comp.out_path): os.remove(comp.out_path)
                            comp.delete_lines()
                            comp.delete_output()
                            comp.check_files()
                            if comp.status == "aborted" or empty_output:
                                comp.isregistered = False
                                comp.isfinished   = False
                                comp.isgood       = False
                                if options.want_submit and comp.verify_state(comp.qc_data.istate):
                                    comp.run(global_env, options, debug=debug)
                                    if debug > 0: print(f"RUN_TASK: Same-run retry passed to submission logic: {comp.out_path}")
                                elif options.want_submit:
                                    report += f"Retry pending for {comp.out_path}: initial state '{comp.qc_data.istate}' is missing or lacks geometry.\n"
                                else:
                                    report += f"Retry pending for {comp.out_path}: automatic submission is disabled.\n"
                            else:
                                report += f"Retry pending for {comp.out_path} --> Run SCOPE again to prepare and submit the next attempt.\n"
                            if debug > 0: print(f"RUN_TASK: Errors handled for {comp.out_path}")
                        else:
                            if not comp.has_update: report += f"Registration failed for current attempt: {comp.out_path}. Inspect the output or rerun with automatic error handling (-e).\n"
                            if debug > 0: print(f"RUN_TASK: Error registering {comp.out_path}")
                    if debug > 0: print("")

            # Aggregate Computation status into the Job, Workflow and Branch.
            this_job.register(debug=0)
        this_workflow.register(debug=0)
        this_branch.register(debug=0)

        if getattr(sys,'should_save',False):
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
