"""Alternative task runner with separate execution, recovery and continuation construction.

Use run_task_v2() explicitly; CLI execution is unchanged. Shared preflight checks
and native SCOPE registration/submission methods are reused. All continuation
selection and construction below is independent of Job's continuation setter.

Under Testing
"""
import os

import numpy as np
from scope.classes_input import interpret_workflow_from_user, set_input_data
from scope.read_write import load_binary
from scope.run_task import verify_files, verify_paths, verify_status

########
def run_task_v2(sys_path: str, inp_paths: list, global_env: str | object, handle_errors: bool=False, debug: int=0) -> str:
    """Processes tasks for a given System

    Parameters:
        sys_path (str): Path of a stored System file.
        inp_paths (list): SCOPE task files, processed in order.
        global_env (str | object): Environment file or loaded Environment.
        handle_errors (bool): Archives failed outputs and allow same-run retries.
        debug (int): Verbosity; called functions receive debug - 1.

    Returns:
        str: Unresolved failures or exhausted recovery strategies.

    New Computations might be created, and are processed during the same invocation.
    """
    # Verifies that files exist
    if not verify_files(sys_path, global_env, inp_paths, debug=max(debug-1,0)): raise ValueError('RUN_TASK_V2: File verification failed')
    # Checks for the presence of termination files (e.g. ISOLATED_FINISHED), which skip loading the System file 
    if not verify_status(sys_path, inp_paths, debug=max(debug-1,0)): return ''

    # Loads System
    system = load_binary(sys_path)

    # Loads Environment
    if isinstance(global_env,str): global_env = load_binary(global_env)
    if not verify_paths(global_env, system, debug=max(debug-1,0)): raise ValueError('RUN_TASK_V2: Path verification failed')
    report = ''


    unavailable_branches = set()

    for inp_path in inp_paths:
        _, options, job_data, qc_data = set_input_data(inp_path, debug=max(debug-1,0))
        
        # 0) Adds "local" task information to the global Environment 
        global_env.read_job_specs(inp_path, debug=max(debug-1,0))
        if global_env.scheduler == 'local':
            if debug > 0: print('RUN_TASK_V2: Local environment; submission disabled')
            options._mod_attr('want_submit',False)
            options._mod_attr('overwrite_inputs',False)
            options._mod_attr('overwrite_logs',False)

        # 1a) Resolves the requested branch
        exists, branch = system.find_branch(job_data.branch)
        if not exists: branch = system.add_branch(job_data.branch)
        if not os.path.isdir(branch.path):
            if branch.name not in unavailable_branches:
                message = (
                    f"Skipped {system.name}/{branch.name}: branch folder not found: {branch.path}\n"
                    "Check whether it was moved or renamed; restore the folder to resume.\n"
                    "To restart this branch, run in Python:\n"
                    "    from scope.read_write import load_binary\n"
                    f"    system = load_binary({sys_path!r})\n"
                    f"    system.remove_branch({branch.name!r})\n"
                    "    system.save()\n"
                    "Then relaunch the SCOPE tasks.\n"
                )
                print(message,end='')
                report += message
                unavailable_branches.add(branch.name)
            continue

        # 1b) Resolves the requested workflows
        for name in interpret_workflow_from_user(system, job_data, debug=max(debug-1,0)):
            exists, workflow = branch.find_workflow(name)
            if not exists:
                workflow = branch.add_workflow(name)
                system.should_save = True
            if workflow is None: 
                raise ValueError(f'RUN_TASK_V2: Could not create workflow for source {name!r}')

            # 1c) Resolves the requested job
            exists, job = workflow.find_job(job_data=job_data)
            if not exists:
                job = workflow.add_job(job_data)
                system.should_save = True
            else: job.check_job_data(inp_path=inp_path, debug=max(debug-1,0))
            if not job.check_requisites(debug=max(debug-1,0)): continue

            # 1d) Resolves the requested computations
            # Create missing computations without duplicating existing ones.
            # The loop also processes continuations created (and appended) during this invocation.
            report += job.set_computations_from_setup(qc_data, debug=max(debug-1,0)) or ''
            for comp in job.computations:
                if debug > 0: print(f'RUN_TASK_V2: {system.name}/{branch.name}/{workflow.name}/{job.name}: step {comp.step}, run {comp.run_number}, keyword={comp.keyword!r}')
                report += process_computation(comp, inp_path, global_env, options, handle_errors=handle_errors, debug=max(debug-1,0))

            # 2) Update Job, Workflow and Branch status from their children.
            job.register(debug=max(debug-1,0))
            workflow.register(debug=max(debug-1,0))
        branch.register(debug=max(debug-1,0))
        if getattr(system,'should_save',False): system.save()
    return report

########
def process_computation(comp: object, inp_path: str, environment: object, options: object, handle_errors: bool=False, debug: int=0) -> str:
    """Processes one Computation using the current task input and available files.

    Refreshes QC settings, checks for later runs, and locates input/output files.
    Missing output resets stored completion flags and triggers a submission
    check. With both files present, registers new results or reuses the existing
    registration, then evaluates continuations and failed-output handling.
    Output without its input is reported and left for inspection.

    Parameters:
        inp_path (str): Task file supplying the current QC settings.
        environment (object): Environment used for queue checks and submission.
        options (object): Submission and overwrite settings from the task.
        handle_errors (bool): Allows failed-output archiving and same-run retries.

    Returns:
        str: Messages requiring attention, or '' when there is nothing to report.
    """
    system = comp._job._workflow._branch._system
    comp.check_qc_data(inp_path=inp_path, debug=max(debug-1,0))
    comp.check_updates(debug=max(debug-1,0))
    comp.check_files()

    # Deleted outputs invalidate cached results without erasing stored States.
    if not comp.output_exists:
        if comp.isregistered or getattr(comp,'isfinished',False) or getattr(comp,'isgood',False) or hasattr(comp,'output') or hasattr(comp,'output_lines'):
            comp.isregistered = False
            comp.isfinished   = False
            comp.isgood       = False
            comp.delete_lines()
            comp.delete_output()
            system.should_save = True
        return handle_pending_computation(comp, environment, options, debug=max(debug-1,0))
    if not comp.input_exists:
        if comp.has_update: return ''
        return f'Missing input for current attempt: {comp.out_path} --> Inspect the computation files.\n'

    worked = True
    if not comp.isregistered:
        worked = comp.register(debug=max(debug-1,0))
        system.should_save = True
    return handle_continuation(comp, environment, options, registration_worked=worked, handle_errors=handle_errors, debug=max(debug-1,0))


def handle_pending_computation(comp: object, environment: object, options: object, debug: int=0) -> str:
    """Checks whether a Computation can be submitted.

    Calls comp.run() for resource selection and submission.
    Skips when submission is disabled or a later run exists. 
    Preceding steps with the same keyword must have their
    input/output files and be registered and finished. 
    A queued/running attempt waits unless options.ignore_submitted=True. 

    Returns:
        str: A report message, which can be '' if no issues are found. 
    """
    if not options.want_submit or comp.has_update: return ''

    # Later steps wait for preceding steps with the same keyword to finish and register.
    if comp.step > 1:
        previous_steps = {}
        for previous in comp._job.computations:
            if previous.keyword == comp.keyword and previous.step < comp.step:
                latest = previous_steps.get(previous.step)
                if latest is None or previous.run_number > latest.run_number: previous_steps[previous.step] = previous
        waiting_step = None if comp.step-1 in previous_steps else comp.step-1
        for previous in previous_steps.values():
            previous.check_files()
            if not previous.input_exists or not previous.output_exists or not previous.isregistered or not getattr(previous,'isfinished',False) or previous.has_update:
                waiting_step = previous.step if waiting_step is None else min(waiting_step,previous.step)
        if waiting_step is not None:
            if debug > 0: print(f'HANDLE_PENDING_COMPUTATION: Waiting for step {waiting_step}, keyword={comp.keyword!r}')
            return ''

    # In general, check whether the current attempt is queued or running before submission.
    comp.check_submission_status(environment, debug=max(debug-1,0))
    if comp.isrunning and not options.ignore_submitted: return ''

    # The Computation requires an initial State (istate) with geometry.
    if not comp.verify_state(comp.qc_data.istate): 
        return f"Retry pending for {comp.out_path} --> Initial state '{comp.qc_data.istate}' is missing or lacks geometry.\n"

    # Submit the current attempt, which creates input files and queues the job.
    comp.run(environment, options, debug=max(debug-1,0))
    # Marks the System for saving, even when submission fails, so that the new run number is recorded.
    comp._job._workflow._branch._system.should_save = True
    return ''


def handle_continuation(comp: object, environment: object, options: object, registration_worked: bool, handle_errors: bool=False, debug: int=0) -> str:
    """Evaluates recovery, optimization continuations and failed-output handling.

    Selects SCF recovery when comp.status is no_scf_convergence. Other recovery
    functions check whether they apply and supply their own reports.
    Creates or finds the recovery run before deciding whether to archive the
    failed output. If that run exists, keeps the previous run's output under
    its original filename. Aborted Computations create no new run; 
    their output may be archived before retrying the same run.

    Parameters:
        registration_worked (bool): Results were registered, or registration was
            already present. This does not mean comp.isgood=True.
        handle_errors (bool): Allows archiving and same-run retries after failed
            registration. SCF/optimization continuations can be created even
            when this is False.

    Returns:
        str: A report message, which can be '' if no issues are found. 
    """

    report = ''
    continuation_stopped = False

    # 1) SCF recovery
    if comp.status == 'no_scf_convergence':
        new_comp, message = recover_scf(comp, debug=max(debug-1,0))
        report += message
        continuation_stopped = new_comp is None
    # 2) Optimization recovery
    report += recover_optimization(comp, debug=max(debug-1,0))
    # 3) Repetitive optimization and minimum searches
    report += continue_rep_opt(comp, registration_worked, debug=max(debug-1,0))
    if registration_worked and comp._job.job_setup == 'find_minimum':
        report += comp._job.continue_minimum_search(comp, debug=max(debug-1,0))
    # 4) Failed-output handling
    report += handle_failed_output(comp, environment, options, handle_errors=handle_errors, continuation_stopped=continuation_stopped, registration_worked=registration_worked, debug=max(debug-1,0))
    return report


def recover_scf(comp: object, debug: int=0) -> tuple[object | None,str]:
    """Creates or reuses an SCF retry for the supplied Computation.

    Expects the caller to have identified an SCF failure. Delegates QC settings,
    attempt limits, construction and reporting to set_continuation_computation().

    Returns:
        tuple[Computation | None, str]: Continuation Computation and report message
        Continuation is None when no further automatic SCF retry is available.
    """
    new_comp, message = set_continuation_computation(comp, 'scf', debug=max(debug-1,0))
    return new_comp, message

def handle_failed_output(comp: object, environment: object, options: object, handle_errors: bool=False, continuation_stopped: bool=False, registration_worked: bool=False, debug: int=0) -> str:
    """Preserves, archives or retries an output whose registration failed.

    - This function is called as a last step when evaluation a potential continuation Computation. 
    - Preserves the file if a later run exists or recovery has stopped.
    - With error handling enabled, also preserves queued or running aborted
      attempts unless options.ignore_submitted=True.
    - Otherwise, it moves a nonempty output file to comp.out_path + '_part', or removes an
      empty output, and discards the stored lines and parser object.

    Aborted/empty attempts can be submitted again with the same step/run number
    when submission is enabled and the initial geometry exists. Other failures
    remain pending for a later invocation. No new Computation is created here.

    Parameters:
        handle_errors (bool): Permit file archiving/removal and automatic retry.
        continuation_stopped (bool): Recovery declined another attempt; retain
            the output for inspection.
        registration_worked (bool): Skips all file handling when True.

    Returns:
        str: A failure or pending-retry message, or '' if handled/preserved.
    """

    # Early exit when registration worked or a later run exists.
    if registration_worked: return ''

    # Early exit when a later run exists or recovery has stopped.
    if comp.has_update:
        comp.delete_lines()
        comp.delete_output()
        return ''
    if continuation_stopped: return ''

    # Deals with the current attempt's output when registration failed and no later run exists.
    if not handle_errors: return f'Registration failed for current attempt: {comp.out_path} --> Inspect the output or rerun with automatic error handling (-e).\n'

    # Preserves queued or running aborted attempts unless options.ignore_submitted=True.
    if comp.status == 'aborted' and options.want_submit:
        comp.check_submission_status(environment, debug=max(debug-1,0))
        if comp.isrunning and not options.ignore_submitted: return ''

    # Removes empty output or archives nonempty output to comp.out_path + '_part'.
    comp.read_lines()
    empty_output = len(comp.output_lines) == 0
    if empty_output:
        if os.path.isfile(comp.out_path): os.remove(comp.out_path)
    else: comp.store(debug=max(debug-1,0))
    # Deletes the stored lines and parser object, so that the next invocation will read the output file again.
    comp.delete_lines()
    comp.delete_output()
    comp.check_files()
    comp._job._workflow._branch._system.should_save = True

    # Attempts to submit the same run again (when submission is enabled) and the initial geometry exists.
    if comp.status == 'aborted' or empty_output:
        comp.isregistered = False
        comp.isfinished   = False
        comp.isgood       = False
        if options.want_submit and comp.verify_state(comp.qc_data.istate):
            comp.run(environment, options, debug=max(debug-1,0))
            return ''
        if options.want_submit: 
            return f"Retry pending for {comp.out_path} --> Initial state '{comp.qc_data.istate}' is missing or lacks geometry.\n"
        return f'Retry pending for {comp.out_path} --> Automatic submission is disabled.\n'
    return f'Retry pending for {comp.out_path} --> Run SCOPE again to prepare and submit the next attempt.\n'


def recover_optimization(comp: object, debug: int=0) -> str:
    """Selects recovery for an unsuccessful optimization or transition-state search.

    Acts only when job.must_be_good=True and the setup is not rep_opt.
    Leaves SCF failures to recover_scf() and aborted attempts to same-run retry.
    opt/relax/vc-relax reuse or create another run; 
    Creates/reuses Computations through set_continuation_computation(); files
    are left in place and no submission is made here.

    Returns:
        str: A run-count or failed-TS message, or '' when no report is needed.
    """

    # SCF failures are handled elsewhere
    if comp.status in ['no_scf_convergence','aborted']: 
        return ''

    # Rep_opt setups are handled elsewhere
    if comp.isgood or not comp._job.must_be_good or comp._job.job_setup == 'rep_opt': 
        return ''

    # Handles optimization-type jobs, including TS searches
    if comp.type in ['opt','relax','vc-relax']:
        _, message = set_continuation_computation(comp, 'opt', debug=max(debug-1,0))
        return message
    elif comp.type == 'ts':
        _, message = set_continuation_computation(comp, 'ts', debug=max(debug-1,0))
        return message
    return ''


def continue_rep_opt(comp: object, registration_worked: bool, debug: int=0) -> str:
    """Records a rep_opt attempt's energy and decides whether another step is needed.

    Does nothing for other Job setups beyond 'rep_opt'. 
    Stores available energy in job.energies in Hartree, even when this attempt cannot advance. 
    Only finished, registered attempts with no later run and status=worked/exhausted/stucked can advance.
    Job.check_convergence() checks the energy, forces and optimization status. 
    If not converged, creates the next step unless it exists or max_steps is reached. 
    Job.register() invoked later updates the overall Job status.

    Parameters:
        registration_worked (bool): Whether this attempt's registration worked
            or was already available; required before advancing.

    Returns:
        str: A step-limit message, or ''. Does not submit the next step.
    """

    job = comp._job

    # Early exit for other job_setups
    if job.job_setup != 'rep_opt': 
        return ''

    # Initializes energy convergence criterion, if it does not exist. 
    energy = getattr(comp,'energy',None)
    if energy is not None:
        if comp.step > len(job.energies): job.energies = np.pad(job.energies,(0,comp.step-len(job.energies)))
        job.energies[comp.step-1] = energy
    if not registration_worked or not comp.isfinished or comp.status not in ['worked','exhausted','stucked'] or comp.has_update: return ''
    if job.check_convergence(comp=comp, debug=max(debug-1,0)):
        if debug > 0: print(f'CONTINUE_REP_OPT: Converged at step {comp.step}, keyword={comp.keyword!r}')
        return ''
    exists, _ = job.find_computation(keyword=comp.keyword, step=comp.step+1, run_number=1)
    if comp.step < job.job_data.max_steps:
        if not exists:
            _, message = set_continuation_computation(comp, 'rep_opt', debug=max(debug-1,0))
            return message
        return ''
    return f'Repetitive optimization exhausted {job.job_data.max_steps} steps without convergence: {comp.out_path} --> Check energy, forces and optimization status.\n'


def set_continuation_computation(comp: object, typ: str, debug: int=0) -> tuple[object | None,str]:
    """Reuses or constructs the Computation following the supplied attempt.

    Parameters:
        comp (object): Preceding Computation, already attached to its Job.
        typ (str): scf/opt/ts keep the step and advance run_number; rep_opt
            advances the step and starts run 1. All retain the keyword.

    Selects a recovery function for the software and requested recovery type.
    That function receives job_data.recovery (default 'standard') and chooses
    QC settings, starting State, attempt limits and report messages. Unsupported
    software/type combinations raise ValueError.

    Reuses existing continuations without re-evaluating their strategy and
    returns an empty message. Otherwise constructs the selected continuation
    and returns its strategy's message. Repetitive optimizations use the
    software's optimization strategy for each new step.

    A new run marks the preceding run has_update=True; a rep_opt step does not.
    New Computations receive origin='<software>_<typ>_recovery_<strategy>'.
    Each strategy has its own attempt count; returning to a previous strategy
    resumes its count. Older origins without a strategy suffix do not count.
    Existing Computations retain their origin when reused or resubmitted.
    New Computations mark Job/Workflow status for registration and the System
    for saving. Input/output files are not generated or moved here, and no job
    is submitted. The Job continuation method is not called.

    Returns:
        tuple[Computation | None,str]: Existing/new successor and an optional
            report message. None indicates that the recovery strategy or
            rep_opt step limit permits no further attempt.
    """
    if typ not in ['scf','opt','ts','rep_opt']: raise ValueError(f'SET_CONTINUATION_COMPUTATION: Unknown type {typ!r}')
    job    = comp._job
    system = job._workflow._branch._system

    # 1) Selects the right function depending on software and typ variables; it selects the requested strategy.
    if comp.software == 'g16':
        if typ == 'scf':
            from scope.software.gaussian.g16_recovery import recover_g16_scf as get_recovery_settings
        elif typ in ['opt','rep_opt'] and comp.type == 'opt':
            from scope.software.gaussian.g16_recovery import recover_g16_opt as get_recovery_settings
        elif typ == 'ts' and comp.type == 'ts':
            from scope.software.gaussian.g16_recovery import recover_g16_ts as get_recovery_settings
        else:
            raise ValueError(f'SET_CONTINUATION_COMPUTATION: No {typ!r} recovery function for software={comp.software!r}, computation type={comp.type!r}')
    elif comp.software == 'qe':
        if typ == 'scf':
            from scope.software.quantum_espresso.qe_recovery import recover_qe_scf as get_recovery_settings
        elif typ in ['opt','rep_opt'] and comp.type in ['opt','relax','vc-relax']:
            from scope.software.quantum_espresso.qe_recovery import recover_qe_opt as get_recovery_settings
        else:
            raise ValueError(f'SET_CONTINUATION_COMPUTATION: No {typ!r} recovery function for software={comp.software!r}, computation type={comp.type!r}')
    else:
        raise ValueError(f'SET_CONTINUATION_COMPUTATION: No recovery functions exist for software={comp.software!r}')

    # 2) Reuses an existing run or step before evaluating new settings.
    if typ == 'rep_opt' and comp.step >= getattr(job.job_data,'max_steps',5): 
        return None, ''

    # 3) Sets Step and Run numbers
    step = comp.step+1 if typ == 'rep_opt' else comp.step
    run  = 1           if typ == 'rep_opt' else comp.run_number+1

    exists, new_comp = job.find_computation(keyword=comp.keyword, step=step, run_number=run)
    if exists:
        if typ != 'rep_opt' and not comp.has_update:
            comp.has_update = True
            system.should_save = True
        return new_comp, ''

    # 3) Obtains QC settings and a report from the requested strategy.
    recovery = getattr(job.job_data,'recovery','standard')
    qc_data, message = get_recovery_settings(comp, recovery=recovery, debug=max(debug-1,0))
    if qc_data is None:
        if comp.has_update:
            comp.has_update = False
            system.should_save = True
        return None, message

    # 4) Constructs the Computation and updates hierarchy registration flags.
    new_comp = job.add_computation(qc_data, step, comp.path, comp_keyword=comp.keyword, is_update=typ != 'rep_opt', debug=max(debug-1,0))
    new_comp.origin = f'{comp.software}_{typ}_recovery_{recovery}'
    new_comp.set_name()
    new_comp.set_paths()
    if typ != 'rep_opt': comp.has_update = True
    job.isregistered = False
    job._workflow.isregistered = False
    system.should_save = True
    if debug > 0: print(f'SET_CONTINUATION_COMPUTATION: {new_comp.name}, step {new_comp.step}, run {new_comp.run_number}')
    return new_comp, message
