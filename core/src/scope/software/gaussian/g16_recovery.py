"""Defines selectable Gaussian recovery strategies.

Each function receives the preceding Computation and the strategy name from
job_data.recovery. It returns (qc_data, message): copied settings permit another
attempt; qc_data=None declines another continuation. The runner constructs and
submits Computations and handles output files separately.

get_attempt_number() supplies the next recovery attempt number for this
function's origin, counting created Computations in the same keyword/step up
to and including the current run. Other recovery types and later runs do not
count. Resubmission consumes no new attempt. Untagged older Computations do
not count. Origins include the strategy name, so each strategy has its own
attempt history. Switching back resumes that strategy's previous count.

A strategy edits a deepcopy of comp.qc_data through _mod_attr(), which keeps
attributes and the stored input dictionary consistent. istate selects the
starting State; fstate selects the State receiving the results. The current
Computation and its settings remain unchanged. Messages may be empty; report
messages end with a newline and separate paths from text with ' --> '.

A Job passes the same recovery name to every recovery function it needs.
A custom optimization strategy may therefore also need an SCF strategy with
that name, so it can handle electronic convergence failures during optimization.

To add a strategy, insert an elif branch before a function's final else and
update its list of available names. For example, a two-attempt TS strategy uses
the attempt_number calculated at the beginning of recover_g16_ts()::

    elif recovery == 'two_attempts':
        qc_data = deepcopy(comp.qc_data)
        if attempt_number == 1:
            qc_data._mod_attr('fctype','calcall')
        elif attempt_number == 2:
            qc_data._mod_attr('fctype','calcall')
            qc_data._mod_attr('tight_opt',True)
        else:
            return None, f'{comp.out_path} --> Two recovery attempts exhausted.\\n'
        return qc_data, ''

The input then selects recovery='two_attempts'. The example retains istate.
Strategies run only when a new continuation is needed;
existing continuations are reused without calling these functions.
"""
from copy import deepcopy
from scope.classes_workflow import get_attempt_number


def recover_g16_scf(comp: object, recovery: str='standard', debug: int=0) -> tuple[object | None,str]:
    """Chooses Gaussian SCF retry settings using the requested strategy.

    'standard' repeats the QC settings without a fixed attempt limit. It uses
    fstate as istate if that State has usable geometry, otherwise keeps istate.
    A warning from the tenth SCF recovery attempt permits further attempts.

    Parameters:
        comp (object): Current Computation in its Job's continuation history.
        recovery (str): Selects a branch below; currently 'standard'.
        debug (int): Verbosity level.

    Returns:
        tuple[Input_data | None,str]: Independent QC settings and report message.
            Changes no current settings and creates no Computation.
    """
    attempt_number = get_attempt_number(comp, f'g16_scf_recovery_{recovery}', debug=max(debug-1,0))

    if recovery == 'standard':
        # Retains all electronic-structure and SCF settings.
        # Only istate may change: uses fstate as the starting geometry when
        # that State has coordinates; otherwise keeps the current istate.

        # Repeated policy: every SCF recovery attempt uses the same settings.
        # 0) Copies the settings without modifying the preceding Computation.
        qc_data = deepcopy(comp.qc_data)

        # 1) Continues from the last geometry, when available.
        if hasattr(qc_data,'fstate') and comp.verify_state(qc_data.fstate,target='opt'):
            qc_data._mod_attr('istate',qc_data.fstate)

        # 2) Warns about repeated failures; this is not an attempt limit.
        next_run = comp.run_number+1
        message = ''
        if attempt_number >= 10:
            message = f'Repeated G16 SCF failure: {comp.out_path} --> Recovery continues with run {next_run} (SCF recovery attempt {attempt_number}); review the SCF settings.\n'
        return qc_data, message

    else:
        raise ValueError(f"RECOVER_G16_SCF: Unknown recovery strategy {recovery!r}; available: 'standard'")


def recover_g16_opt(comp: object, recovery: str='standard', debug: int=0) -> tuple[object | None,str]:
    """Chooses Gaussian optimization continuation settings.

    'standard' retains the QC settings and continues from the final State when
    it has usable geometry, otherwise from the initial State. Recovery runs
    have no fixed limit and warn from the tenth optimization recovery attempt.
    For rep_opt Jobs, the Job's convergence and max_steps determine whether a new step is
    needed; this function supplies its settings without a run-count warning.

    Parameters:
        comp (object): Preceding optimization Computation.
        recovery (str): Selects a branch below; currently 'standard'.
        debug (int): Verbosity level.

    Returns:
        tuple[Input_data | None,str]: Independent QC settings and report message.
            Changes no current settings and creates no Computation.
    """

    attempt_number = get_attempt_number(comp, f'g16_opt_recovery_{recovery}', debug=max(debug-1,0))

    if recovery == 'standard':
        # Retains the optimizer, convergence thresholds and electronic settings.
        # Only istate may change: uses fstate as the starting geometry when
        # that State has coordinates; otherwise keeps the current istate.

        # Repeated policy: every optimization recovery follows the same recipe.
        # 0) Copies the settings for another optimization.
        qc_data = deepcopy(comp.qc_data)

        # 1) Continues from the last geometry, when available.
        if hasattr(qc_data,'fstate') and comp.verify_state(qc_data.fstate,target='opt'):
            qc_data._mod_attr('istate',qc_data.fstate)

        # 2) Warns for recovery runs; rep_opt advances to a new step and run 1.
        next_run = comp.run_number+1
        is_rep_opt = comp._job.job_setup == 'rep_opt'
        message = ''
        if not is_rep_opt and attempt_number >= 10:
            message = f'Optimization remains unconverged: {comp.out_path} --> Recovery continues with run {next_run} (optimization recovery attempt {attempt_number}); review the optimization settings.\n'
        return qc_data, message

    else:
        raise ValueError(f"RECOVER_G16_OPT: Unknown recovery strategy {recovery!r}; available: 'standard'")


def recover_g16_ts(comp: object, recovery: str='standard', debug: int=0) -> tuple[object | None,str]:
    """Chooses Gaussian transition-state recovery settings.

    'standard' has one fixed recovery attempt: fctype='calcall', retaining
    istate. It stops after that attempt or if calcall was already selected.
    Interleaved SCF recoveries do not consume TS attempts.

    Parameters:
        comp (object): Preceding transition-state Computation.
        recovery (str): Selects a branch below; currently 'standard'.
        debug (int): Verbosity level.

    Returns:
        tuple[Input_data | None,str]: Independent QC settings and report message,
            or None and an explanation when no further retry is available.
    """
    attempt_number = get_attempt_number(comp, f'g16_ts_recovery_{recovery}', debug=max(debug-1,0))

    if recovery == 'standard':
        # Sets fctype='calcall', so Gaussian recalculates the force constants
        # at every geometry-optimization step. Keeps istate and all other
        # QC settings unchanged.

        # Fixed sequence: one TS recovery attempt, independently of SCF retries.
        qc_data = deepcopy(comp.qc_data)
        message = f'Automatic G16 TS recovery exhausted: {comp.out_path} --> No further automatic TS recovery is available (one calcall attempt).\n'

        if attempt_number == 1:
            # Recovery attempt 1: recomputes the force constants at every step.
            if qc_data.fctype == 'calcall': return None, message
            qc_data._mod_attr('fctype','calcall')
            if debug > 0: print('RECOVER_G16_TS: Attempt 1: using calcall and retaining the initial geometry')
            return qc_data, ''

        else:
            # No second TS recovery attempt is defined by this strategy.
            return None, message

    else:
        raise ValueError(f"RECOVER_G16_TS: Unknown recovery strategy {recovery!r}; available: 'standard'")
