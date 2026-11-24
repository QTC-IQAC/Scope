"""Defines selectable Quantum Espresso recovery strategies.

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
update its list of available names. For example, a two-attempt SCF strategy uses
the attempt_number calculated at the beginning of recover_qe_scf()::

    elif recovery == 'two_attempts':
        qc_data = deepcopy(comp.qc_data)
        if attempt_number == 1:
            qc_data._mod_attr('mix_beta',0.2)
        elif attempt_number == 2:
            qc_data._mod_attr('mix_beta',0.1)
            qc_data._mod_attr('mix_ndim',12)
        else:
            return None, f'{comp.out_path} --> Two recovery attempts exhausted.\\n'
        return qc_data, ''

The input then selects recovery='two_attempts'. The example retains istate.
Strategies run only when a new continuation is needed;
existing continuations are reused without calling these functions.
"""
from copy import deepcopy
from scope.classes_workflow import get_attempt_number


def recover_qe_scf(comp: object, recovery: str='standard', debug: int=0) -> tuple[object | None,str]:
    """Chooses QE SCF retry settings using the requested strategy.

    'standard' starts with an assessment of the final SCF cycle. 
    Oscillation status triggers reduced mixing_beta, increased mixing_ndim 
    and capped iterations; 
    Slow convergence increases the iteration limit. 
    Unknown behavior uses conservative mixing. It switches to local-TF at the 
    tenth SCF recovery attempt and onwards, or earlier when proposed
    settings repeat current ones. It stops after five local-TF SCF recoveries
    in the same keyword/step, or when settings already repeat in local-TF mode.
    It uses the final State when it has usable geometry, otherwise retains the initial State.

    'diverse' proposes three SCF recovery trials: plain mixing with beta=0.05
    and ndim=20; local-TF with beta=0.05 and ndim=20; then local-TF with
    beta=0.05 and ndim=4. It retains the preceding elec_maxstep allowance,
    with a minimum of 200. If the first recovery is slowly converging, it instead
    retains the mixing settings and extends the iteration allowance by the larger
    value between of 100 iterations or 50%, capped at 500; 
    an existing larger allowance is retained. Subsequent failures
    use trials 2 and 3. It stops after three SCF recovery Computations and uses
    the same final-State geometry selection as 'standard'.

    Parameters:
        comp (object): Preceding Computation. Its output supplies the final SCF
            behavior; its Job supplies the attempts with the same keyword/step.
        recovery (str): Selects 'standard' or 'diverse'.
        debug (int): Verbosity; the SCF assessment receives debug - 1.

    Returns:
        tuple[Input_data | None,str]: Independent QC settings and report message,
            or None and an explanation when recovery is exhausted. Creates the
            output parser when needed, but changes no current QC settings and
            creates no Computation.
    """
    origin = f'qe_scf_recovery_{recovery}'
    attempt_number = get_attempt_number(comp, origin, debug=max(debug-1,0))

    if recovery == 'standard':
        # Changes the mixing and SCF iteration settings as follows:
        # - Oscillating SCF: halves mix_beta, rounds it to two decimals and
        #   keeps it between 0.05 and 0.20. Adds 4 to mix_ndim, keeping the
        #   result between 12 and 20. Caps elec_maxstep at 100.
        # - Slowly converging SCF: increases elec_maxstep by the larger of
        #   100 iterations or 50%, with a ceiling of 500. Keeps mixing unchanged.
        # - Unclear SCF behavior: keeps mix_beta between 0.05 and 0.10 and
        #   raises mix_ndim to at least 12. Keeps elec_maxstep unchanged.
        # - Alternative mixing: sets mix_mode='local-TF', caps mix_beta at
        #   0.10 and raises mix_ndim to at least 12. Keeps elec_maxstep unchanged.
        # Uses fstate as istate when that State has coordinates; otherwise
        # keeps istate. All other QC settings remain unchanged.

        # 0) Defines the attempt limits for this strategy.
        switch_to_tf_on_attempt = 10
        max_local_tf_attempts   = 5
        exhausted = f'Automatic QE SCF recovery exhausted: {comp.out_path} --> Review the electronic and SCF settings before retrying.\n'

        # 1) Counts local-TF SCF recoveries through the current run in this chain.
        local_tf_attempts = 0
        for candidate in comp._job.computations:
            if getattr(candidate,'origin',None) != origin:
                continue
            if candidate.keyword != comp.keyword or candidate.step != comp.step:
                continue
            if candidate.run_number > comp.run_number:
                continue
            if getattr(candidate.qc_data,'mix_mode','plain').lower() == 'local-tf':
                local_tf_attempts += 1
        if local_tf_attempts >= max_local_tf_attempts:
            return None, exhausted

        # 2) Reads the convergence behavior and the current settings.
        if not hasattr(comp,'output'): 
            comp.create_output(debug=max(debug-1,0))

        behavior = comp.output.assess_scf_convergence(debug=max(debug-1,0))
        old_mix_beta = getattr(comp.qc_data,'mix_beta',0.7)
        old_mix_mode = getattr(comp.qc_data,'mix_mode','plain')
        old_mix_ndim = getattr(comp.qc_data,'mix_ndim',8)
        old_maxstep  = getattr(comp.qc_data,'elec_maxstep',200)
        new_mix_beta = old_mix_beta
        new_mix_mode = old_mix_mode
        new_mix_ndim = old_mix_ndim
        new_maxstep  = old_maxstep

        # 3) Chooses the settings for the next attempt.
        if attempt_number >= switch_to_tf_on_attempt and old_mix_mode.lower() != 'local-tf':
            # Uses local-TF from SCF recovery attempt 10, regardless of other retries.
            new_mix_beta = min(old_mix_beta,0.10)
            new_mix_mode = 'local-TF'
            new_mix_ndim = max(old_mix_ndim,12)

        elif behavior == 'oscillating':
            # Reduces mixing, increases its history, and caps iterations at 100.
            new_mix_beta = min(old_mix_beta*0.5,0.20)
            new_mix_beta = max(0.05,round(new_mix_beta,2))
            new_mix_ndim = max(old_mix_ndim+4,12)
            new_mix_ndim = min(new_mix_ndim,20)
            new_maxstep  = min(old_maxstep,100)

        elif behavior == 'slow':
            # Adds at least 100 iterations or 50%, with a ceiling of 500.
            new_maxstep = max(old_maxstep+100,int(old_maxstep*1.5))
            new_maxstep = min(new_maxstep,500)

        else:
            # Uses conservative mixing when the behavior cannot be classified.
            new_mix_beta = min(old_mix_beta,0.10)
            new_mix_beta = max(0.05,new_mix_beta)
            new_mix_ndim = max(old_mix_ndim,12)

        # 4) Avoids repeating unchanged settings: switches mode early or stops.
        old_settings = (old_mix_beta,old_mix_mode,old_mix_ndim,old_maxstep)
        new_settings = (new_mix_beta,new_mix_mode,new_mix_ndim,new_maxstep)
        if new_settings == old_settings:
            if old_mix_mode.lower() == 'local-tf':
                return None, exhausted
            new_mix_beta = min(old_mix_beta,0.10)
            new_mix_mode = 'local-TF'
            new_mix_ndim = max(old_mix_ndim,12)

        # 5) Copies the settings and applies only the selected changes.
        qc_data = deepcopy(comp.qc_data)
        qc_data._mod_attr('mix_beta',new_mix_beta)
        qc_data._mod_attr('mix_mode',new_mix_mode)
        qc_data._mod_attr('mix_ndim',new_mix_ndim)
        qc_data._mod_attr('elec_maxstep',new_maxstep)

        # 6) Continues from the last geometry, when available.
        if hasattr(qc_data,'fstate') and comp.verify_state(qc_data.fstate,target='opt'):
            qc_data._mod_attr('istate',qc_data.fstate)
        if debug > 0: print(f'RECOVER_QE_SCF: {behavior}; beta={new_mix_beta:g}, mode={new_mix_mode}, ndim={new_mix_ndim}, maxstep={new_maxstep}')
        return qc_data, ''

    elif recovery == 'diverse':
        # Trial 1: sets plain mixing, mix_beta=0.05 and mix_ndim=20.
        # Trial 2: sets local-TF mixing, mix_beta=0.05 and mix_ndim=20.
        # Trial 3: sets local-TF mixing, mix_beta=0.05 and mix_ndim=4.
        # The first trial keeps mixing unchanged for a slowly converging SCF
        # and instead extends elec_maxstep by at least 100 iterations or 50%,
        # capped at 500. An existing allowance above 500 is retained.
        # Every trial keeps the preceding elec_maxstep allowance at least 200.
        # Uses the last usable geometry; all other QC settings remain unchanged.
        # These are proposed trials; improved convergence is not established.

        # 0) Stops after three SCF recovery attempts.
        if attempt_number > 3:
            message = f"QE SCF recovery strategy 'diverse' exhausted: {comp.out_path} --> No further automatic SCF recovery after three attempts.\n"
            return None, message

        qc_data = deepcopy(comp.qc_data)
        old_maxstep = getattr(comp.qc_data,'elec_maxstep',200)
        new_maxstep = max(old_maxstep,200)

        # 1) Applies the selected recipe explicitly.
        if attempt_number == 1:
            if not hasattr(comp,'output'):
                comp.create_output(debug=max(debug-1,0))
            behavior = comp.output.assess_scf_convergence(debug=max(debug-1,0))

            if behavior == 'slow':
                # Keeps all mixing settings and extends only the iteration limit.
                extended_maxstep = max(old_maxstep+100,int(old_maxstep*1.5))
                extended_maxstep = min(extended_maxstep,500)
                new_maxstep = max(new_maxstep,extended_maxstep)
            else:
                qc_data._mod_attr('mix_mode','plain')
                qc_data._mod_attr('mix_beta',0.05)
                qc_data._mod_attr('mix_ndim',20)

        elif attempt_number == 2:
            qc_data._mod_attr('mix_mode','local-TF')
            qc_data._mod_attr('mix_beta',0.05)
            qc_data._mod_attr('mix_ndim',20)

        elif attempt_number == 3:
            qc_data._mod_attr('mix_mode','local-TF')
            qc_data._mod_attr('mix_beta',0.05)
            qc_data._mod_attr('mix_ndim',4)

        # 2) Preserves the iteration allowance and the final usable geometry.
        qc_data._mod_attr('elec_maxstep',new_maxstep)
        if hasattr(qc_data,'fstate') and comp.verify_state(qc_data.fstate,target='opt'):
            qc_data._mod_attr('istate',qc_data.fstate)
        if debug > 0: print(f'RECOVER_QE_SCF: diverse attempt {attempt_number}; beta={qc_data.mix_beta:g}, mode={qc_data.mix_mode}, ndim={qc_data.mix_ndim}, maxstep={qc_data.elec_maxstep}')
        return qc_data, ''

    else:
        raise ValueError(f"RECOVER_QE_SCF: Unknown recovery strategy {recovery!r}; available: 'standard', 'diverse'")


def recover_qe_opt(comp: object, recovery: str='standard', debug: int=0) -> tuple[object | None,str]:
    """Chooses QE optimization settings for relax and vc-relax continuations.

    'standard' retains the QC settings and continues from the final State when
    it has usable geometry, otherwise from the initial State. Recovery runs
    have no fixed limit and warn from the tenth optimization recovery attempt.
    For rep_opt Jobs, the Job's convergence and max_steps determine whether a new step is
    needed; this function supplies its settings without a run-count warning.

    'diverse' uses this same optimization policy; its different recipes apply
    only to SCF recovery.

    Parameters:
        comp (object): Preceding relax or vc-relax Computation.
        recovery (str): Selects 'standard' or 'diverse'; both use the same optimization policy.
        debug (int): Verbosity level.

    Returns:
        tuple[Input_data | None,str]: Independent QC settings and report message.
            Changes no current settings and creates no Computation.
    """
    attempt_number = get_attempt_number(comp, f'qe_opt_recovery_{recovery}', debug=max(debug-1,0))

    if recovery in ['standard','diverse']:
        # Retains the geometry/cell optimizer, convergence thresholds, cutoffs
        # and SCF settings. Only istate may change: uses fstate as the starting
        # geometry when that State has coordinates; otherwise keeps istate.

        # Repeated policy: every optimization recovery follows the same recipe until max_attempts.
        max_attempts = 10

        # 0) Copies the settings for another optimization.
        qc_data = deepcopy(comp.qc_data)

        # 1) Continues from the last geometry, when available.
        if hasattr(qc_data,'fstate') and comp.verify_state(qc_data.fstate,target='opt'):
            qc_data._mod_attr('istate',qc_data.fstate)

        # 2) Warns for recovery runs; rep_opt advances to a new step and run 1.
        next_run = comp.run_number+1
        is_rep_opt = comp._job.job_setup == 'rep_opt'
        message = ''
        if not is_rep_opt and attempt_number >= max_attempts:
            message = f'Optimization remains unconverged: {comp.out_path} --> Recovery continues with run {next_run} (optimization recovery attempt {attempt_number}); review the optimization settings.\n'
        return qc_data, message

    else:
        raise ValueError(f"RECOVER_QE_OPT: Unknown recovery strategy {recovery!r}; available: 'standard', 'diverse'")
