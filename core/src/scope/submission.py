"""Builds scheduler directives and Bash commands for SCOPE submission files."""


def _wall_time_seconds(minutes: int | float | str | None) -> int | None:
    """Returns whole seconds for stored minutes or a legacy clock, rounding fractional seconds up."""
    from scope.classes_queue import _read_queue_time
    if minutes is None: return None
    if isinstance(minutes,str): return _read_queue_time(minutes,'slurm')
    if isinstance(minutes,int) or minutes.is_integer(): return int(minutes)*60
    seconds = int(minutes*60)
    # Compare in minutes so an exact whole second survives the conversion back from minutes.
    if seconds/60 < minutes: seconds += 1
    return seconds


def get_slurm_signal_time(time_limit: int | float | str | None, debug: int=0) -> int:
    """Returns advance-warning seconds for a SLURM wall-time limit.

    Uses 10% of a finite limit, between 10 and 60 seconds. Accepts minutes or
    a Slurm duration string. An unavailable, unlimited or zero limit uses 60 seconds.
    """
    seconds = _wall_time_seconds(time_limit)
    if seconds is None or seconds <= 0: return 60
    return min(60,max(10,seconds//10))


def gen_job_priority_directive(environment: object) -> str:
    """Returns a Slurm or SGE submission directive, or an empty string at normal priority."""
    # Defines the Penalty for low and lowest priority
    job_priority_penalties = {'normal': 0, 'low': 50, 'lowest': 100}

    penalty = job_priority_penalties[getattr(environment,'job_priority','normal')]
    if penalty == 0: return ''
    if environment.scheduler == 'slurm': return f'#SBATCH --nice={penalty}'
    if environment.scheduler == 'sge':   return f'#$ -p -{penalty}'
    return ''


def gen_slurm_time_directives(queue: object) -> str:
    """Returns Slurm's wall-time and warning directives. An unknown limit produces only the warning directive."""
    from scope.classes_queue import get_submission_time_limit
    time_limit = get_submission_time_limit(queue)
    if isinstance(time_limit,str):
        # Older Queues can store a clock with seconds. Slurm grants whole minutes.
        seconds = _wall_time_seconds(time_limit)
        time_limit = None if seconds is None else (seconds+59)//60
    directives = []
    if time_limit is not None:
        days, remainder = divmod(time_limit,24*60)
        hours, minutes  = divmod(remainder,60)
        day_prefix = f'{days}-' if days else ''
        directives.append(f'#SBATCH --time={day_prefix}{hours:02d}:{minutes:02d}:00')
    directives.append(f'#SBATCH --signal=B:USR1@{get_slurm_signal_time(time_limit)}')
    return '\n'.join(directives)


def gen_sge_time_directives(queue: object, debug: int=0) -> tuple[str,str]:
    """Returns SGE wall-time directives and the signal that triggers result copying.

    Uses the task's time_limit, the queue's finite maximum, or 30 days for an
    unlimited queue. Omits h_rt when the hard limit is unknown. Older Queues collect
    missing settings when possible; unavailable settings use -notify/USR2.

    Without a finite soft limit, -notify requests USR2 before termination; USR1
    then denotes suspension. A finite s_rt instead uses USR1 for termination
    warning and omits -notify to avoid that ambiguity. Its value is reduced when
    necessary to leave copying time before h_rt. The queue's notify interval
    controls the grace period; zero or unavailable intervals cannot guarantee it.
    Assumes SGE's standard warning signals.
    """
    from scope.classes_queue import get_submission_time_limit
    if hasattr(queue,'update_settings'): queue.update_settings(debug=max(debug-1,0))

    # 1) Converts the requested and configured wall times from minutes to seconds.
    time_limit = _wall_time_seconds(get_submission_time_limit(queue))
    soft_limit = _wall_time_seconds(getattr(queue,'soft_time_limit',None))
    notify     = getattr(queue,'notify_time',None)

    # 2) Selects the appropriate SGE warning without confusing it with suspension.
    directives = []
    if soft_limit is None:
        directives.append('#$ -notify')
        warning_signal = 'USR2'
    else:
        warning_signal = 'USR1'
        if time_limit is not None:
            margin = int(notify) if notify is not None and notify > 0 else min(60,max(1,time_limit//10))
            soft_limit = min(soft_limit,max(1,time_limit-margin))

    # 3) Writes hours without a day prefix, as required by SGE.
    for resource, seconds in [('h_rt',time_limit),('s_rt',soft_limit)]:
        if seconds is not None:
            hours, remainder = divmod(seconds,3600)
            minutes, seconds = divmod(remainder,60)
            directives.append(f'#$ -l {resource}={hours:02d}:{minutes:02d}:{seconds:02d}')
    if debug > 0: print(f'GEN_SGE_TIME_DIRECTIVES: {queue.name}: warning={warning_signal}, notify={notify} s')
    return '\n'.join(directives), warning_signal


def gen_run_block(run_command: str, warning_signal: str='USR1', notify_time: float | None=None, debug: int=0) -> str:
    """Returns Bash commands that stop a calculation before copying partial results.

    Runs the command and its local children in a separate process group.
    Allows at most two seconds after TERM before forcing that group to stop,
    using at most half of notify_time when that interval is known. Sets
    run_status=124 on interruption, otherwise 0 as in existing submissions.
    The caller copies results, cleans scratch and exits with run_status.
    warning_signal is USR1 for SLURM and SGE soft limits, or USR2 for SGE
    -notify. With USR2, suspension warnings on USR1 are ignored. TERM/INT
    also trigger copying. SIGKILL cannot be handled; scheduler grace periods
    must allow enough time to stop the calculation and copy its files.
    """
    if warning_signal not in ('USR1','USR2'): raise ValueError('GEN_RUN_BLOCK: warning_signal must be USR1 or USR2')
    signals       = f'{warning_signal} TERM INT'
    stop_timeout  = 2 if notify_time is None else min(2,max(0,notify_time/2))
    stop_checks   = int(stop_timeout*10)
    suspension    = "# Leave SGE suspension warnings to the scheduler\ntrap '' USR1\n\n" if warning_signal == 'USR2' else ''
    return suspension + f'''# Stop the calculation on a time warning or termination request
run_status=0
calculation_pid=
stop_calculation() {{
    trap '' {signals}
    run_status=124
    echo "SCOPE: Interrupted; returning partial output"
    if [ -n "$calculation_pid" ]; then
        kill -TERM -- "-$calculation_pid" 2>/dev/null

        # Allow up to {stop_checks/10:g} seconds for the calculation and local children to stop
        for ((attempt=0; attempt<{stop_checks}; attempt++)); do
            if ! kill -0 -- "-$calculation_pid" 2>/dev/null; then break; fi
            sleep 0.1
        done
        if kill -0 -- "-$calculation_pid" 2>/dev/null; then
            echo "SCOPE: Calculation did not stop; forcing termination"
            kill -KILL -- "-$calculation_pid" 2>/dev/null
        fi
        wait "$calculation_pid" 2>/dev/null
    fi
}}
trap stop_calculation {signals}

# Start the calculation in its own process group, then wait for it
set -m
{run_command} &
calculation_pid=$!
set +m
if [ "$run_status" -ne 0 ]; then stop_calculation; fi  # Warning during startup
wait "$calculation_pid" 2>/dev/null

# The calculation has stopped; allow the result-copy blocks below to finish
trap '' {signals}
'''
