# QE repetitive optimization

For `job_setup='rep_opt'`, `max_steps` is the maximum number of optimization steps in the SCOPE Job, including its first step. Each step is a separate QE calculation. QE's `nstep` (`qc_data.opt_maxstep` in SCOPE) independently limits the structural steps inside each calculation.

## Convergence

`Job.check_convergence(comp=None, debug=0)` checks a selected Computation or, without `comp`, the latest attempt of each keyword chain. A successful repetitive optimization requires:

- the current calculation has finished and registered successfully;
- QE explicitly reports BFGS convergence, rather than exhaustion or a stalled optimizer;
- consecutive steps of the same keyword have finite registered energies whose absolute difference meets `job_data.energy_thres` (Hartree, for the full calculated system);
- the current maximum absolute Cartesian force component is below `qc_data.forc_conv` (Ry/Bohr);

The first step cannot establish energy stability. Zero energy is a valid value. Missing, malformed, incomplete or nonfinite force data prevent convergence. The force metric is the largest absolute Cartesian component, matching QE BFGS, without division by atom count. SCOPE currently generates unconstrained atomic positions; this check is not a constraint-projected force criterion for externally edited constrained inputs.

`Computation.energy` (Hartree) and `Computation.forces` (QE: Ry/Bohr) retain the parsed results of each attempt. `Job.check_convergence()` checks the maximum absolute force component directly against `forc_conv`. Geometry States can be overwritten by subsequent steps, so repetitive convergence never recovers an earlier step's energy or force from the shared State. `Job.energies` remains available as a per-step summary. Previously registered Computations without the stored energy/force attributes need registration from their outputs before these criteria can be evaluated.

QE's explicit BFGS convergence also retains its native cell criterion for `vc-relax`. SCOPE does not independently re-evaluate the refreshed stress tensor in this change. Other backends retain their existing energy/native-optimization checks; a Gaussian force threshold is not inferred from QE's units.

## Energy, forces and geometry

After `vc-relax`, QE may calculate energy and forces again using a plane-wave basis built for the optimized unit cell. This final SCF has no new `ATOMIC_POSITIONS` card. `QE_output.get_energy_forces_block()` selects it when present, while geometry and cell readers retain the preceding block containing coordinates. If that SCF fails or is incomplete, its energy and forces are unavailable; older values are not substituted. Without the optional final SCF, the existing complete-block reader is used.

Generic registration stores available results and maintains State links. Existing State results are preserved when a calculation produces no replacement. The convergence check uses each Computation's own results, so an older value remaining on a shared State cannot establish convergence.

`parse_opt_status()` distinguishes `finished`, `exhausted`, `stucked` (the existing status spelling), and `not finished`. The generic `End of BFGS Geometry Optimization` marker alone does not prove convergence. The latest explicit outcome wins when multiple outcomes occur in the same output. Likewise, a later SCF failure overrides an earlier SCF success.

## Continuations and completion

Registered, terminated QE optimization attempts with usable geometry can create the next repetitive step even if the optimizer was exhausted or stuck. SCF failures retain the existing SCF recovery path. Aborted computations retain the existing same-run retry policy.

Existing next steps are reused. No step beyond `max_steps` is created. Exhaustion is reported for manual attention, and the Job remains not good; normal prerequisites requiring a good Job cannot advance. A finished execution is distinct from a converged optimization. Cached convergence/status changes are saved even when no output needed reparsing.

`Job.check_convergence(comp=comp)` checks the specified step against its predecessor; without `comp`, it checks the latest step of every keyword chain. It returns a decision without modifying Job status. `Job.register()` stores `isconverged` and `isgood`. Changed status marks the parent Workflow for aggregate registration and sets `System.should_save=True`; Workflow and Branch status updates use the same flag. These registration methods have no return value.

`run_task()` uses `System.should_save` as its only save flag. Computation processing, submissions, additions and tracked status changes mark it directly; no local `updated` flag is maintained. Rechecking cached, unchanged tasks performs no save. `System.save()` clears `should_save` in memory and in the serialized object after successful writing; failed writes leave the in-memory flag True. Source and Branch additions also mark the System for saving. This flag tracks explicit SCOPE operations; arbitrary manual attribute changes are not automatically detected.

Computation energy and force attributes are created during registration, rather than initialized for every new Computation.

References: [QE BFGS termination and criteria](https://github.com/QEF/q-e/blob/qe-7.0/Modules/bfgs_module.f90), [QE final refreshed SCF](https://github.com/QEF/q-e/blob/qe-7.0/PW/src/run_pwscf.f90).
