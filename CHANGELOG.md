# Changelog

## [1.0.1][1.0.1]

Changes since `v1.0.0`, covering 43 commits through `823b1ed` (8 November 2026). This update primarily affects `scope-qc`; the add-ons have citation updates but no functional code changes in this range.

### Quantum ESPRESSO

- Added optional wavefunction storage and reuse through `job_data.save_wf`: `False` disables it, `'job'` shares one evolving wavefunction per Computation keyword chain within a Job, and `'workflow'` allows sharing between Jobs of a Workflow. `True` is an alias for `'job'`. The directory of the wavefunction is checked before reuse; incomplete new data do not replace a previously stored wavefunction. Persistent copies omit pseudopotentials identified in the saved XML. ([61cb4f7][61cb4f7], [017108e][017108e], [dadc18b][dadc18b])
- Replaced random SCF mixing changes with recovery based on the observed convergence trend. Continuations adjust mixing parameters and electronic iteration limits, can switch to `local-TF`, and stop after limited recovery attempts. ([3940355][3940355])
- Corrected energy and force registration for the final SCF that QE may perform using the optimized unit cell after `vc-relax`. Failed or incomplete final SCFs no longer silently supply earlier energy/force values. Optimization status distinguishes BFGS convergence from exhausted or stalled optimization. ([c89a409][c89a409])
- Exposed mixing and optimization controls, including `mix_mode`, `mix_ndim`, `opt_maxstep` (`nstep`), `bfgs_ndim`, `trust_rad_ini`, `trust_rad_min`, `upscale`, `ion_dynamics`, and `cell_dynamics`. Several defaults changed; see below. ([61cb4f7][61cb4f7], [f061d1b][f061d1b], [bd19a60][bd19a60])
- Generated inputs now include `nosym = .true.` in `&SYSTEM`. Cutoff selection enforces `ecutrho >= 10 * ecutwfc` after considering all species. ([ece74b5][ece74b5], [823b1ed][823b1ed])
- Corrected coordinate output to preserve the requested State's atom order while matching source spin assignments. Ambiguous reordered structures with different spins on identical elements are rejected. ([9ff7c95][9ff7c95])

### Workflow execution and HPC files

- QE repetitive optimization now requires energy stability, forces below `forc_conv`, and explicit optimizer convergence. Each Computation retains its own registered energy and forces, so later State updates cannot overwrite earlier convergence data. Existing next steps are reused, `max_steps` includes the first step, and unconverged Jobs block dependencies that require success. ([c89a409][c89a409], [132414b][132414b])
- Pending `rep_opt` steps wait for finished, registered outputs from the latest runs of earlier steps in their keyword chain. Deleting an output invalidates cached completion and parser data, preventing retained continuation steps from submitting together during a restart. ([2160573][2160573])
- With automatic error handling and submission enabled, aborted or empty-output attempts can be retried during the same invocation under the same run number. Failed non-empty outputs are archived when error handling applies; superseded outputs are preserved. Corrected a regression in this recovery path. ([0da8166][0da8166], [5089977][5089977], [9ace563][9ace563])
- Missing folders belonging to stored Branches produce an explanatory warning and are skipped, with instructions to restore the folder or remove the Branch before restarting. ([832224f][832224f])
- QE and Gaussian Slurm/SGE scripts use isolated work directories under `TMPDIR` and remove them after copying outputs and requested wavefunction/checkpoint files. Copy failures retain scratch and report its location. Scheduler stdout/stderr files remain in the computation folder, and stdout reports the execution directory. ([6e68efb][6e68efb], [ff46d9e][ff46d9e], [b49549f][b49549f])
- Added `System.should_save` to track changes requiring persistence. Successful saving clears the flag in memory and on disk; failed saving leaves it set. `run_task()` uses this shared flag instead of a separate local "update" flag. ([891f76b][891f76b], [132414b][132414b])
- The CLI summary focuses on unresolved issues requiring manual attention, avoiding obsolete resubmission messages for attempts already handled or superseded. ([be07e73][be07e73], [9ace563][9ace563])

### Gaussian

- Added `qc_data.print_forces` for SCF calculations and `qc_data.unrestricted` to request unrestricted calculations for singlet systems. Both default to `False`. ([78316ac][78316ac], [d5815eb][d5815eb])
- Corrected input-header assembly in submission scripts to prevent a race between `cat` and `mv` that could produce “Route card not found” errors. ([65e6d3c][65e6d3c])

### Normal-mode displacements and sampling

- Corrected mass weighting and unit conversion in `geom_sampling_from_vnm()` and `displace_coords_with_vnm()`. Cartesian geometries use Angstrom, normal coordinates use `sqrt(electron mass) * Bohr`, and harmonic potential increments use Hartree. Stored VNMs are not modified during these operations. ([adcccd2][adcccd2])
- `displace_coords_with_vnm()` now interprets `amplitude` as a fixed normal-coordinate displacement per selected mode, with default `10`, replacing the previous frequency-dependent displacement rule. Old amplitude values are not numerically equivalent. ([adcccd2][adcccd2])
- Removed the implicit cap on sampling widths and disabled frequency damping and the lower taper by default. Added the optional smooth upper taper `freq_top_limit` to both sampling functions. Docstrings explain how attenuation, connectivity rejection, additional sampling rounds and furthest-point sampling affect the distribution. ([adcccd2][adcccd2], [6f555b0][6f555b0])
- Corrected cumulative-coordinate handling and candidate selection across sampling rounds, including empty/rejected-candidate cases. Furthest-point sampling selects distinct indices and validates candidate counts and finite distances. Returned harmonic energies describe the final cumulative normal coordinates. ([adcccd2][adcccd2])

### File handling and documentation

- Moved `read_lines_file()` to `scope.read_write`. It removes NUL characters and returns a Python list instead of a fixed-width NumPy string array, avoiding excessive memory allocation for corrupted files containing exceptionally long lines. ([84d968a][84d968a])
- Added `overwrite=False` to `write_xyz()`. Existing files are preserved unless `overwrite=True` or `append=True` is supplied; combining both options raises an error. ([df3c634][df3c634])
- Fixed `Input_data` comparisons to include public names containing underscores and dictionary-valued settings, with exact NumPy array comparisons. ([c2aa7e5][c2aa7e5])
- Clarified that `main` is the development branch and that releases/PyPI provide stable versions. Updated manuscript citations to DOI `10.1039/D6DD00236F` and added the Zenodo concept DOI `10.5281/zenodo.23337435`. Added repetitive-optimization documentation. ([f9a73cc][f9a73cc], [930ec22][930ec22], [89ee646][89ee646])

### Changed defaults and compatibility

| Setting                                                     | 1.0.0            | 1.0.1                                    |
| ----------------------------------------------------------- | ---------------- | ---------------------------------------- |
| QE `qc_data.elec_conv` (`conv_thr`, Ry)                 | `1e-5`         | `1e-8`                                 |
| QE `qc_data.forc_conv` (`forc_conv_thr`, Ry/Bohr)       | `1e-5`         | `1e-4`                                 |
| QE structural steps per calculation (`nstep`)             | Fixed at `300` | `qc_data.opt_maxstep`, default `100` |
| QE `upscale`                                              | Fixed at `10`  | `qc_data.upscale`, default `100`     |
| `rep_opt` `job_data.max_steps`                          | `10`           | `5`                                    |
| Sampling `sigma_damp_factor`                              | `1`            | `0`                                    |
| `geom_sampling_from_vnm()` `freq_bottom_limit` (cm^-1)  | `80`           | `0`                                    |
| `State.sample_geometries()` `freq_bottom_limit` (cm^-1) | `50`           | `0`                                    |

- Review numerical results when changing versions: QE convergence thresholds, generated inputs and normal-mode displacement conventions have changed.
- Update external imports of `read_lines_file()` from `scope.parse_general` to `scope.read_write`, and adapt code that expects a NumPy array.
- Use `overwrite=True` for XYZ-writing calls that intentionally replace existing files.
- Launcher-check functions were renamed to `check_qe_launcher()` and `check_g16_launcher()`. ([61cb4f7][61cb4f7])
- Review positional calls to the sampling functions, which now accept `freq_top_limit`, and to `Job.check_convergence()`, which now accepts an optional Computation before `debug`. Prefer keyword arguments for these options.
- Regenerate submission files to obtain scratch cleanup and other script-generation changes.

[1.0.1]: https://github.com/QTC-IQAC/Scope/compare/v1.0.0...823b1ede5ade61ea4bf868bfcab09bccf6cbd625
[017108e]: https://github.com/QTC-IQAC/Scope/commit/017108e
[0da8166]: https://github.com/QTC-IQAC/Scope/commit/0da8166
[132414b]: https://github.com/QTC-IQAC/Scope/commit/132414b
[2160573]: https://github.com/QTC-IQAC/Scope/commit/2160573
[3940355]: https://github.com/QTC-IQAC/Scope/commit/3940355
[5089977]: https://github.com/QTC-IQAC/Scope/commit/5089977
[61cb4f7]: https://github.com/QTC-IQAC/Scope/commit/61cb4f7
[65e6d3c]: https://github.com/QTC-IQAC/Scope/commit/65e6d3c
[6e68efb]: https://github.com/QTC-IQAC/Scope/commit/6e68efb
[6f555b0]: https://github.com/QTC-IQAC/Scope/commit/6f555b0
[78316ac]: https://github.com/QTC-IQAC/Scope/commit/78316ac
[823b1ed]: https://github.com/QTC-IQAC/Scope/commit/823b1ed
[832224f]: https://github.com/QTC-IQAC/Scope/commit/832224f
[84d968a]: https://github.com/QTC-IQAC/Scope/commit/84d968a
[891f76b]: https://github.com/QTC-IQAC/Scope/commit/891f76b
[89ee646]: https://github.com/QTC-IQAC/Scope/commit/89ee646
[930ec22]: https://github.com/QTC-IQAC/Scope/commit/930ec22
[9ace563]: https://github.com/QTC-IQAC/Scope/commit/9ace563
[9ff7c95]: https://github.com/QTC-IQAC/Scope/commit/9ff7c95
[adcccd2]: https://github.com/QTC-IQAC/Scope/commit/adcccd2
[b49549f]: https://github.com/QTC-IQAC/Scope/commit/b49549f
[bd19a60]: https://github.com/QTC-IQAC/Scope/commit/bd19a60
[be07e73]: https://github.com/QTC-IQAC/Scope/commit/be07e73
[c2aa7e5]: https://github.com/QTC-IQAC/Scope/commit/c2aa7e5
[c89a409]: https://github.com/QTC-IQAC/Scope/commit/c89a409
[dadc18b]: https://github.com/QTC-IQAC/Scope/commit/dadc18b
[d5815eb]: https://github.com/QTC-IQAC/Scope/commit/d5815eb
[df3c634]: https://github.com/QTC-IQAC/Scope/commit/df3c634
[ece74b5]: https://github.com/QTC-IQAC/Scope/commit/ece74b5
[f061d1b]: https://github.com/QTC-IQAC/Scope/commit/f061d1b
[f9a73cc]: https://github.com/QTC-IQAC/Scope/commit/f9a73cc
[ff46d9e]: https://github.com/QTC-IQAC/Scope/commit/ff46d9e
