#!/usr/bin/env python3
import os
from scope.classes_state       import find_state
from scope.classes_environment import set_user
from scope.operations.dicts_and_lists import where_in_array
from scope.other               import get_metal_idxs, get_metal_species
from scope.elementdata         import ElementData 
elemdatabase = ElementData()
QE_PP_RELEASES_URL = "https://github.com/QTC-IQAC/Scope/releases/tag/PP_Libraries"

#######################
def gen_qe_input(comp: object, debug: int=0):
    """
    Generate a Quantum Espresso input file for a computation.

    Parameters:
        comp (object):    Computation object containing the QC settings, source, and file paths.
        debug (int):      Verbosity level.

    Returns:
        None

    The input starts from an atomic guess. Continuation batch scripts may switch
    their scratch copy to file initialization after preparing a valid saved wf.
    Optimization inputs include qc_data.bfgs_ndim in &ions (default 1, standard BFGS).

    Warning:
        Atomic coordinates and labels are taken from the requested state
        (`comp.qc_data.istate`) in State order. Spins come from comp.source. Reordered
        States require uniform spins per element, or an explicit atom mapping when
        chemically identical sites have different intended spins.
    """ 

    ######################
    ### INITIAL CHECKS ###
    ######################
    ## 1-Verify Information in initial state
    assert hasattr(comp.qc_data,"istate"), f"istate = {comp.qc_data.istate} not found in comp.qc_data"
    exists, istate = find_state(comp.source, comp.qc_data.istate)
    assert exists,                         f"istate = {comp.qc_data.istate} does not exist"
    assert hasattr(istate,"labels"),       f"istate = {comp.qc_data.istate} doesn't have labels"
    assert hasattr(istate,"coord"),        f"istate = {comp.qc_data.istate} doesn't have coordinates"

    ## 2-Cell or Molecule?
    if   comp.source.object_type == "cell":   system_type = "cell"
    elif comp.source.object_type == "specie": system_type = "molecule"
    else: print(f"GEN_QE_input: unrecognised {comp.source.object_type=}")

    ## 3-Determines the PP_Library path
    if not hasattr(comp.qc_data,"pp_library"): f"PP_Library could not be found. Please set it in the 'qc_data' section of the input. Options are: vanderbilt, efficiency and precision"
    import scope ## Trick to retrieve the relative path
    PP_root = os.path.abspath(scope.__file__).replace("__init__.py","software/quantum_espresso/PP_Libraries/")
    if   comp.qc_data.pp_library.lower() == "efficiency":  PP_path = PP_root + "Efficiency/"
    elif comp.qc_data.pp_library.lower() == "precision":   PP_path = PP_root + "Precision/"
    elif comp.qc_data.pp_library.lower() == "vanderbilt":  PP_path = PP_root + "Vanderbilt_USPP/"
    else: raise ValueError(f"GEN_QE_INPUT: could not recognise PP_library: {comp.qc_data.pp_library}")

    if comp.qc_data.pp_library.lower() in ["efficiency", "precision"]:
        if not os.path.isdir(PP_path) or not os.path.isfile(os.path.join(PP_path, "config.json")):
            lib_name = "Efficiency" if comp.qc_data.pp_library.lower() == "efficiency" else "Precision"
            raise FileNotFoundError(
                f"GEN_QE_INPUT: PP library '{lib_name}' could not be found in:\n"
                f"  {PP_path}\n"
                f"The scope-qc package only ships the Vanderbilt_USPP library.\n"
                f"Please download the '{lib_name}' pseudopotential archive from:\n"
                f"  {QE_PP_RELEASES_URL}\n"
                f"and extract the '{lib_name}' folder inside:\n"
                f"  {PP_root}\n"
                f"so that the file\n"
                f"  {os.path.join(PP_path, 'config.json')}\n"
                f"is present."
            )

    ## Checks requisites
    if system_type == "cell":     
        assert hasattr(istate,"cell_vector"),  f"GEN_QE_input: {comp.qc_data.istate=} doesn't have cell vectors"
        assert istate.charge == int(0),        f"GEN_QE_input: {istate.charge=} must be 0. The STATE charge must be 0"

    if debug >= 1 : print("GEN_QE_INPUT: system_type:", system_type)
    if debug >= 1 : print("GEN_QE_INPUT: creating input in path:", comp.inp_path)

    #########################
    ### DETERMINE SPECIES ###
    #########################
    metal_indices            = get_metal_idxs(istate.labels)           ## Indices of metal atoms in the state
    metal_species            = get_metal_species(istate.labels)        ## 
    species                  = get_label_spin_pairs(comp.source)       ## Spin information is taken from the source
    nspecies                 = len(species)

    # Coordinates follow State order; intended spins must describe those same atoms.
    source_labels = [atom.label for atom in comp.source.atoms]
    spins         = list(comp.source.atomic_spins)
    if len(istate.labels) != len(istate.coord) or sorted(istate.labels) != sorted(source_labels):
        raise ValueError("GEN_QE_INPUT: State and source must describe the same atoms")
    if list(istate.labels) != source_labels:
        spins_by_label = {}
        for label, spin in zip(source_labels, spins): spins_by_label.setdefault(label, set()).add(spin)
        if any(len(values) != 1 for values in spins_by_label.values()):
            raise ValueError("GEN_QE_INPUT: Reordered State with mixed site spins requires an explicit atom mapping")
        spins = [next(iter(spins_by_label[label])) for label in istate.labels]

    if debug >= 0: 
        print("GEN_QE_INPUT: Received cell with this data:")
        print(f"    species:  {species}")
        print(f"    metal_species: {metal_species}")
        print(f"    metal_indices: {metal_indices}")
        print(f"    ismagnetic:    {istate.ismagnetic}")

    ####################
    ### WRITES INPUT ###
    ####################
    with open(comp.inp_path, 'w+') as inp:
        print(" &control", file=inp)
        print(f"    calculation='{comp.qc_data.comp_type}'", file=inp)
        print( "    restart_mode='from_scratch'", file=inp)
        print(f"    pseudo_dir ='{PP_path}'", file=inp)
        print( "    disk_io='low'", file=inp)
        prefix = os.path.splitext(comp.wf_name)[0] if getattr(comp._job.job_data, 'save_wf', False) else comp.name
        print(f"    prefix='{prefix}'", file=inp)

        ## Keep QE scratch files in the working directory chosen by the batch script.
        print("    outdir='./'", file=inp)

        ## Print forces for finite differences computation.
        if hasattr(comp.qc_data,"print_forces"):
            if comp.qc_data.print_forces and comp.qc_data.comp_type == "scf": print(f"    tprnfor = .true.", file=inp)

        ## Keywords for Optimization-related computations 
        if comp.qc_data.comp_type == "opt" or comp.qc_data.comp_type == "relax" or comp.qc_data.comp_type == "vc-relax":
            print(f"    wf_collect = .true.", file=inp)
            print(f"    tprnfor = .true.", file=inp)
            print(f"    nstep = {comp.qc_data.opt_maxstep}", file=inp)
            print(f"    forc_conv_thr = {comp.qc_data.forc_conv}", file=inp)
        print("/", file=inp)

        #/////////////////////
        #// System control ///
        #/////////////////////
        min_cowfc = 25
        min_corho = 250
        ### Determines cutoffs using the info contained in the Pseudopotentials.
        ### Minimum values are also set
        for idx, spec in enumerate(species):
            pp, cutoff_wfc, cutoff_rho = get_pp(spec[0], PP_path)
            if cutoff_wfc > min_cowfc:      min_cowfc = cutoff_wfc      # updates
            if cutoff_rho > min_corho:      min_corho = cutoff_rho      # updates
            if cutoff_rho < 10*cutoff_wfc:  min_corho = 10*cutoff_wfc   # updates

        if comp.qc_data.comp_type == "vc-relax": min_cowfc *= 2; min_corho *= 2  ## In vc-relax it is convenient to minimize pulay stress

        print(" &system", file=inp)

        ## For isolated molecules
        if system_type == "molecule": 
            if not hasattr(comp.qc_data,"cubeside"): 
                raise ValueError(f"GEN_QE_INPUT: For molecules, the 'cubeside' parameter must be set in the qc_data section of the input. It is missing.") 
            print(f"    ibrav=1, celldm(1)={comp.qc_data.cubeside}", file=inp)

        ## For Crystals
        elif system_type == "cell":     
            print(f"    ibrav=0,", file=inp)

        print(f"    nat={istate.natoms}, ntyp={nspecies}, ecutwfc={int(min_cowfc)}, ecutrho={float(min_corho)}", file=inp)
        print(f"    nspin=2,", file=inp)
        print(f"    tot_charge={istate.charge}", file=inp)
        print(f"    tot_magnetization={istate.spin}", file=inp)
  
        ## Starting Magnetization
        for idx, spec in enumerate(species):
            if spec[1] != 0: print(f"    starting_magnetization({idx+1})={spec[1]}", file=inp)

        ## Hubbard: Here it is assumed that the Hubbard U term will apply to all metals
        if comp.qc_data.is_hubbard and comp.qc_data.version <= 7.0:
            print("    lda_plus_u=.true.,", file=inp)
            for idx, spec in enumerate(species):
                if spec[0] in metal_species: print(f"    Hubbard_U({idx+1})={comp.qc_data.uterm}", file=inp)

        if comp.qc_data.is_grimme:
            if comp.qc_data.grimme_type == "d2":
                print("    vdw_corr='grimme-d2'", file=inp)
                print("    dftd3_version=2", file=inp)
            elif comp.qc_data.grimme_type == "d3":
                print("    vdw_corr='grimme-d3'", file=inp)
                print("    dftd3_version=3", file=inp)
            elif comp.qc_data.grimme_type == "d3bj":
                print("    vdw_corr='grimme-d3'", file=inp)
                print("    dftd3_version=4", file=inp)
            elif comp.qc_data.grimme_type == "d3m":
                print("    vdw_corr='grimme-d3'", file=inp)
                print("    dftd3_version=5", file=inp)
            elif comp.qc_data.grimme_type == "d3mbj":
                print("    vdw_corr='grimme-d3'", file=inp)
                print("    dftd3_version=6", file=inp)
            
        if system_type == "molecule": print("    assume_isolated='mp'", file=inp)
        print("/", file=inp)
        
        #////////////////////////
        #// Electrons control ///
        #////////////////////////
        print(f" &electrons", file=inp)
        print(f"    diagonalization='david'", file=inp)
        print(f"    electron_maxstep= {comp.qc_data.elec_maxstep}", file=inp)
        print(f"    conv_thr={comp.qc_data.elec_conv}", file=inp)
        print(f"    mixing_beta={comp.qc_data.mix_beta}", file=inp)
        print(f"    mixing_mode='{comp.qc_data.mix_mode}'", file=inp)
        print(f"    mixing_ndim={comp.qc_data.mix_ndim}", file=inp)
        print(f"    startingpot='atomic'", file=inp)
        print(f"    startingwfc='atomic+random'", file=inp)
        print("/", file=inp)
   
        #///////////////////
        #// Ions control ///
        #///////////////////
        if comp.qc_data.comp_type == "opt" or comp.qc_data.comp_type == "relax" or comp.qc_data.comp_type == "vc-relax":
            print(" &ions", file=inp)
            print(f"    upscale={comp.qc_data.upscale}", file=inp)
            print(f"    ion_dynamics='{comp.qc_data.ion_dynamics}'", file=inp)
            print(f"    bfgs_ndim={comp.qc_data.bfgs_ndim}", file=inp)
            print(f"    trust_radius_ini= {comp.qc_data.trust_rad_ini}", file=inp)
            print(f"    trust_radius_min= {comp.qc_data.trust_rad_min}", file=inp)
            print("/", file=inp)
        
        #///////////////////
        #// Cell control ///
        #///////////////////
        if comp.qc_data.comp_type == "vc-relax":
            print(" &cell", file=inp)
            print(f"    cell_dynamics='{comp.qc_data.cell_dynamics}'", file=inp)
            print(f"    cell_dofree='all'", file=inp)
            print(f"    cell_factor= 1.2D0", file=inp)
            print(f"    press= {comp.qc_data.pressure}", file=inp)
            if comp.qc_data.pressure == 0: print("    press_conv_thr= 0.5D0", file=inp)
            else:                          print("    press_conv_thr= 0.01D0", file=inp)
            print("/", file=inp)

        #//////////////////
        #// Cell Params ///
        #//////////////////
        if system_type == "cell":
            print("CELL_PARAMETERS angstrom", file=inp)
            print(f"{istate.cell_vector[0][0]:14.6f} {istate.cell_vector[0][1]:14.6f} {istate.cell_vector[0][2]:14.6f}", file=inp)
            print(f"{istate.cell_vector[1][0]:14.6f} {istate.cell_vector[1][1]:14.6f} {istate.cell_vector[1][2]:14.6f}", file=inp)
            print(f"{istate.cell_vector[2][0]:14.6f} {istate.cell_vector[2][1]:14.6f} {istate.cell_vector[2][2]:14.6f}", file=inp)
         
        #///////////////////
        #// Atom Species ///
        #///////////////////
        print("ATOMIC_SPECIES", file=inp)
        for idx, spec in enumerate(species):
            if spec[1] != 0: label = f"{spec[0]}{spec[1]}"
            else:            label = spec[0]
            weight = elemdatabase.elementweight[spec[0]]
            pp, cutoff_wfc, cutoff_rho = get_pp(spec[0], PP_path)
            print(f"{label} {weight:6.4f} {pp}", file=inp)
        
        #/////////////////////////////////////////////////////////////////////
        #// Atom Coords, taken from istate, but spins taken from the source //
        #/////////////////////////////////////////////////////////////////////
        print("ATOMIC_POSITIONS angstrom", file=inp)
        for label, coord, spin in zip(istate.labels, istate.coord, spins):
            if spin != 0: label = f"{label}{spin}"
            print(f"{label:4}        {coord[0]:12.6f}   {coord[1]:12.6f}   {coord[2]:12.6f}", file=inp)
        print("K_POINTS gamma", file=inp)

        #//////////////////////////////////////////
        #// HUBBARD data for versions above 7.0 ///
        #//////////////////////////////////////////
        if comp.qc_data.is_hubbard and comp.qc_data.version > 7.0:
            print("HUBBARD (atomic)", file=inp)
            for idx, spec in enumerate(species):
                if spec[1] != 0: print(f"U {spec[0]}-3d {comp.qc_data.uterm}", file=inp)

###################################################
def gen_qe_subfile(comp: object, queue: object, module: str, procs: int=1, exe: str="pw.x", version: float=7.0, launcher: str=None, debug: int=0):
    """Generate a QE script with isolated scratch and optional evolving wf storage.

    save_wf='job' (also True) reuses data only for continuations. 'workflow'
    allows reuse from the first computation of every participating Job. The
    shared path is checked at runtime; unavailable/incomplete data leaves the
    atomic guess unchanged. Job requisites must enforce execution order.
    Scheduler stdout/stderr use absolute paths in the computation folder.
    """
    import shlex
    wf_prepare, wf_store = gen_qe_wf_blocks(comp, debug=debug)
    stderr_path = shlex.quote(os.path.abspath(os.path.join(comp.path, comp.name+'.stderr')))
    stdout_path = shlex.quote(os.path.abspath(os.path.join(comp.path, comp.name+'.stdout')))

    # Retrieve the preferred launcher from the queue's environment if not provided.
    if launcher is None: launcher = 'srun' if queue._environment.scheduler == 'slurm' else 'mpirun'
    launcher_args = shlex.split(launcher)
    if not launcher_args: raise ValueError('GEN_QE_SUBFILE: Launcher cannot be empty')
    launcher_name = os.path.basename(launcher_args[0])
    if launcher_name == 'srun': launcher_args += ['-n', str(procs)]
    elif launcher_name in ['mpirun', 'mpiexec']: launcher_args += ['-np', str(procs)]
    launch_command = f'{shlex.join(launcher_args)} {shlex.quote(exe)}'

    if debug > 0: print(f"GEN_QE_SUBFILE: Creating submission file for {comp.name}")
    if debug > 0: print(f"GEN_QE_SUBFILE: Using QE Module From Environment: '{module}'") 
    if debug > 0: print(f"GEN_QE_SUBFILE: Requesting Cores={procs}")
    version = str(version)

    with open(comp.sub_path, 'w+') as sub:
        if queue._environment.scheduler == 'slurm':
            print(f"#!/bin/bash", file=sub)
            print(f"#SBATCH -J {comp.name}", file=sub)
            print(f"#SBATCH -e {stderr_path}", file=sub)
            print(f"#SBATCH -o {stdout_path}", file=sub)
            print(f"#SBATCH -p {queue.name}", file=sub)
            print(f"#SBATCH --nodes=1", file=sub)
            print(f"#SBATCH --ntasks={procs}", file=sub)
            print(f"#SBATCH --mem-per-cpu=1900MB", file=sub)
            print(f"#SBATCH --time={queue.time_limit}", file=sub)

        elif queue._environment.scheduler == 'sge':
            print(f"#!/bin/bash", file=sub)
            print(f"#$ -N {comp.name}", file=sub)
            print(f"#$ -o {stdout_path}", file=sub)
            print(f"#$ -e {stderr_path}", file=sub)
            print(f"#$ -q {queue.name}" , file=sub)
            print(f"#$ -pe smp {procs}", file=sub)
            print(f"#$ -cwd", file=sub)

        # 1) Prepare an isolated directory inside the job scratch, for every QE run.
        print('', file=sub)
        print(f'module load {module}', file=sub)
        print('set OMP_NUM_THREADS=1', file=sub)
        print('ulimit -l unlimited', file=sub)
        print('', file=sub)
        print('# Prepare job scratch', file=sub)
        print('JOBDIR=$PWD', file=sub)
        print('SCRATCHDIR=$(mktemp -d "${TMPDIR:?}/scope.XXXXXX")', file=sub)
        print('cd "$SCRATCHDIR"', file=sub)
        print('echo "SCOPE: Running in $PWD"', file=sub)
        print(f'cp "$JOBDIR"/{shlex.quote(comp.inp_name)} .', file=sub)

        # 2) Reuse the keyword chain's wf for continuations or workflow-sharing Jobs.
        if wf_prepare: print(wf_prepare, file=sub)

        # 3) Run QE. Saved data will be copied even if QE does not converge.
        print('', file=sub)
        print('# Run QE', file=sub)
        if queue._environment.scheduler == 'slurm' and procs >= 128: print(f"{launch_command} < {comp.inp_name} > {comp.out_name} -pd .true.", file=sub)
        else:                                                        print(f"{launch_command} < {comp.inp_name} > {comp.out_name}", file=sub)

        # 4) Preserve the wf. Stop only this copy block on failure, then return the output.
        if wf_store: print(wf_store, file=sub)

        # 5) Return the text output even when preserving the wf failed.
        print('', file=sub)
        print('# Return the text output', file=sub)
        print(f'if ! cp -pr {shlex.quote(comp.out_name)} "$JOBDIR"; then', file=sub)
        print('    exit 1', file=sub)
        print('fi', file=sub)
        if wf_store: print('exit "$wf_status"', file=sub)

        os.chmod(comp.sub_path, 0o777)

#####################
## Other Functions ##
#####################
def gen_qe_wf_blocks(comp: object, debug: int=0) -> tuple[str, str]:
    """Return the wf preparation and storage shell blocks for a QE submission.

    Disabled storage returns two empty strings. Job mode prepares a saved wf
    only for continuations; workflow mode also prepares it for initial runs.
    Runtime helpers use the submitting Python environment. Storage stages and
    checks new data before replacing the previous wf, independently of QE exit status.
    Pseudopotential copies named in the XML are omitted from persistent storage.
    """
    import shlex
    import sys
    save_wf = getattr(comp._job.job_data, 'save_wf', False)
    if not save_wf: return '', ''
    reuse_wf = save_wf == 'workflow' or comp.is_update or comp.step > 1 or comp.run_number > 1

    wf_prepare = ''
    if reuse_wf:
        wf_prepare_code = f'import sys; from scope.software.quantum_espresso.qe_input import prepare_qe_wf; prepare_qe_wf(sys.argv[1], sys.argv[2], debug={max(debug - 1,0)})'
        wf_prepare = f'\n# Prepare the saved wf; keep the atomic guess if unavailable\n{shlex.quote(sys.executable)} -c {shlex.quote(wf_prepare_code)} {shlex.quote(comp.wf_path)} {shlex.quote(comp.inp_name)}'

    wf_check_code     = f'import sys; from scope.software.quantum_espresso.qe_input import check_qe_wf; sys.exit(not check_qe_wf(sys.argv[1], debug={max(debug - 1,0)}))'
    wf_check_command  = f'{shlex.quote(sys.executable)} -c {shlex.quote(wf_check_code)}'
    wf_remove_code    = f'import sys; from scope.software.quantum_espresso.qe_input import remove_qe_pseudopotentials; remove_qe_pseudopotentials(sys.argv[1], debug={max(debug - 1,0)})'
    wf_remove_command = f'{shlex.quote(sys.executable)} -c {shlex.quote(wf_remove_code)}'
    wf_store = f'''
# Save the wavefunction directory, replacing any previous saved copy
WF_NAME={shlex.quote(comp.wf_name)}
WF_DIR={shlex.quote(os.path.dirname(comp.wf_path))}
wf_status=0
if [ -d "$SCRATCHDIR/$WF_NAME" ]; then
    (
        set -e  # Stop this block if a file operation fails

        # Copy the new data before touching the previous copy
        wf_staging=$(mktemp -d "$WF_DIR/.scope-wf.XXXXXX")
        cp -pr "$SCRATCHDIR/$WF_NAME" "$wf_staging/"

        # Retain the previous wf if the new saved data is incomplete
        if ! {wf_check_command} "$wf_staging/$WF_NAME"; then
            rm -rf "$wf_staging"
            exit 0
        fi

        # SCOPE supplies pseudopotentials through pseudo_dir; omit local copies
        {wf_remove_command} "$wf_staging/$WF_NAME"

        # Keep the previous copy until its replacement is installed
        if [ -e "$WF_DIR/$WF_NAME" ]; then
            mv "$WF_DIR/$WF_NAME" "$wf_staging/previous"
        fi
        mv "$wf_staging/$WF_NAME" "$WF_DIR/$WF_NAME"
        rm -rf "$wf_staging"
    )
    wf_status=$?
else
    echo "SCOPE: QE save directory not found: $SCRATCHDIR/$WF_NAME" >&2
fi'''
    return wf_prepare, wf_store

def remove_qe_pseudopotentials(wf_path: str, debug: int=0) -> None:
    """Remove pseudopotential copies named in QE output XML from a staged .save folder.

    Only local basenames listed under output/atomic_species are removed, regardless
    of extension. Missing copies are ignored; XML and other files are retained.
    Call on the staged copy, leaving QE scratch intact. File/XML errors propagate
    so the submission script retains the previous saved wf if cleanup fails.
    """
    from xml.etree import ElementTree

    root = ElementTree.parse(os.path.join(wf_path, 'data-file-schema.xml')).getroot()
    removed = 0
    for item in root.findall('./{*}output/{*}atomic_species/{*}species/{*}pseudo_file'):
        name = os.path.basename((item.text or '').strip())
        if not name: continue
        path = os.path.join(wf_path, name)
        if os.path.isfile(path) or os.path.islink(path):
            os.remove(path)
            removed += 1
    if debug > 0: print(f"REMOVE_QE_PSEUDOPOTENTIALS: removed {removed} local copies from {wf_path}")

def check_qe_wf(wf_path: str, debug: int=0) -> bool:
    """Check the basic completeness of a collected QE .save directory.

    Read data-file-schema.xml and require readable, nonempty density and orbital
    files for every k-point/spin channel, in one consistent .dat or .hdf5 format.
    Require occup.txt for DFT+U and paw.txt for PAW when identified in the metadata.
    Pseudopotential copies are optional; SCOPE supplies them through pseudo_dir.
    Legacy/distributed saves are unsupported. Return only a Boolean; explanations
    are printed when debug > 0. Files are opened read-only and only one byte of
    each binary file is read. This does not validate binary contents, SCF convergence,
    the electronic state or compatibility with a receiving computation/QE build.
    """
    from xml.etree import ElementTree

    if wf_path is None or not os.path.isdir(wf_path):
        if debug > 0: print(f"CHECK_QE_WF: save directory not found: {wf_path}")
        return False

    try:
        root  = ElementTree.parse(os.path.join(wf_path, 'data-file-schema.xml')).getroot()
        bands = root.find('./{*}output/{*}band_structure')
        if root.tag.rsplit('}',1)[-1] != 'espresso' or bands is None: raise ValueError('QE output band metadata is missing')
        nks       = int(bands.findtext('{*}nks', default='0'))
        lsda      = bands.findtext('{*}lsda', default='').strip().lower()
        collected = bands.findtext('{*}wf_collected', default='').strip().lower()
        if nks < 1: raise ValueError('Invalid or missing k-point count')
        if lsda not in ['true', 'false', '1', '0']: raise ValueError('Invalid or missing spin-channel flag')
        if collected not in ['true', '1']: raise ValueError('Metadata does not identify collected wavefunctions')
        if len(bands.findall('{*}ks_energies')) != nks: raise ValueError('K-point metadata is incomplete')

        # QE XML counts physical k-points; collinear LSDA needs up and down files for each.
        prefixes = ['wfcup', 'wfcdw'] if lsda in ['true', '1'] else ['wfc']
        binaries = ['charge-density'] + [f'{prefix}{idx}' for prefix in prefixes for idx in range(1,nks+1)]
        auxiliary = []
        if root.find('./{*}output/{*}dft/{*}dftU') is not None: auxiliary.append('occup.txt')
        paw = root.findtext('./{*}output/{*}algorithmic_info/{*}paw', default='false').strip().lower()
        if paw not in ['true', 'false', '1', '0']: raise ValueError('Invalid PAW flag')
        if paw in ['true', '1']: auxiliary.append('paw.txt')
        for name in auxiliary:
            with open(os.path.join(wf_path, name), 'rb') as data:
                if not data.read(1): raise ValueError(f'{name} is empty')

        # A QE build reads one binary format; accepting a mixture would hide missing files.
        failures = []
        for extension in ['.dat', '.hdf5']:
            try:
                for name in binaries:
                    with open(os.path.join(wf_path, name+extension), 'rb') as data:
                        if not data.read(1): raise ValueError(f'{name+extension} is empty')
            except (OSError, ValueError) as exc:
                failures.append(str(exc))
                continue
            if debug > 0: print(f"CHECK_QE_WF: required files present and nonempty ({nks} k-points, {extension}): {wf_path}")
            return True
        if debug > 0: print(f"CHECK_QE_WF: no complete binary file set: {'; '.join(failures)}")
        return False
    except (OSError, ValueError, ElementTree.ParseError) as exc:
        if debug > 0: print(f"CHECK_QE_WF: {exc}")
        return False

def prepare_qe_wf(wf_path: str, inp_path: str, debug: int=0) -> bool:
    """Copy a saved QE wf beside a generated scratch input and enable file initialization.

    Return True only after copying, checking the copy and replacing the scratch
    input. On failure retain the original input and return False. Expected input
    uses the saved directory's prefix, outdir='./', and SCOPE's atomic defaults.
    Whitespace around assignment signs is allowed and preserved.
    The persistent wf is read-only; the scratch destination must not already exist.
    This checks file completeness, not binary contents or calculation compatibility.
    """
    import shutil
    import tempfile

    if not check_qe_wf(wf_path, debug=max(debug - 1,0)):
        if debug > 0: print(f"PREPARE_QE_WF: no complete saved wf; keeping the atomic guess: {wf_path}")
        return False
    inp_path = os.path.abspath(inp_path)
    wf_name  = os.path.basename(os.path.normpath(wf_path))
    target   = os.path.join(os.path.dirname(inp_path), wf_name)
    if os.path.lexists(target):
        if debug > 0: print(f"PREPARE_QE_WF: scratch destination already exists: {target}")
        return False

    copy_started = False
    temporary    = None
    try:
        with open(inp_path, 'r') as inp: text = inp.read()
        expected = {'prefix': f"'{os.path.splitext(wf_name)[0]}'", 'outdir': "'./'", 'startingpot': "'atomic'", 'startingwfc': "'atomic+random'"}
        lines = text.splitlines(keepends=True)
        settings = {}
        for idx, line in enumerate(lines):
            key, separator, value = line.partition('=')
            key = key.strip()
            if separator and key in expected:
                if key in settings or value.strip() != expected[key]: raise ValueError('Input does not use the expected prefix, outdir and atomic defaults')
                settings[key] = idx
        if settings.keys() != expected.keys(): raise ValueError('Input does not use the expected prefix, outdir and atomic defaults')
        copy_started = True
        shutil.copytree(wf_path, target)
        if not check_qe_wf(target, debug=max(debug - 1,0)): raise ValueError('Copied wf is incomplete')
        for key in ['startingpot', 'startingwfc']: lines[settings[key]] = lines[settings[key]].replace(expected[key], "'file'")
        text = ''.join(lines)
        with tempfile.NamedTemporaryFile(mode='w', dir=os.path.dirname(inp_path), delete=False) as inp:
            temporary = inp.name
            inp.write(text)
        os.replace(temporary, inp_path)
    except (OSError, ValueError) as exc:
        if copy_started: shutil.rmtree(target, ignore_errors=True)
        if debug > 0: print(f"PREPARE_QE_WF: {exc}; keeping the atomic guess")
        return False
    finally:
        if temporary is not None and os.path.isfile(temporary): os.remove(temporary)
    if debug > 0: print(f"PREPARE_QE_WF: initialized from {wf_path}")
    return True

def get_label_spin_pairs(source: object, debug: int=0):
    if not hasattr(source,"atoms"): source.set_atoms(debug=debug)
    pairs = []
    for at in source.atoms:
        if tuple([at.label,at.spin]) not in pairs: pairs.append(tuple([at.label,at.spin]))
    if debug > 0: print(f"GET_LABEL_SPIN_PAIRS. Pairs of label-spin found: {pairs}")
    return pairs

######
def get_pp(elem: str, path: str):
    import json
    try: 
        with open(path+"config.json", 'r') as cfg:
            data = json.load(cfg)
            elemdata = data[elem]
            pp_name = elemdata["filename"]
            cutoff_wfc = elemdata["cutoff_wfc"]
            cutoff_rho = elemdata["cutoff_rho"]
    except Exception as exc: 
        raise ValueError(f"GET_PP: error reading JSON or parsing element {elem}. Printing Exception \n {exc}")
    return pp_name, cutoff_wfc, cutoff_rho

################
### Check QE ###
################
def check_qe_launcher(module: str, scheduler: str='slurm', debug: int=0) -> dict:
    """Inspect QE's executable, MPI libraries, and possible external launchers.

    Returns:
        Static availability, probe results, MPI clues, launcher candidates, and
        warnings and a preferred launcher. Preference is not runtime validation.

    Open MPI prioritizes mpirun, then mpiexec, then direct Slurm launch. Intel
    MPI prioritizes Slurm PMI-2 when available. Other/unknown implementations
    prefer mpirun/mpiexec before direct Slurm launch; site settings may override.

    Only inspection commands run. Static executables and wrapper scripts can
    hide MPI dependencies; launcher version output does not establish a match
    with pw.x. No allocation, MPI ranks, or QE calculations are started.
    """
    from scope.classes_environment import run_module_checks

    commands = {
        'executable': 'command -v pw.x',
        'libraries': 'qe_executable=$(command -v pw.x) && command -v ldd >/dev/null && ldd "$qe_executable"',
        'mpirun': 'command -v mpirun && mpirun --version',
        'mpiexec': 'command -v mpiexec && mpiexec --version',
    }
    if scheduler == 'slurm': commands['slurm_mpi'] = 'srun --mpi=list'
    checks = run_module_checks(module, commands, debug=max(debug - 1, 0))
    # 0) Check that the module loads and pw.x is available.
    module_check  = checks.get('module')
    executable    = checks.get('executable')
    libraries     = checks.get('libraries')
    slurm_mpi     = checks.get('slurm_mpi')
    status        = 'inconclusive'
    mpi_libraries = []
    mpi_family    = None
    launchers     = []
    warnings      = []
    if module_check is None or not module_check.ok:
        warnings.append('Module check failed; verify Bash module setup.')
    else:
        if executable is not None:
            status = 'available' if executable.ok else 'unavailable'
            if not executable.ok: warnings.append('pw.x not found.')

        # 1) Identify MPI clues from dependencies, not the module name.
        if libraries is not None and libraries.ok:
            for line in libraries.stdout.splitlines():
                if any(word in line.lower() for word in ['libmpi', 'libmpich', 'libpmi']): mpi_libraries.append(line.strip())
            library_text = '\n'.join(mpi_libraries).lower()
            if 'openmpi' in library_text or 'open-mpi' in library_text: mpi_family = 'Open MPI'
            elif 'intel' in library_text and 'libmpi' in library_text: mpi_family = 'Intel MPI'
            elif 'mpich' in library_text: mpi_family = 'MPICH-family'
            if 'not found' in libraries.stdout: warnings.append('Unresolved shared libraries.')
            if not mpi_libraries: warnings.append('MPI libraries not visible; build type unknown.')
            elif mpi_family is None: warnings.append('MPI implementation unknown.')
        else: warnings.append('Dependencies unavailable; check ldd, static build, or wrapper.')

        # 2) Collect launcher candidates. Compatibility remains untested.
        for launcher in ['mpirun', 'mpiexec']:
            launcher_check = checks.get(launcher)
            if launcher_check is not None and launcher_check.ok: launchers.append(launcher)
        if slurm_mpi is not None and slurm_mpi.ok:
            plugins = slurm_mpi.stdout.split()
            if 'pmix' in plugins and mpi_family in [None, 'Open MPI', 'MPICH-family']: launchers.append('srun --mpi=pmix')
            if 'pmi2' in plugins and mpi_family in [None, 'Intel MPI', 'MPICH-family']: launchers.append('srun --mpi=pmi2')
        elif scheduler == 'slurm': warnings.append('Slurm MPI plugins unknown.')
        if not launchers: warnings.append('No launcher candidates found.')
    if 'session' in checks: warnings.append('Inspection incomplete.')

    # 3) Rank candidates and choose a default only when pw.x is available.
    if mpi_family == 'Intel MPI' and 'srun --mpi=pmi2' in launchers:
        launchers.remove('srun --mpi=pmi2')
        launchers.insert(0, 'srun --mpi=pmi2')
    preferred_launcher = launchers[0] if status == 'available' and launchers else None
    if debug > 0:
        for name, check in checks.items(): print(f'CHECK_QE: {name}: {"OK" if check.ok else "FAILED"}')
    return {'module': module, 'status': status, 'runtime_validated': False, 'checks': checks, 'mpi_libraries': mpi_libraries, 'mpi_family': mpi_family, 'launcher_candidates': launchers, 'launcher': preferred_launcher, 'warnings': warnings}
