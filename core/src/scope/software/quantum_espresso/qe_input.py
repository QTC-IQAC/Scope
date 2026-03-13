#!/usr/bin/env python3
import os
from scope.classes_state       import find_state
from scope.classes_environment import set_user
from scope.operations.dicts_and_lists import where_in_array
from scope.other               import get_metal_idxs, get_metal_species
from scope.elementdata         import ElementData 
from scope.software.quantum_espresso.qe_functions import get_pp, get_qe_data
elemdatabase = ElementData()

#######################
def gen_qe_input(comp: object, debug: int=0):

    ######################
    ### INITIAL CHECKS ###
    ######################
    ## 1-Verify Information in initial state
    assert hasattr(comp.qc_data,"istate"), f"istate = {comp.qc_data.istate} not found in comp.qc_data"
    exists, istate    = find_state(comp.source, comp.qc_data.istate)
    assert exists,                         f"istate = {comp.qc_data.istate} does not exist"
    assert hasattr(istate,"labels"),       f"istate = {comp.qc_data.istate} doesn't have labels"
    assert hasattr(istate,"coord"),        f"istate = {comp.qc_data.istate} doesn't have coordinates"
    if not hasattr(istate,"atoms"):        istate.get_atoms()

    ## 2-Cell or Molecule?
    if   comp.source.object_type == "cell":   system_type = "cell"
    elif comp.source.object_type == "specie": system_type = "molecule"
    else: print(f"GEN_QE_input: unrecognised {comp.source.object_type=}")

    ## 3-Determines the PP_Library path
    if not hasattr(comp.qc_data,"pp_library"): f"PP_Library could not be found. Please set it in the qc_data.section of the Job"
    import scope ## Trick to retrieve the relative path
    PP_path = os.path.abspath(scope.__file__).replace("__init__.py","software/quantum_espresso/PP_Libraries/")
    if   comp.qc_data.pp_library.lower() == "efficiency":  PP_path += "Efficiency/"
    elif comp.qc_data.pp_library.lower() == "precision":   PP_path += "Precision/"
    elif comp.qc_data.pp_library.lower() == "vanderbilt":  PP_path += "Vanderbilt_USPP/"
    else: raise ValueError(f"GEN_QE_INPUT: could not recognise PP_library: {comp.qc_data.pp_library}")

    ## Checks requisites
    if system_type == "cell":     
        assert hasattr(istate,"cell_vector"),  f"GEN_QE_input: {comp.qc_data.istate=} doesn't have cell vectors"
        assert istate.charge == int(0),        f"GEN_QE_input: {istate.charge=} must be 0. The STATE charge must be 0"

    if debug >= 1 : print("GEN_QE_INPUT: system_type:", system_type)
    if debug >= 1 : print("GEN_QE_INPUT: creating input in path:", comp.inp_path)

    #########################
    ### DETERMINE SPECIES ###
    #########################
    metal_indices            = get_metal_idxs(istate.labels)     ## Indices of metal atoms in the structure 
    metal_species            = get_metal_species(istate.labels)  ## 
    species                  = get_qe_data(istate)
    nspecies                 = len(species)

    if debug >= 1: 
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
        print(f"    calculation='{comp.qc_data.jobtype}'", file=inp)
        print( "    restart_mode='from_scratch'", file=inp)
        print(f"    pseudo_dir ='{PP_path}'", file=inp)
        print( "    disk_io='low'", file=inp)
        print(f"    prefix='{comp.name}'", file=inp)

        ## Print forces for finite differences computation.
        if hasattr(comp.qc_data,"print_forces"):
            if comp.qc_data.print_forces and comp.qc_data.jobtype == "scf": print(f"    tprnfor = .true.", file=inp)

        ## Keywords for Optimization-related computations 
        if comp.qc_data.jobtype == "opt" or comp.qc_data.jobtype == "relax" or comp.qc_data.jobtype == "vc-relax":
            print(f"    wf_collect = .true.", file=inp)
            print(f"    tprnfor = .true.", file=inp)
            print(f"    nstep = 300", file=inp)
            print(f"    forc_conv_thr = {comp.qc_data.forc_conv}", file=inp)
        print("/", file=inp)

        #/////////////////////
        #// System control ///
        #/////////////////////
        min_cowfc = 25
        min_corho = 200
        ### Determines cutoffs using the info contained in the Pseudopotentials
        for idx, spec in enumerate(species):
            pp, cutoff_wfc, cutoff_rho = get_pp(spec[0], PP_path)
            if cutoff_wfc > min_cowfc: min_cowfc = cutoff_wfc # updates
            if cutoff_rho > min_corho: min_corho = cutoff_rho # updates

        if comp.qc_data.jobtype == "vc-relax": min_cowfc *= 2; min_corho *= 2  ## In vc-relax it is convenient to minimize pulay stress

        print(" &system", file=inp)
        if   system_type == "molecule": print(f"    ibrav=1, celldm(1)={comp.qc_data.cubeside}", file=inp)
        elif system_type == "cell":     print(f"    ibrav=0,", file=inp)
        print(f"    nat={istate.natoms}, ntyp={nspecies}, ecutwfc={int(min_cowfc)}, ecutrho={float(min_corho)}", file=inp)
        print(f"    nspin=2,", file=inp)
        print(f"    tot_charge={istate.charge}", file=inp)
        print(f"    tot_magnetization={istate.spin}", file=inp)
  
        ## Starting Magnetization
        for idx, spec in enumerate(species):
            if spec[1] != 0: print(f"    starting_magnetization({idx})={spec[1]}", file=inp)

        ## Hubbard: Here it is assumed that the Hubbard U term will apply to all metals
        if comp.qc_data.is_hubbard and comp.qc_data.version <= 7.0:
            print("    lda_plus_u=.true.,", file=inp)
            for idx, spec in enumerate(species):
                if spec[1] != 0: print(f"    Hubbard_U({idx})={comp.qc_data.uterm}", file=inp)

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
        print(f"    conv_thr = {comp.qc_data.elec_conv}", file=inp)
        print(f"    mixing_beta = {comp.qc_data.mix_beta}", file=inp)
        print("/", file=inp)
   
        #///////////////////
        #// Ions control ///
        #///////////////////
        if comp.qc_data.jobtype == "opt" or comp.qc_data.jobtype == "relax" or comp.qc_data.jobtype == "vc-relax":
            print(" &ions", file=inp)
            print("    upscale=10", file=inp)
            print("    ion_dynamics='bfgs'", file=inp)
            print("    trust_radius_ini= 0.1D0", file=inp)
            print("/", file=inp)
        
        #///////////////////
        #// Cell control ///
        #///////////////////
        if comp.qc_data.jobtype == "vc-relax":
            print(" &cell", file=inp)
            print("    cell_dynamics='bfgs'", file=inp)
            print("    cell_dofree='all'", file=inp)
            print("    cell_factor= 1.2D0", file=inp)
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
        
        #///////////////////////////////////////////////////////////////////////////
        #// Atom Coords, taken from the state object. Using the provided istate  ///
        #///////////////////////////////////////////////////////////////////////////
        print("ATOMIC_POSITIONS angstrom", file=inp)
        # First, it prints the metal atoms
        for idx, at in enumerate(istate.atoms):
            if idx in metal_indices:
                if at.spin == 0: l = at.label
                else:            l = at.get_decorated_label(typ="spin") 
                print(f"{l:4}        {at.coord[0]:12.6f}   {at.coord[1]:12.6f}   {at.coord[2]:12.6f}", file=inp)
        # Then the rest
        for idx, at in enumerate(istate.atoms):
            if idx not in metal_indices:
                if at.spin == 0: l = at.label
                else:            l = at.get_decorated_label(typ="spin") 
                print(f"{l:4}        {at.coord[0]:12.6f}   {at.coord[1]:12.6f}   {at.coord[2]:12.6f}", file=inp)
        print("K_POINTS gamma", file=inp)

        #/////////////////////////////
        #// HUBBARD data after 7.1 ///
        #/////////////////////////////
        if comp.qc_data.is_hubbard and comp.qc_data.version > 7.0:
            print("HUBBARD (atomic)", file=inp)
            for idx, spec in enumerate(species):
                if spec[1] != 0: print(f"U {spec[0]}-3d {comp.qc_data.uterm}", file=inp)

###################################################
def gen_qe_subfile(comp: object, queue: object, module: str, procs: int=1, exe: str="pw.x", version: float=7.0): 
    version = str(version)
    with open(comp.sub_path, 'w+') as sub:
        if queue._environment.scheduler == 'slurm':
            print(f"#!/bin/bash", file=sub)
            print(f"#SBATCH -J {comp.name}", file=sub)
            print(f"#SBATCH -e {comp.name}.stderr", file=sub)
            print(f"#SBATCH -o {comp.name}.stdout", file=sub)
            print(f"#SBATCH -p {queue.name}", file=sub)
            print(f"#SBATCH --nodes=1", file=sub)
            print(f"#SBATCH --ntasks={procs}", file=sub)
            print(f"#SBATCH --mem-per-cpu=1900MB", file=sub)
            print(f"#SBATCH --time={queue.time_limit}", file=sub)
            print(f"", file=sub)
            print(f"module load {module}", file=sub)
            print(f"", file=sub)
            print(f"set OMP_NUM_THREADS=1", file=sub)
            print(f"ulimit -l unlimited", file=sub)
            print(f"", file=sub)
            print(f"JOBDIR=$PWD", file=sub)
            print(f"cd $TMPDIR", file=sub)
            print(f"cp $JOBDIR/{comp.inp_name} .", file=sub)
            if procs >= 128: print(f"srun pw.x < {comp.inp_name} > {comp.out_name} -pd .true.", file=sub)
            else:            print(f"srun pw.x < {comp.inp_name} > {comp.out_name}", file=sub)
            print(f"cp -pr {comp.out_name} $JOBDIR", file=sub)

        elif queue._environment.scheduler == 'sge':
            print(f"#!/bin/bash", file=sub)
            print(f"#$ -N {comp.name}", file=sub)
            print(f"#$ -o {comp.name}.stdout", file=sub)
            print(f"#$ -e {comp.name}.stderr", file=sub)
            print(f"#$ -q {queue.name}" , file=sub)
            print(f"#$ -pe smp {procs}", file=sub)
            print(f"#$ -cwd", file=sub)
            print(f"", file=sub)
            #print(f"source /etc/profile.d/modules.csh", file=sub)
            #print(f"source $HOME/.bashrc", file=sub)
            #print(f". /etc/profile", file=sub)
            print(f"module load {module}", file=sub)
            print(f"", file=sub)
            print(f"set OMP_NUM_THREADS=1", file=sub)
            print(f"ulimit -l unlimited", file=sub)
            print(f"", file=sub)
            print(f"JOBDIR=$PWD", file=sub)
            print(f"cd $TMPDIR", file=sub)
            print(f"cp $JOBDIR/{comp.inp_name} .", file=sub)
            print(f"mpirun -np {procs} pw.x < {comp.inp_name} > {comp.out_name}", file=sub)
            print(f"cp -pr {comp.out_name} $JOBDIR", file=sub)

        os.chmod(comp.sub_path, 0o777)

###################################################
