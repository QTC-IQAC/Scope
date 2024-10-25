#!/usr/bin/env python3
import sys
import os
import numpy as np
import pickle
from collections import Counter
from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal

from Scope.Classes_State import find_state
from Scope.Environment import set_user, set_cluster
from Scope.Elementdata import ElementData 
from Scope.Other import get_metal_idxs, get_metal_species, where_in_array
elemdatabase = ElementData()

#######################
def get_pp(ppfolder, elem):

    Found = False
    while not Found:
        for file in os.listdir(ppfolder):
            with open(ppfolder+file, 'rb') as pp:
                splitpp = file.split(".")
                if (len(splitpp) == 3):
                    element=splitpp[0]
                    functional=splitpp[1]
                    typ=splitpp[2]
                elif (len(splitpp) == 4):
                    element=splitpp[0]
                    functional=splitpp[1]
                    typ1=splitpp[2]
                    typ2=splitpp[3]
                else:
                    print("Can't interpret pp name for", file)

                if (element == elem):
                    Found = True
                    sendpp = file
                    break

        if file == os.listdir(ppfolder)[-1] and not Found:
            sendpp = "Not found"
            break

    return sendpp

#######################
def gen_QE_input(comp: object, environment: object, debug: int=0):

    gmol = comp._job._recipe.subject

    ## 1-Change some variable names to simplify calls
    if not hasattr(environment,"PP_Library"): f"PP_Library could not be found. Please set it in the environment"
    PP_Library = environment.PP_Library
    cluster    = environment.cluster

    ## 2a-Verify Information in State
    assert hasattr(comp.qc_data,"istate"), f"istate = {comp.qc_data.istate} not found in comp.qc_data"
    exists, istate    = find_state(gmol, comp.qc_data.istate)
    assert exists, f"istate = {comp.qc_data.istate} does not exist"
    assert hasattr(istate,"labels"), f"istate = {comp.qc_data.istate} doesn't have labels"
    assert hasattr(istate,"coord"),  f"istate = {comp.qc_data.istate} doesn't have coordinates"

    ## 2b-Cell or Molecule?
    if hasattr(gmol,"cellparam"): system_type = "cell"
    else:                         system_type = "molecule"
    
    ## Charge of System
    if system_type == "cell":
        if hasattr(gmol,"type"): 
            if gmol.type == "perxyz": system_charge = gmol.totcharge     
            else:                     system_charge = int(0)
        else:                         system_charge = int(0)
    else:                             system_charge = gmol.totcharge

    if system_type == "cell":     
        assert hasattr(istate,"cellvec"),  f"istate = {comp.qc_data.istate} doesn't have cell vectors"
        if hasattr(istate,"natoms"): natoms = istate.natoms
        else:                        natoms = len(istate.labels) 
    if system_type == "molecule":    natoms = istate.natoms

    if debug >= 1 : print("GEN_QE_INPUT: system_type:", system_type)
    if debug >= 1 : print("GEN_QE_INPUT: creating input in path:", comp.inp_path)

    #########################
    ### DETERMINE SPECIES ###
    #########################
    comp.spin_config.get_QE_data()
    metal_indices = get_metal_idxs(gmol.labels)
    metal_species = get_metal_species(gmol.labels)
    elems = comp.spin_config.elems
    nelems = comp.spin_config.nelems

    species = []
    for idx, elem in enumerate(elems):
        if elem in metal_species:
            added = []
            for un in comp.spin_config.magn_uniques:
                if elem == un[0] and un[1] not in added:
                    species.append(str(un[0])+str(un[1]))
                    added.append(un[1])
        else:
            species.append(elem)
    nspecies = len(species)

    if debug >= 1: 
        print("GEN_QE_INPUT: Received cell with this data:")
        print(f"    elems:    {elems}")
        print(f"    species:  {species}")
        print(f"    metal_species: {metal_species}")
        print(f"    metal_indices: {metal_indices}")
        if len(comp.spin_config.atomic_spins) > 0: 
            print("GEN_QE_INPUT: Also received spin configuration with this data:")
            print(f"    magn_pairs:    {comp.spin_config.magn_pairs}")
            print(f"    magn_uniques:  {comp.spin_config.magn_uniques}")
        print(f"GEN_QE_INPUT: ismagnetic:    {comp.spin_config.ismagnetic}")
        print(f"GEN_QE_INPUT: magnetization: {comp.spin_config.total_magnetization}")
    
    ########################### 
    ### WRITES INPUT: PART 1 ##
    ########################### 
    with open(comp.inp_path, 'w+') as inp:
        print(" &control", file=inp)
        print(f"    calculation='{comp.qc_data.jobtype}'", file=inp)
        print( "    restart_mode='from_scratch'", file=inp)
        print(f"    pseudo_dir ='{PP_Library}'", file=inp)
        print( "    disk_io='low'", file=inp)

        if 'uam' in cluster: print( "    max_seconds=61200", file=inp)
        if 'portal' in cluster: print(f"    outdir='/scratch/g4vela/QE_WFC'", file=inp)
        elif 'csuc' in cluster or 'login' in cluster: print(f"    outdir='/scratch/svela/QE_WFC/'", file=inp) 

        print(f"    prefix='{comp.refcode}{comp.suffix}'", file=inp)

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
        print(" &system", file=inp)
        if system_type == "molecule": print(f"    ibrav=1, celldm(1)={comp.qc_data.cubeside}", file=inp)
        elif system_type == "cell":   print(f"    ibrav=0,", file=inp)
        print(f"    nat={natoms}, ntyp={nspecies}, ecutwfc={int(comp.qc_data.cutoff)}, ecutrho={float(comp.qc_data.cutoff)*8}", file=inp)
        if not hasattr(comp.spin_config,"ismagnetic"): comp.spin_config.get_total_magnetization()
        if comp.spin_config.ismagnetic: print(f"    nspin=2,", file=inp)
        else:          print(f"    nspin=1,", file=inp)

        print(f"    tot_charge={system_charge}", file=inp)
        
        #if system_type == "molecule" and hasattr(gmol, 'totcharge'):
        #    print(f"    tot_charge={gmol.totcharge}", file=inp)
        #else:
        #    print("    tot_charge=0", file=inp)
        
        ## Total Magnetization
        #tot_magn = 0      
        #if len(spin_config) != 0:
        #    for t in spin_config:
        #        tot_magn += t[1]
        if comp.spin_config.ismagnetic: print(f"    tot_magnetization={comp.spin_config.total_magnetization}", file=inp)
  
        ## Starting Magnetization
        for idx, u in enumerate(comp.spin_config.magn_uniques):
            if u[1] != 0: print(f"    starting_magnetization({where_in_array(elems,u[0])[0]+1})={u[1]}", file=inp)

        ## Hubbard: Here it is assumed that the Hubbard U term will apply to all metals
        if comp.qc_data.is_hubbard and comp.qc_data.version <= 7.0:
            print("    lda_plus_u=.true.,", file=inp)
            for idx, u in enumerate(comp.spin_config.magn_uniques):
                print(f"    Hubbard_U({where_in_array(elems,u[0])[0]+1})={comp.qc_data.uterm}", file=inp)

        if comp.qc_data.is_grimme:
            print("    vdw_corr='grimme-d3'", file=inp)
            print("    dftd3_version=4", file=inp)
            
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
            print(f"{istate.cellvec[0][0]:14.6f} {istate.cellvec[0][1]:14.6f} {istate.cellvec[0][2]:14.6f}", file=inp)
            print(f"{istate.cellvec[1][0]:14.6f} {istate.cellvec[1][1]:14.6f} {istate.cellvec[1][2]:14.6f}", file=inp)
            print(f"{istate.cellvec[2][0]:14.6f} {istate.cellvec[2][1]:14.6f} {istate.cellvec[2][2]:14.6f}", file=inp)
         
        #///////////////////
        #// Atom Species ///
        #///////////////////
        #if len(spin_config) == 0: elems = list(set(gmol.labels))   ## This list of elements must be updated to include the different magnetization labels (eg. Fe1, Fe2)
        #else:                     elems = list(t[0] for t in uniques)
        print("ATOMIC_SPECIES", file=inp)
        for idx, spec in enumerate(species):
            if spec[-1].isdigit(): label = spec[:-1]
            else: label = spec
            weight = elemdatabase.elementweight[label]
            pp = get_pp(PP_Library, label)
            print(f"{spec} {weight:6.4f} {pp}", file=inp)
        
        #///////////////////////////////////////////////////////////////////////////
        #// Atom Coords, taken from the state object. Using the provided istate  ///
        #///////////////////////////////////////////////////////////////////////////
        print("ATOMIC_POSITIONS angstrom", file=inp)
        for idx, l in enumerate(istate.labels):
            pointer = 0
            if idx in metal_indices and len(comp.spin_config.atomic_spins) != 0: 
                label = comp.spin_config.atomic_spins[pointer].get_mod_label() 
                pointer += 1
            else: 
                label = l
            print(f"{label:4}        {istate.coord[idx][0]:12.6f}   {istate.coord[idx][1]:12.6f}   {istate.coord[idx][2]:12.6f}", file=inp)
        print("K_POINTS gamma", file=inp)

        #/////////////////////////////
        #// HUBBARD data after 7.1 ///
        #/////////////////////////////
        if comp.qc_data.is_hubbard and comp.qc_data.version > 7.0:
            print("HUBBARD (atomic)", file=inp)
            for idx, u in enumerate(comp.spin_config.magn_uniques):
                print(f"U {u[0]}-3d {comp.qc_data.uterm}", file=inp)

###################################################
def gen_QE_subfile(comp: object, queue: object, procs: int=1, exe: str="pw.x", version: str="6.4.1"): 

    version = str(version)
    cluster = queue._environment.cluster
    user    = queue._environment.user
    group   = queue._environment.group

    if 'portal' in cluster:
        with open(comp.sub_path, 'w+') as sub:
            print(f"#!/bin/bash", file=sub)
            print(f"#$ -N {comp.refcode}{comp.suffix}", file=sub) 
            print(f"#$ -pe smp {procs}", file=sub)
            print(f"#$ -cwd", file=sub)
            print(f"#$ -o /home/{user}/x-stds/{comp.refcode}{comp.suffix}.stdout", file=sub)
            print(f"#$ -e /home/{user}/x-stds/{comp.refcode}{comp.suffix}.stderr", file=sub)
            print(f"#$ -S /bin/bash", file=sub)
            print(f"#$ -q {queue.name}" , file=sub)
            print(f"", file=sub)
            print(f"source /etc/profile.d/modules.csh", file=sub)

            #if version == "6.4.1": print(f"module load espresso/6.4.1_ompi", file=sub)
            ## Only one version implemented 
            print(f"module load espresso/6.4.1_ompi", file=sub)

            print(f"", file=sub)
            print(f"set OMP_NUM_THREADS=1", file=sub)
            print(f"ulimit -l unlimited", file=sub)
            print(f"", file=sub)
            print(f"WORKDIR=$PWD", file=sub)
            print(f"cd $TMPDIR", file=sub)
            print(f"cp $WORKDIR/{comp.inp_name} .", file=sub)
            print(f"mpirun -np {procs} pw.x < {comp.inp_name} > {comp.out_name}", file=sub)
            print(f"cp -pr {comp.out_name} $WORKDIR", file=sub)
            os.chmod(comp.sub_path, 0o777)

    elif cluster == "login2" or cluster == "login3":
        if   cluster == "login2": storage = scratch
        elif cluster == "login3": storage = data
        with open(comp.sub_path, 'w+') as sub:
            print(f"#!/bin/bash", file=sub)
            print(f"#SBATCH -J {comp.refcode}{comp.suffix}", file=sub)
            print(f"#SBATCH -e /{storage}/{user}/std_files/{comp.refcode}{comp.suffix}.stderr", file=sub)
            print(f"#SBATCH -o /{storage}/{user}/std_files/{comp.refcode}{comp.suffix}.stdout", file=sub)
            print(f"#SBATCH -p {queue.name}", file=sub)
            print(f"#SBATCH --nodes=1", file=sub)
            print(f"#SBATCH --ntasks={procs}", file=sub)
            print(f"#SBATCH --time=7-0", file=sub)
            print(f"", file=sub)
            if   version == "6.4.1": print(f"module load apps/quantumespresso/6.4.1", file=sub)
            elif version == "7.0":   print(f"module load quantumespresso/7.0", file=sub)
            elif version == "7.2":   print(f"module load quantumespresso/7.2", file=sub)
            else:                    print(f"UNKNOWN VERSION",version)
            print(f"", file=sub)
            print(f"set OMP_NUM_THREADS=1", file=sub)
            print(f"ulimit -l unlimited", file=sub)
            print(f"", file=sub)
            print(f"WORKDIR=$PWD", file=sub)
            print(f"cd $TMPDIR", file=sub)
            print(f"cp $WORKDIR/{comp.inp_name} .", file=sub)
            print(f"srun pw.x < {comp.inp_name} > {comp.out_name}", file=sub)
            print(f"cp -pr {comp.out_name} $WORKDIR", file=sub)
            os.chmod(comp.sub_path, 0o777)
     
    elif 'uam' in cluster:
        project = 'ub100'
        with open(comp.sub_path, 'w+') as sub:
            print(f"#!/bin/bash", file=sub)
            print(f"#SBATCH -J {comp.refcode}{comp.suffix}", file=sub)
            print(f"#SBATCH -e /scratch/{project}/std_files/{comp.refcode}{comp.suffix}.stderr", file=sub)
            print(f"#SBATCH -o /scratch/{project}/std_files/{comp.refcode}{comp.suffix}.stdout", file=sub)
            print(f"#SBATCH -p {queue.name}", file=sub)
            print(f"#SBATCH --exclude=cibeles3-05", file=sub)
            print(f"#SBATCH -A ub100_serv", file=sub)
            print(f"#SBATCH --nodes=1", file=sub)
            print(f"#SBATCH --ntasks={procs}", file=sub)
            print(f"", file=sub)
            print(f"module load espresso/6.5", file=sub)
            print(f"export SLURM_MPI_TYPE=pmi2", file=sub)
            print(f"", file=sub)
            print(f"set OMP_NUM_THREADS=1", file=sub)
            print(f"ulimit -l unlimited", file=sub)
            print(f"", file=sub)
            print(f"JOBDIR=$PWD", file=sub)
            print(f'export RUNDIR="/temporal/{user}/jobs/$SLURM_JOBID"', file=sub)
            print(f"mkdir -p $RUNDIR", file=sub)
            print(f"cd $RUNDIR", file=sub)
            print(f"", file=sub)
            print(f"cp -i $JOBDIR/{comp.inp_name} .", file=sub)
            print(f"timeout 71h srun pw.x < {comp.inp_name} > {comp.out_name}", file=sub)
            print(f"cp -pr {comp.out_name} $JOBDIR/", file=sub)
            os.chmod(comp.sub_path, 0o777)

    ## Cesga
    elif 'login' in cluster and group == 'csic':
        with open(comp.sub_path, 'w+') as sub:
            print(f"#!/bin/bash", file=sub)
            print(f"#SBATCH -J {comp.refcode}{comp.suffix}", file=sub)
            print(f"#SBATCH -e /scratch/{user}/std_files/{comp.refcode}{comp.suffix}.stderr", file=sub)
            print(f"#SBATCH -o /scratch/{user}/std_files/{comp.refcode}{comp.suffix}.stdout", file=sub)
            print(f"#SBATCH -p {queue.name}", file=sub)
            print(f"#SBATCH --nodes=1", file=sub)
            print(f"#SBATCH --ntasks={procs}", file=sub)
            print(f"#SBATCH --time=10-0", file=sub)
            print(f"#SBATCH -x pirineus", file=sub)
            print(f"", file=sub)
            print(f"module load quantumespresso/6.5", file=sub)
            print(f"", file=sub)
            print(f"set OMP_NUM_THREADS=1", file=sub)
            print(f"ulimit -l unlimited", file=sub)
            print(f"", file=sub)
            print(f"WORKDIR=$PWD", file=sub)
            print(f"cd $TMPDIR", file=sub)
            print(f"cp $WORKDIR/{comp.inp_name} .", file=sub)
            print(f"srun pw.x < {comp.inp_name} > {comp.out_name}", file=sub)
            print(f"cp -pr {comp.out_name} $WORKDIR", file=sub)
            os.chmod(comp.sub_path, 0o777)

    else: print("WRITE_QE_INPUTS: Cluster not recognized")
        
###################################################







