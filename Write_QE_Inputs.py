#!/usr/bin/env python3
import sys
import os
import numpy as np
import pickle
from collections import Counter
from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal

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
def gen_QE_input(comp, debug: int=0):

    gmol = comp._job._recipe.subject

    PP_Library = comp.qc_data.pseudo
    cluster = set_cluster()

    if hasattr(gmol,"cellparam"): system_type = "cell"
    else: system_type = "molecule"

    if system_type == "cell":     
        if hasattr(gmol,"natoms"): natoms = gmol.natoms
        else:                      natoms = len(gmol.labels) 
    if system_type == "molecule":  natoms = gmol.natoms

    if debug >= 1 : print(system_type)
    if debug >= 1 : print("Creating input in path:", comp.inp_path)

    # Corrects PP_Library path if necessary
    if PP_Library[-1] != '/': PP_Library += '/'

    #########################
    ### DETERMINE SPECIES ###
    #########################
    comp.spin_config.get_QE_data()
    comp.spin_config.get_info()
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
        print("GEN_QE_INPUT: ismagnetic:    {comp.spin_config.ismagnetic}")
        print("GEN_QE_INPUT: magnetization: {comp.spin_config.total_magnetization}")
    
    ########################### 
    ### WRITES INPUT: PART 1 ##
    ########################### 
    with open(comp.inp_path, 'w+') as inp:
        print(" &control", file=inp)
        print(f"    calculation='{comp.qc_data.jobtype}'", file=inp)
        print( "    restart_mode='from_scratch'", file=inp)
        print(f"    pseudo_dir ='{PP_Library}'", file=inp)
        print( "    disk_io='low'", file=inp)

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
            print(f"    forc_conv_thr = 1.0D-5", file=inp)

        print("/", file=inp)

        #/////////////////////
        #// System control ///
        #/////////////////////
        print(" &system", file=inp)
        if system_type == "molecule": print(f"    ibrav=1, celldm(1)={cubeside}", file=inp)
        elif system_type == "cell":   print(f"    ibrav=0,", file=inp)
        print(f"    nat={natoms}, ntyp={nspecies}, ecutwfc={int(comp.qc_data.cutoff)}, ecutrho={float(comp.qc_data.cutoff)*8}", file=inp)
        if not hasattr(comp.spin_config,"ismagnetic"): comp.spin_config.get_total_magnetization()
        if comp.spin_config.ismagnetic: print(f"    nspin=2,", file=inp)
        else:          print(f"    nspin=1,", file=inp)
        
        if system_type == "molecule" and hasattr(gmol, 'totcharge'):
            print(f"    tot_charge={gmol.totcharge}", file=inp)
        else:
            print("    tot_charge=0", file=inp)
        
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
        if comp.qc_data.is_hubbard:
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
        print(" &electrons", file=inp)
        print("    diagonalization='david'", file=inp)
        print("    electron_maxstep=250", file=inp)
        print("    conv_thr = 1.0e-5", file=inp)
        print("    mixing_beta = 0.20", file=inp)
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
            print(f"{gmol.cellvec[0][0]:14.6f} {gmol.cellvec[0][1]:14.6f} {gmol.cellvec[0][2]:14.6f}", file=inp)
            print(f"{gmol.cellvec[1][0]:14.6f} {gmol.cellvec[1][1]:14.6f} {gmol.cellvec[1][2]:14.6f}", file=inp)
            print(f"{gmol.cellvec[2][0]:14.6f} {gmol.cellvec[2][1]:14.6f} {gmol.cellvec[2][2]:14.6f}", file=inp)
         
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
        #// Atom Coords, taken from the mol object. Using the provided coord_tag ///
        #///////////////////////////////////////////////////////////////////////////
        print("ATOMIC_POSITIONS angstrom", file=inp)
        cor_tag = getattr(gmol,comp.qc_data.coord_tag)
        for idx, l in enumerate(gmol.labels):
            pointer = 0
            if idx in metal_indices and len(comp.spin_config.atomic_spins) != 0: 
                label = comp.spin_config.atomic_spins[pointer].get_mod_label() 
                pointer += 1
            else: 
                label = l
            print(f"{label:4}        {cor_tag[idx][0]:12.6f}   {cor_tag[idx][1]:12.6f}   {cor_tag[idx][2]:12.6f}", file=inp)
        print("K_POINTS gamma", file=inp)

###################################################
def gen_QE_subfile(comp: object, procs: int=1, queue: str="iqtc09", cluster: str=set_cluster(), user: str=set_user(), exe: str="pw.x", version: str="6.4.1"): 
    if 'portal' in cluster:
        with open(comp.sub_path, 'w+') as sub:
            print(f"#!/bin/bash", file=sub)
            print(f"#$ -N {comp.refcode}{comp.suffix}", file=sub) 
            print(f"#$ -pe smp {procs}", file=sub)
            print(f"#$ -cwd", file=sub)
            print(f"#$ -o /home/{user}/x-stds/{comp.refcode}{comp.suffix}.stdout", file=sub)
            print(f"#$ -e /home/{user}/x-stds/{comp.refcode}{comp.suffix}.stderr", file=sub)
            print(f"#$ -S /bin/bash", file=sub)
            print(f"#$ -q {queue}.q" , file=sub)
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

    elif 'login' in cluster or 'csuc' in cluster:
        with open(comp.sub_path, 'w+') as sub:
            print(f"#!/bin/bash", file=sub)
            print(f"#SBATCH -J {comp.refcode}{comp.suffix}", file=sub)
            print(f"#SBATCH -e /scratch/{user}/std_files/{comp.refcode}{comp.suffix}.stderr", file=sub)
            print(f"#SBATCH -o /scratch/{user}/std_files/{comp.refcode}{comp.suffix}.stdout", file=sub)
            print(f"#SBATCH -p {queue}", file=sub)
            print(f"#SBATCH --nodes=1", file=sub)
            print(f"#SBATCH --ntasks={procs}", file=sub)
            print(f"#SBATCH --time=10-0", file=sub)
            print(f"#SBATCH -x pirineus", file=sub)
            print(f"", file=sub)
            print(f"module load apps/quantumespresso/6.4.1", file=sub)
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







