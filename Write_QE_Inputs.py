#!/usr/bin/env python3
import sys
import os
import numpy as np
import pickle
from collections import Counter
from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal
from cell2mol.elementdata import ElementData 

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
def gen_QE_input(gmol: object, spin_config: list, path: str, name: str, suffix: str, extension: str, PP_Library: str, job_type: str="scf", isHubbard: bool=False, isGrimme: bool=True, U: str="2.27", cutoff: str="70", cubeside: str="70.0", debug: int=0):
 
    if hasattr(gmol,"cellparam"): system_type = "cell"
    else: system_type = "molecule"

    if system_type == "cell": natoms = len(gmol.labels) 
    if system_type == "molecule": natoms = gmol.natoms

    if debug > 0 : print(system_type)
    if debug > 0 : print("Creating input in path:", path+name+suffix+extension)

    ######################################
    ### NECESSARY TO DETERMINE SPECIES ###
    ######################################
    if len(spin_config) == 0:
        elems = list(set(gmol.labels)) 
        nspecies = len(list(set(gmol.labels))) 
    else:
        pairs = list(set(spin_config))  ## tuples of element_labels + magnetization
        elems = list(set(t[0] for t in pairs))
        uniques = list(set(pairs)) 
        nspecies = len(pairs)

    metal_species = []
    for idx, l in enumerate(elems):
        if l[-1].isdigit(): label = l[:-1]
        else: label = l
        if (elemdatabase.elementblock[label] == 'd' or elemdatabase.elementblock[label] == 'f') and l not in metal_species: metal_species.append(l)

    #############################
    ### DECIDES MAGNETIZATION ###  It will be printed later
    #############################
    # CASE 1: If no instructions are given, everything is LS
    if len(spin_config) == 0: 
        magnetization = 0

    # CASE 2: If instructions are given. Note that only HS and LS is allowed. No IS possible yet
    elif len(spin_config) == len(gmol.labels): 
        magnetization = 0
        for idx, tupl in enumerate(spin_config):
            magnetization += tupl[1]
            #gmol.labels[idx] = tupl[0]
    else: print("QE_input: len(spin_config) must be 0, or equal to len(gmol.labels)")
    
    ########################### 
    ### WRITES INPUT: PART 1 ##
    ########################### 
    with open(path+name+suffix+extension, 'w') as inp:
        print(" &control", file=inp)
        print(f"    calculation='{job_type}'", file=inp)
        print("    restart_mode='from_scratch'", file=inp)
        print(f"    pseudo_dir ='{PP_Library}'", file=inp)
        print("    disk_io='low'", file=inp)
        print(f"    outdir='./WFC'", file=inp)
        #print(f"    outdir='/scratch/g4vela/QE_WFC/'", file=inp)
        print(f"    prefix='{name}{suffix}'", file=inp)

        if job_type == "opt" or job_type == "relax" or job_type == "vc-relax":
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
        print(f"    nat={natoms}, ntyp={nspecies}, ecutwfc={int(cutoff)}, ecutrho={float(cutoff)*8}", file=inp)
        if magnetization != 0: print(f"    nspin=2,", file=inp)
        else:                  print(f"    nspin=1,", file=inp)
        
        if system_type == "molecule" and hasattr(gmol, 'totcharge'):
            print(f"    tot_charge={gmol.totcharge}", file=inp)
        else:
            print("    tot_charge=0", file=inp)
        
        ## Starting Magnetization
        if len(spin_config) == 0:
            print(f"    tot_magnetization={magnetization}", file=inp)
        else:
            for idx, u in enumerate(uniques):
                if t[1] != 0: print(f"    starting_magnetization({idx+1})={t[1]}", file=inp)

        ## Hubbard: Here it is assumed that the Hubbard U term will apply to all metals
        if isHubbard:
            print("    lda_plus_u=.true.,", file=inp)
            for idx, el in enumerate(elems):
                if el in metal_species: print(f"    Hubbard_U({idx+1})={U}", file=inp)

        if isGrimme:
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
        if job_type == "opt" or job_type == "relax" or job_type == "vc-relax":
            print(" &ions", file=inp)
            print("    upscale=10", file=inp)
            print("    ion_dynamics='bfgs'", file=inp)
            print("    trust_radius_ini= 0.1D0", file=inp)
            print("/", file=inp)
        
        #///////////////////
        #// Cell control ///
        #///////////////////
        if job_type == "vc-relax":
            print(" &cell", file=inp)
            print("    cell_dynamics='bfgs'", file=inp)
            print("    cell_dofree='all'", file=inp)
            print("    cell_factor= 1.2D0", file=inp)
            print("    press= 0.D0", file=inp)
            print("    press_conv_thr= 0.5D0", file=inp)
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
        elems = list(set(gmol.labels))   ## This list of elements must be updated to include the different magnetization labels (eg. Fe1, Fe2)
        print("ATOMIC_SPECIES", file=inp)
        for idx, elem in enumerate(elems):
            if elem[-1].isdigit(): label = elem[:-1]
            else: label = elem
            weight = elemdatabase.elementweight[label]
            pp = get_pp(PP_Library, label)
            print(f"{elem} {weight:6.4f} {pp}", file=inp)
        
        #//////////////////
        #// Atom Coords ///
        #//////////////////
        print("ATOMIC_POSITIONS angstrom", file=inp)
        #    for a in gmol.atoms:
        #        print("%s  %.6f  %.6f  %.6f" % (a.label, a.coord[0], a.coord[1], a.coord[2]), file=inp)
        for idx, l in enumerate(gmol.labels):
            print(f"{l}        {gmol.pos[idx][0]:.6f}   {gmol.pos[idx][1]:.6f}   {gmol.pos[idx][2]:.6f}", file=inp)
        print("K_POINTS gamma", file=inp)

###################################################
def gen_QE_subfile(path: str, name: str, suffix: str, extension: str="", procs: int=1, queue: str="iqtc09", exe: str="pw.x", version: str="6.4.1"):
    with open(path+name+suffix+extension, 'w+') as sub:
        print(f"#!/bin/bash", file=sub)
        print(f"#$ -N {name}{suffix}", file=sub) 
        print(f"#$ -pe smp {procs}", file=sub)
        print(f"#$ -cwd", file=sub)
        print(f"#$ -o {name}{suffix}.stdout", file=sub)
        print(f"#$ -e {name}{suffix}.stderr", file=sub)
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
        print(f"mpirun -np {procs} {exe} < {name}{suffix}.inp > {name}{suffix}.out", file=sub)
    os.chmod(path+name+suffix+extension, 0o777)
        
###################################################

#def get_spin_config(cell, ):

#    #########################
#    ### IDENTIFIES METALS ###
#    #########################
#    elems = list(set(gmol.labels))
#    metal_indices = []
#    for idx, l in enumerate(gmol.labels):
#        if (elemdatabase.elementblock[l] == 'd' or elemdatabase.elementblock[l] == 'f'): metal_indices.append(idx)
#
#    if debug > 0 : print(len(metal_indices), "metals")
#    if debug > 0 : print(metal_indices)


#   CASE 3: GMOL is a molecule, and instructions are not
#  elif typ == "molecule":
#          if hasattr(mol, 'spin'):
#              spin = mol.spin
#              magnetization =  spin - 1
#          else:
#              if ((mol.eleccount + mol.totcharge) % 2 == 0):
#                  spin = 1            #Singlet
#                  magnetization = 0   #Singlet
#              else:
#                  spin = 2             #Doublet
#                  magnetization = 1    #Doublet


#class specie(object):
#    def __init__(self, label)

