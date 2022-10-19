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
def gen_QE_iso_input(mol, path, name, suffix="", extension="", PP_Library="", typ="scf", isHubbard=False, isGrimme=True, U=2.27, cutoff=70, cubeside=70.0):
 
    #IDENTIFIES METALS
    elems = list(set(mol.labels)) 
    thereismetal = []
    for idx, l in enumerate(elems): 
        if (elemdatabase.elementblock[l] == 'd' or elemdatabase.elementblock[l] == 'f'): thereismetal.append(True)
        else: thereismetal.append(False)
    if any(thereismetal):
        nummetal = sum(thereismetal)
        trues = [i for i, x in enumerate(thereismetal) if x]
        
    # IS THERE SPIN?
    if hasattr(mol, 'spin'):
        spin = mol.spin
        magnetization = spin - 1
    else:
        if ((mol.eleccount + mol.totcharge) % 2 == 0):
            spin = 1         #Singlet
            magnetization = 0   #Singlet 
        else:
            spin = 2         #Doublet
            magnetization = 1   #Doublet
    
    with open(path+name+suffix+extension, 'w') as inp:
        print(" &control", file=inp)
        print(f"    calculation='{typ}'", file=inp)
        print("    restart_mode='from_scratch'", file=inp)
        print(f"    pseudo_dir ='{PP_Library}'", file=inp)
        print("    disk_io='low'", file=inp)
        print(f"    outdir='./WFC'", file=inp)
        #print(f"    outdir='/scratch/g4vela/QE_WFC/'", file=inp)
        print(f"    prefix='{name}{suffix}'", file=inp)

        if typ == "opt" or typ == "relax" or typ == "vc-relax":
            print(f"    wf_collect = .true.", file=inp)
            print(f"    tprnfor = .true.", file=inp)
            print(f"    nstep = 300", file=inp)
            print(f"    forc_conv_thr = 1.0D-5", file=inp)

        print("/", file=inp)

        ######################
        ### System control ###
        ######################
        print(" &system", file=inp)
        print(f"    ibrav=1, celldm(1)={cubeside}", file=inp)
        print(f"    nat={mol.natoms}, ntyp={len(list(set(mol.labels)))}, ecutwfc={cutoff}, ecutrho={cutoff*8}", file=inp)
        print("    nspin=2,", file=inp)
        
        if hasattr(mol, 'totcharge'):
            print(f"    tot_charge={mol.totcharge}", file=inp)
        else:
            print("    tot_charge=0", file=inp)
        
        if any(thereismetal): 
            totalmagnetization = 0
            for idx, met in enumerate(thereismetal):
                if met: print(f"    starting_magnetization({idx+1})={magnetization}", file=inp)
                if met: totalmagnetization += magnetization
            else:
                pass
        else: 
            totalmagnetization = 0
        print(f"    tot_magnetization={totalmagnetization}", file=inp)
        
        if isHubbard:
            if any(thereismetal): 
                print("    lda_plus_u=.true.,", file=inp)
                for idx, met in enumerate(thereismetal):
                    if met: print(f"    Hubbard_U({idx+1})={U}", file=inp)

        if isGrimme:
            print("    vdw_corr='grimme-d3'", file=inp)
            print("    dftd3_version=4", file=inp)
            
        print("    assume_isolated='mp'", file=inp)
        print("/", file=inp)
        
        #########################
        ### Electrons control ###
        #########################
        print(" &electrons", file=inp)
        print("    diagonalization='david'", file=inp)
        print("    electron_maxstep=150", file=inp)
        print("    conv_thr = 1.0e-5", file=inp)
        print("    mixing_beta = 0.20", file=inp)
        print("/", file=inp)
   
        ### Ions control ###
        if typ == "opt" or typ == "relax" or typ == "vc-relax":
            print(" &ions", file=inp)
            print("    upscale=10", file=inp)
            print("    ion_dynamics='bfgs'", file=inp)
            print("    trust_radius_ini= 0.1D0", file=inp)
            print("/", file=inp)
        
        ### Cell control ###
        if typ == "vc-relax":
            print(" &cell", file=inp)
            print("    cell_dynamics='bfgs'", file=inp)
            print("    cell_dofree='all'", file=inp)
            print("    cell_factor= 1.2D0", file=inp)
            print("    press= 0.D0", file=inp)
            print("    press_conv_thr= 0.5D0", file=inp)
            print("/", file=inp)

        print("ATOMIC_SPECIES", file=inp)
        for idx, elem in enumerate(elems):
            label = elem
            weight = elemdatabase.elementweight[elem]
            pp = get_pp(PP_Library, elem)
            print(label, weight, pp, file=inp)
        
        print("ATOMIC_POSITIONS angstrom", file=inp)
        for a in mol.atoms:
            print("%s  %.6f  %.6f  %.6f" % (a.label, a.coord[0], a.coord[1], a.coord[2]), file=inp)
        print("K_POINTS gamma", file=inp)

###################################################
def gen_QE_subfile(path, name, suffix, extension="", cores=1, queue="iqtc09", exe="pw.x"):
    with open(path+name+suffix+extension, 'w+') as sub:
        print(f"#!/bin/bash", file=sub)
        print(f"#$ -N {name}{suffix}", file=sub) 
        print(f"#$ -pe smp {cores}", file=sub)
        print(f"#$ -cwd", file=sub)
        print(f"#$ -o {name}{suffix}.stdout", file=sub)
        print(f"#$ -e {name}{suffix}.stderr", file=sub)
        print(f"#$ -S /bin/bash", file=sub)
        print(f"#$ -q {queue}.q" , file=sub)
        print(f"", file=sub)
        print(f"source /etc/profile.d/modules.csh", file=sub)
        print(f"module load espresso/6.4.1_ompi", file=sub)
        print(f"", file=sub)
        print(f"set OMP_NUM_THREADS=1", file=sub)
        print(f"ulimit -l unlimited", file=sub)
        print(f"", file=sub)
        print(f"mpirun -np {cores} {exe} < {name}{suffix}.inp > {name}{suffix}.out", file=sub)
    os.chmod(path+name+suffix+extension, 0o777)
        
###################################################
