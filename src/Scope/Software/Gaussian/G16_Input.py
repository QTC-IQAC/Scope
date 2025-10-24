#!/usr/bin/env python3
import sys
import os
import numpy as np
import pickle
from collections import Counter

from scope.other               import get_metal_idxs
from scope.classes_State       import find_state
from scope.classes_Environment import set_user
from scope.elementdata         import ElementData
elemdatabase = ElementData()

#######################
def gen_g16_input(comp, debug: int=0):
    ## 1-Change some variable names to simplify calls
    source     = comp.source
    jobtype    = comp.qc_data.jobtype
    functional = comp.qc_data.functional
    basis      = comp.qc_data.basis
    loose_opt  = comp.qc_data.loose_opt
    tight_opt  = comp.qc_data.tight_opt

    assert hasattr(comp.qc_data,"istate"), f"istate = {comp.qc_data.istate} not found in comp.qc_data"
    exists, istate    = find_state(source, comp.qc_data.istate)
    assert exists, f"istate = {comp.qc_data.istate} does not exist in source" 
    assert hasattr(istate,"labels"), f"istate = {comp.qc_data.istate} doesn't have labels"
    assert hasattr(istate,"coord"),  f"istate = {comp.qc_data.istate} doesn't have coordinates"

    ### 2-Starts printing input 
    with open(comp.inp_path, 'w') as inp:
        print(f"%chk={comp.name}.chk", file=inp)
        
        ## 2.1-General keywords
        commandline = []
        commandline.append("#p")
        commandline.append(" nosymm")
        commandline.append(" scf=(maxconventionalcycles=200,xqc)")

        ## 2.2-Functional
        if   functional == "pbe":     commandline.append(" UPBEPBE")
        elif functional == "b3lyp":   commandline.append(" UB3LYP")
        elif functional == "b3lyp*":  commandline.append(" UB3LYP IOp(3/76=1000002000) IOp(3/77=0720008000) IOp(3/78=0810010000)")
        elif functional == "b3lyp**": commandline.append(" UB3LYP IOp(3/76=1000001500) IOp(3/77=0720008500) IOp(3/78=0810010000)")
        else: print("G16_INPUT: functional", functional, "not recognized")

        ## 2.3-Basis
        if   basis == "def2sv":    commandline.append(" def2SV")
        elif basis == "def2svp":   commandline.append(" def2SVP")
        elif basis == "def2tzv":   commandline.append(" def2TZV")
        elif basis == "def2tzvp":  commandline.append(" def2TZVP")
        elif basis == "def2tzvpp": commandline.append(" def2TZVPP")
        elif basis == "sto-3g":    commandline.append(" STO-3G")
        else: print("G16_INPUT: basis", basis, "not recognized")

        ## 2.4-Jobtype
        if jobtype == "opt" or jobtype == "opth" or jobtype == "opt&freq": 
            if   loose_opt: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral,loose)")
            elif tight_opt: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral,VeryTight)")
            else:           commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral)")
        elif jobtype == "opt&freq" or jobtype == "freq": commandline.append(" freq(hpmodes)")
        else: print("G16_INPUT: jobtype", jobtype, "not recognized")

        ## 2.5-Grimme
        if comp.qc_data.is_grimme: commandline.append(" EmpiricalDispersion=GD3BJ")

        ## 2.6-Integral Grid
        if tight_opt: commandline.append(" Int=UltraFine")

        ## 3-Commandline is put together
        commandline = ''.join(commandline)
        print(f"{commandline}", file=inp) 
        print("", file=inp) 
        print(f"Title Card", file=inp) 
        print("", file=inp) 
        print(f"{istate.charge} {istate.spin_multiplicity}", file=inp) 

        ##################################################################
        ### Coordinates, which are taken from the initial state object ###
        ##################################################################
        for idx, z in enumerate(zip(istate.labels, istate.coord)):
            if jobtype.lower() != "opth": print("%s  %.6f  %.6f  %.6f" % (z[0], z[1][0], z[1][1], z[1][2]), file=inp)
            else: 
                if a.label == 'H': print("%s %s %.6f  %.6f  %.6f" % (z[0], " 0", z[1][0], z[1][1], s[1][2]), file=inp)
                else:              print("%s $s %.6f  %.6f  %.6f" % (z[0], "-1", z[1][0], z[1][1], s[1][2]), file=inp)
        print("", file=inp) 

###################################################
def gen_g16_subfile(comp: object, queue: object, module: str, procs: int=1, savechk: bool=False):
   
    if    procs*float(1.5) >= float(8): mem=int(procs*1.5)   ## Estimates the available memory as 1.5 GB per CPU
    else:                               mem=int(4)           ## With a minimum of 4GB

    with open(comp.sub_path, 'w+') as sub:
        print(f"#!/bin/bash", file=sub)
        print(f"#SBATCH -J {comp.name}", file=sub)
        print(f"#SBATCH -e {comp.name}.stderr", file=sub)
        print(f"#SBATCH -o {comp.name}.stdout", file=sub)
        print(f"#SBATCH -p {queue.name}", file=sub)
        print(f"#SBATCH --nodes=1", file=sub)
        print(f"#SBATCH --ntasks={procs}", file=sub)
        print(f"#SBATCH --time={queue.time_limit}", file=sub)
        print(f"", file=sub)
        print(f"module load {module}", file=sub)
        print(f"", file=sub)
        print(f"JOBDIR=$PWD", file=sub)
        print(f"cd $TMPDIR", file=sub)
        print(f"export GAUSS_SCRDIR $TMPDIR", file=sub)
        print(f"cp $JOBDIR/{comp.inp_name} .", file=sub)
        print(f"echo '%nprocs={procs}' >  tmp1", file=sub)
        print(f"echo '%mem={mem}gb'    >> tmp1", file=sub)
        print(f"cat tmp1 {comp.inp_name} > tmp2 | mv -f tmp2 {comp.inp_name}", file=sub)
        print(f"g16 < {comp.inp_name} > {comp.out_name}", file=sub)
        print(f"cp -pr *.log $JOBDIR/", file=sub)
        if savechk: print(f"cp -pr *.chk $JOBDIR", file=sub)  ## This option must be tested
        os.chmod(comp.sub_path, 0o777)                       

###################################################
