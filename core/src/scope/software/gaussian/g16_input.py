#!/usr/bin/env python3
import sys
import os
import numpy as np
import pickle
from collections import Counter

from scope.other               import get_metal_idxs
from scope.classes_state       import find_state
from scope.classes_environment import set_user
from scope.elementdata         import ElementData
elemdatabase = ElementData()

#######################
def gen_g16_input(comp, debug: int=0):

    ## 1-Change some variable names to simplify calls
    source     = comp.source
    comp_type  = comp.qc_data.comp_type
    functional = comp.qc_data.functional
    basis      = comp.qc_data.basis
    loose_opt  = comp.qc_data.loose_opt
    tight_opt  = comp.qc_data.tight_opt

    assert hasattr(comp.qc_data,"istate"), f"istate = {comp.qc_data.istate} not found in comp.qc_data"
    exists, istate    = find_state(source, comp.qc_data.istate)
    assert exists,                   f"istate = {comp.qc_data.istate} does not exist in source" 
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

        ## 2.2 Decides Unrestricted or Restricted based on the spin multiplicity of the initial state
        is_unrestricted = False
        if source.spin > 0: 
            is_unrestricted = True
            commandline.append(" U")
        else:
            commandline.append(" ")

        ## 2.2-Functional
        if   functional == "pbe":     commandline.append("PBEPBE")
        elif functional == "b3lyp":   commandline.append("B3LYP")
        elif functional == "b3lyp*":  commandline.append("B3LYP IOp(3/76=1000002000) IOp(3/77=0720008000) IOp(3/78=0810010000)")
        elif functional == "b3lyp**": commandline.append("B3LYP IOp(3/76=1000001500) IOp(3/77=0720008500) IOp(3/78=0810010000)")
        elif functional == "wb97xd":  commandline.append("wB97XD")
        else: 
            raise ValueError("G16_INPUT: functional", functional, "not recognized. Implemented are pbe, wb97xd, b3lyp, b3lyp* and b3lyp**")

        ## 2.3-Basis
        basis = basis.lower()
        if basis in BASIS_ALIASES:
            commandline.append(f" {BASIS_ALIASES[basis]}")
        else: 
            raise ValueError("G16_INPUT: basis", basis, f"not recognized. Implemented are {SUPPORTED_BASIS_NAMES}")

        ## 2.4-Jobtype
        if comp_type == "opt" or comp_type == "opth" or comp_type == "opt&freq": 
            if comp.qc_data.fctype == 'recalcfc':
                if   loose_opt: commandline.append(f" opt=(RecalcFC={comp.qc_data.recalc_steps},cartesian,skipdihedral,loose)")
                elif tight_opt: commandline.append(f" opt=(RecalcFC={comp.qc_data.recalc_steps},cartesian,skipdihedral,VeryTight)")
                else:           commandline.append(f" opt=(RecalcFC={comp.qc_data.recalc_steps},cartesian,skipdihedral)")
            elif comp.qc_data.fctype == 'calcall':
                if   loose_opt: commandline.append(f" opt=(CalcAll,cartesian,skipdihedral,loose)")
                elif tight_opt: commandline.append(f" opt=(CalcAll,cartesian,skipdihedral,VeryTight)")
                else:           commandline.append(f" opt=(CalcAll,cartesian,skipdihedral)")
            elif comp.qc_data.fctype == 'calcfc':
                if   loose_opt: commandline.append(f" opt=(CalcFC,cartesian,skipdihedral,loose)")
                elif tight_opt: commandline.append(f" opt=(CalcFC,cartesian,skipdihedral,VeryTight)")
                else:           commandline.append(f" opt=(CalcFC,cartesian,skipdihedral)")
        
        elif comp_type == "ts":
            if comp.qc_data.fctype == 'recalcfc':
                commandline.append(f" opt=(TS,NoEigenTest,RecalcFC={comp.qc_data.recalc_steps},Cartesian,SkipDihedral,MaxCycles=300)")
            elif comp.qc_data.fctype == 'calcall':
                commandline.append(f" opt=(TS,NoEigenTest,CalcAll,Cartesian,SkipDihedral,MaxCycles=300)")
            elif comp.qc_data.fctype == 'calcfc':
                commandline.append(f" opt=(TS,NoEigenTest,CalcFC,Cartesian,SkipDihedral,MaxCycles=300)")

        elif comp_type == "opt&freq" or comp_type == "ts&freq" or comp_type == "freq": 
            commandline.append(" freq(hpmodes)")

        elif comp_type == "td":  
            if not is_unrestricted:  commandline.append(f" td=({comp.qc_data.td_type},nstates={comp.qc_data.td_nstates})")
            else:                    commandline.append(f" td=(nstates={comp.qc_data.td_nstates})")
        elif comp_type == "tda": 
            if not is_unrestricted:  commandline.append(f" tda=({comp.qc_data.td_type},nstates={comp.qc_data.td_nstates})")
            else:                    commandline.append(f" tda=(nstates={comp.qc_data.td_nstates})")

        elif comp_type == "scf": pass
        else: 
            raise ValueError("G16_INPUT: comp_type", comp_type, "not recognized. Implemented are opt, opt&freq, freq, td, tda and scf")

        ## 2.5-Grimme
        if functional == "wb97xd": comp.qc_data._mod_attr("is_grimme", False)       ## Grimme Disperson Corrections are already implemented in wb97xd
        if comp.qc_data.is_grimme: 
            if   comp.qc_data.grimme_type == 'd3bj': commandline.append(" EmpiricalDispersion=GD3BJ")
            elif comp.qc_data.grimme_type == 'd3':   commandline.append(" EmpiricalDispersion=GD3")
            elif comp.qc_data.grimme_type == 'd2':   commandline.append(" EmpiricalDispersion=GD2")

        ## 2.6-Integral Grid
        if tight_opt: commandline.append(" Int=UltraFine")

        ## 2.7-TD-DFT options for Theodore: 
        if comp_type == "td" or comp_type == "tda": commandline.append(" pop=full iop(9/40=3) GFINPUT Integral=NoXCTest")

        ## 2.8-Other options
        if comp.qc_data.ultrafine_grid: commandline.append(" Int=Ultrafine")

        ## 3-Commandline is put together and file is written
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
            if comp_type.lower() != "opth": print("%s  %.6f  %.6f  %.6f" % (z[0], z[1][0], z[1][1], z[1][2]), file=inp)
            else: 
                if a.label == 'H': print("%s %s %.6f  %.6f  %.6f" % (z[0], " 0", z[1][0], z[1][1], s[1][2]), file=inp)
                else:              print("%s $s %.6f  %.6f  %.6f" % (z[0], "-1", z[1][0], z[1][1], s[1][2]), file=inp)
        print("", file=inp) 

###################################################
def gen_g16_subfile(comp: object, queue: object, module: str, procs: int=1, savechk: bool=False):
   
    if    procs*float(1.5) >= float(8): mem=int(procs*1.5)   ## Estimates the available memory as 1.5 GB per CPU
    else:                               mem=int(4)           ## With a minimum of 4GB

    with open(comp.sub_path, 'w+') as sub:
        if queue._environment.scheduler == 'slurm':
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
            print(f"JOBDIR=$PWD", file=sub)
            print(f"cd $TMPDIR", file=sub)
            print(f"export GAUSS_SCRDIR $TMPDIR", file=sub)
            print(f"cp $JOBDIR/{comp.inp_name} .", file=sub)
            print(f"echo '%nprocs={procs}' >  tmp1", file=sub)
            print(f"echo '%mem={mem}gb'    >> tmp1", file=sub)
            print(f"cat tmp1 {comp.inp_name} > tmp2 | mv -f tmp2 {comp.inp_name}", file=sub)
            print(f"g16 < {comp.inp_name} > {comp.out_name}", file=sub)
            print(f"cp -pr *.log $JOBDIR", file=sub)
            if savechk: print(f"cp -pr *.chk $JOBDIR", file=sub)

        os.chmod(comp.sub_path, 0o777)

################
## BASIS SETS ##
################
BASIS_ALIASES = {
    "def2sv": "def2SV",
    "def2svp": "def2SVP",
    "def2tzv": "def2TZV",
    "def2tzvp": "def2TZVP",
    "def2tzvpp": "def2TZVPP",
    "sto-3g": "STO-3G",
    "sto-3g*": "STO-3G*",
    "3-21g": "3-21G",
    "3-21g*": "3-21G*",
    "6-21g": "6-21G",
    "6-21g*": "6-21G*",
    "6-21g(d)": "6-21G*",
    "6-21g**": "6-21G**",
    "6-21g(d,p)": "6-21G(d,p)",
    "4-31g": "4-31G",
    "4-31g*": "4-31G*",
    "4-31g(d)": "4-31G*",
    "4-31g**": "4-31G**",
    "4-31g(d,p)": "4-31G(d,p)",
    "6-31g": "6-31G",
    "6-31g*": "6-31G(d)",
    "6-31g(d)": "6-31G(d)",
    "6-31g**": "6-31G(d,p)",
    "6-31g(d,p)": "6-31G(d,p)",
    "6-31+g": "6-31+G",
    "6-31+g*": "6-31+G(d)",
    "6-31+g(d)": "6-31+G(d)",
    "6-31+g**": "6-31+G(d,p)",
    "6-31+g(d,p)": "6-31+G(d,p)",
    "6-31++g": "6-31++G",
    "6-31++g*": "6-31++G(d)",
    "6-31++g(d)": "6-31++G(d)",
    "6-31++g**": "6-31++G(d,p)",
    "6-31++g(d,p)": "6-31++G(d,p)",
    "6-311g": "6-311G",
    "6-311g*": "6-311G(d)",
    "6-311g(d)": "6-311G(d)",
    "6-311g**": "6-311G(d,p)",
    "6-311g(d,p)": "6-311G(d,p)",
    "6-311+g": "6-311+G",
    "6-311+g*": "6-311+G(d)",
    "6-311+g(d)": "6-311+G(d)",
    "6-311+g**": "6-311+G(d,p)",
    "6-311+g(d,p)": "6-311+G(d,p)",
    "6-311++g": "6-311++G",
    "6-311++g*": "6-311++G(d)",
    "6-311++g(d)": "6-311++G(d)",
    "6-311++g**": "6-311++G(d,p)",
    "6-311++g(d,p)": "6-311++G(d,p)",
    "cc-pvdz": "cc-pVDZ",
    "cc-pvtz": "cc-pVTZ",
    "cc-pvqz": "cc-pVQZ",
    "cc-pv5z": "cc-pV5Z",
    "cc-pv6z": "cc-pV6Z",
    "aug-cc-pvdz": "aug-cc-pVDZ",
    "aug-cc-pvtz": "aug-cc-pVTZ",
    "aug-cc-pvqz": "aug-cc-pVQZ",
    "aug-cc-pv5z": "aug-cc-pV5Z",
    "aug-cc-pv6z": "aug-cc-pV6Z",
}

SUPPORTED_BASIS_NAMES = ", ".join(BASIS_ALIASES.keys())