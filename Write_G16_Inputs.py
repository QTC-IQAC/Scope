#!/usr/bin/env python3
import sys
import os
import numpy as np
import pickle
from collections import Counter
#from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal

from Scope.Environment import set_user, set_cluster
from Scope.Elementdata import ElementData
from Scope.Other import get_metal_idxs
from Scope.Classes_State import find_state

#######################
def gen_G16_input(comp, debug: int=0):
    gmol = comp._job._recipe.subject
    
    ## 1-Change variable names to simplify calls
    jobtype = comp.qc_data.jobtype
    functional = comp.qc_data.functional
    basis = comp.qc_data.basis
    loose_opt = comp.qc_data.loose_opt
    tight_opt = comp.qc_data.tight_opt

    assert hasattr(comp.qc_data,"istate"), f"istate = {comp.qc_data.istate} not found in gmol"
    exists, istate    = find_state(gmol, comp.qc_data.istate)
    assert exists, f"istate = {comp.qc_data.istate} does not exist" 
    assert hasattr(istate,"labels"), f"istate = {comp.qc_data.istate} doesn't have labels"
    assert hasattr(istate,"coord"),  f"istate = {comp.qc_data.istate} doesn't have coordinates"

    ### 2-Starts printing input 
    with open(comp.inp_path, 'w') as inp:
        print(f"%chk={comp.refcode}{comp.suffix}.chk", file=inp)
        
        ## 2.1-General keywords
        commandline = []
        commandline.append("#p")
        commandline.append(" nosymm")
        commandline.append(" scf=(maxconventionalcycles=200,xqc)")

        ## 2.2-Functional
        if   functional == "pbe":       commandline.append(" UPBEPBE")
        elif functional == "b3lyp":   commandline.append(" UB3LYP")
        elif functional == "b3lyp*":  commandline.append(" UB3LYP IOp(3/76=1000002000) IOp(3/77=0720008000) IOp(3/78=0810010000)")
        elif functional == "b3lyp**": commandline.append(" UB3LYP IOp(3/76=1000001500) IOp(3/77=0720008500) IOp(3/78=0810010000)")
        else: print("G16_INPUT: functional", functional, "not recognized")

        ## 2.3-Basis
        if basis == "def2sv":      commandline.append(" def2SV")
        elif basis == "def2svp":   commandline.append(" def2SVP")
        elif basis == "def2tzv":   commandline.append(" def2TZV")
        elif basis == "def2tzvp":  commandline.append(" def2TZVP")
        elif basis == "def2tzvpp": commandline.append(" def2TZVPP")
        else: print("G16_INPUT: basis", basis, "not recognized")

        ## 2.4-Jobtype
        if jobtype == "opt" or jobtype == "opth" or jobtype == "opt&freq": 
            if   loose_opt: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral,loose)")
            elif tight_opt: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral,VeryTight)")
            else:           commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral)")
        elif jobtype == "opt&freq" or jobtype == "freq": commandline.append(" freq")
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
        print(f"{gmol.totcharge} {comp.spin_config.multiplicity}", file=inp) 

        ##################################################################
        ### Coordinates, which are taken from the initial state object ###
        ##################################################################
        for idx, z in enumerate(zip(istate.labels, istate.coord)):
            if jobtype.lower() != "opth": print("%s  %.6f  %.6f  %.6f" % (z[0], z[1][0], z[1][1], z[1][2]), file=inp)
            else: 
                if a.label == 'H': print("%s %s %.6f  %.6f  %.6f" % (z[0], " 0", z[1][0], z[1][1], s[1][2]), file=inp)
                else:              print("%s $s %.6f  %.6f  %.6f" % (z[0], "-1", z[1][0], z[1][1], s[1][2]), file=inp)

        #for a in gmol.atoms:
        #    ta = getattr(a,igeom)
        #    if jobtype.lower() != "opth": print("%s  %.6f  %.6f  %.6f" % (a.label, ta[0], ta[1], ta[2]), file=inp)
        #    else: 
        #        if a.label == 'H': print("%s %s %.6f  %.6f  %.6f" % (a.label, " 0", ta[0], ta[1], ta[2]), file=inp)
        #        else:              print("%s $s %.6f  %.6f  %.6f" % (a.label, "-1", ta[0], ta[1], ta[2]), file=inp)
        print("", file=inp) 

###################################################
def gen_G16_subfile(comp: object, procs: int=1, queue: str="iqtc09", cluster: str=set_cluster(), user: str=set_user(), savechk: bool=False):
    ## savechk must be tested
   
    if    procs*float(1.5) >= float(8): mem=int(procs*1.5)
    else:                               mem=int(4)

    if 'login' in cluster or 'csuc' in cluster:
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
            print(f"module load apps/gaussian/g16c2", file=sub)
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
            if savechk: print(f"cp -pr *.chk $JOBDIR", file=sub)
            os.chmod(comp.sub_path, 0o777)

    elif 'uam' in cluster:
        project = 'ub100'
        with open(comp.sub_path, 'w+') as sub:
            print(f"#!/bin/bash", file=sub)
            print(f"#SBATCH -J {comp.refcode}{comp.suffix}", file=sub)
            print(f"#SBATCH -e /scratch/{project}/std_files/{comp.refcode}{comp.suffix}.stderr", file=sub)
            print(f"#SBATCH -o /scratch/{project}/std_files/{comp.refcode}{comp.suffix}.stdout", file=sub)
            print(f"#SBATCH -p {queue}", file=sub)
            print(f"#SBATCH --exclude=cibeles3-05", file=sub)
            print(f"#SBATCH -A ub100_serv", file=sub)
            print(f"#SBATCH --nodes=1", file=sub)
            print(f"#SBATCH --ntasks={procs}", file=sub)
            print(f"", file=sub)
            print(f"source ~/.bashrc", file=sub)
            print(f"module load gaussian/gaussian16.b1-SSE4", file=sub)
            print(f"", file=sub)
            print(f"JOBDIR=$PWD", file=sub)
            print(f'export RUNDIR="/temporal/{user}/jobs/$SLURM_JOBID"', file=sub)
            print(f'export GAUSS_SCRDIR="/temporal/{user}/jobs/$SLURM_JOBID"', file=sub)
#            print(f'export RUNDIR="/scratch/{project}/jobs/$SLURM_JOBID"', file=sub)
#            print(f'export GAUSS_SCRDIR="/scratch/{project}/jobs/$SLURM_JOBID"', file=sub)
            print(f"mkdir -p $RUNDIR", file=sub)
            print(f"cd $RUNDIR", file=sub)
            print(f"", file=sub)
            print(f"cp -i $JOBDIR/{comp.inp_name} .", file=sub)
            print(f"echo '%nprocs={procs}' >  tmp1", file=sub)
            print(f"echo '%mem={mem}gb'    >> tmp1", file=sub)
            print(f"cat tmp1 {comp.inp_name} > tmp2", file=sub)
            print(f"rm tmp1", file=sub)
            print(f"mv -f tmp2 {comp.inp_name}", file=sub)
            print(f"timeout 71h g16 < {comp.inp_name} > {comp.out_name}", file=sub)
            print(f"cp -pr *.log $JOBDIR/", file=sub)
            if savechk: print(f"cp -pr *.chk $JOBDIR", file=sub)
            print(f"cd $JOBDIR/", file=sub)
            print(f"rm $RUNDIR/Gau.*", file=sub)
            print(f"rm $RUNDIR/{comp.refcode}*", file=sub)
            os.chmod(comp.sub_path, 0o777)

    elif 'portal' in cluster:
        with open(comp.sub_path, 'w+') as sub:
            print(f"#!/bin/bash", file=sub)
            print(f"#$ -N {comp.refcode}{comp.suffix}", file=sub) 
            print(f"#$ -pe smp {procs}", file=sub)
            print(f"#$ -cwd", file=sub)
            print(f"#$ -o /home/{user}/x-stds/{comp.refcode}{comp.suffix}.stdout", file=sub)
            print(f"#$ -e /home/{user}/x-stds/{comp.refcode}{comp.suffix}.stderr", file=sub)
            print(f"#$ -q {queue}.q" , file=sub)
            print(f"", file=sub)
            print(f"source /etc/profile.d/modules.csh", file=sub)
            print(f"source $HOME/.bashrc", file=sub)
            print(f". /etc/profile", file=sub)
            print(f"module load gaussian/g16b01", file=sub)
            print(f"", file=sub)
            print(f"WORKDIR=$PWD", file=sub)
            print(f"cd $TMPDIR", file=sub)
            print(f"export GAUSS_SCRDIR $TMPDIR", file=sub)
            print(f"cp $WORKDIR/{comp.inp_name} .", file=sub)
            print(f"echo '%nprocs={procs}' >  tmp1", file=sub)
            print(f"echo '%mem={mem}gb'    >> tmp1", file=sub)
            print(f"cat tmp1 {comp.inp_name} > tmp2 | mv -f tmp2 {comp.inp_name}", file=sub)
#            print(f"sed -i '1s/^/%mem={mem}gb\n' {name}{suffix}.com", file=sub)
#            print(f"sed -i '1s/^/%nprocs={procs}\n' {name}{suffix}.com", file=sub)
            print(f"g16 < {comp.inp_name} > {comp.out_name}", file=sub)
            print(f"cp -pr *.log $WORKDIR", file=sub)
            if savechk: print(f"cp -pr *.chk $WORKDIR", file=sub)
            os.chmod(comp.sub_path, 0o777)
        
###################################################
