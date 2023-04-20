#!/usr/bin/env python3
import sys
import os
import numpy as np
import pickle
from collections import Counter
#from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal

from Test_V3.Environment import set_user, set_cluster
from Test_V3.Elementdata import ElementData
from Test_V3.Other import get_metal_idxs

#######################
def gen_G16_input(comp, debug: int=0):
    gmol = comp._job._recipe.subject
    
    ### 2-Starts printing input 
    with open(comp.inp_path, 'w') as inp:
        print(f"%chk={comp.refcode}{comp.suffix}.chk", file=inp)
        
        ## 2.1-General keywords
        commandline = []
        commandline.append("#p")
        commandline.append(" nosymm")
        commandline.append(" scf=(maxconventionalcycles=200,xqc)")

        ## 2.2-Functional
        if hasattr(comp.qc_data,"functional"): functional = comp.qc_data.functional.lower()
        else:                             functional = "B3LYP**"
        if   functional == "pbe":       commandline.append(" UPBEPBE")
        elif functional == "b3lyp":   commandline.append(" UB3LYP")
        elif functional == "b3lyp*":  commandline.append(" UB3LYP IOp(3/76=1000002000) IOp(3/77=0720008000) IOp(3/78=0810010000)")
        elif functional == "b3lyp**": commandline.append(" UB3LYP IOp(3/76=1000001500) IOp(3/77=0720008500) IOp(3/78=0810010000)")
        else: print("G16_INPUT: functional", functional, "not recognized")

        ## 2.3-Basis
        if hasattr(comp.qc_data,"basis"): basis = comp.qc_data.basis.lower()
        else:                        basis = "def2SVP"
        if basis == "def2sv":      commandline.append(" def2SV")
        elif basis == "def2svp":   commandline.append(" def2SVP")
        elif basis == "def2tzv":   commandline.append(" def2TZV")
        elif basis == "def2tzvp":  commandline.append(" def2TZVP")
        elif basis == "def2tzvpp": commandline.append(" def2TZVPP")
        else: print("G16_INPUT: basis", basis, "not recognized")

        ## 2.4-Jobtype
        if hasattr(comp.qc_data,"jobtype"): jobtype = comp.qc_data.jobtype.lower()
        else:                          jobtype = "scf"
        if hasattr(comp.qc_data,"loose_opt"): loose_opt = comp.qc_data.loose_opt
        else:                            loose_opt = False
        if hasattr(comp.qc_data,"tight_opt"): tight_opt = comp.qc_data.tight_opt
        else:                            tight_opt = False
        if jobtype == "opt" or jobtype == "opth" or jobtype == "opt&freq": 
            if   loose_opt: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral,loose)")
            elif tight_opt: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral,VeryTight)")
            else:           commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral)")
        elif jobtype == "opt&freq" or jobtype == "freq": commandline.append(" freq")
        else: print("G16_INPUT: jobtype", jobtype, "not recognized")

        ## 2.5-Grimme
        if hasattr(comp.qc_data,"isGrimme"): isGrimme = comp.qc_data.isGrimme
        else:                           isGrimme = False
        if isGrimme: commandline.append(" EmpiricalDispersion=GD3BJ")

        ## 2.6-Integral Grid
        if tight_opt: commandline.append(" Int=UltraFine")

        ## 3-Commandline is put together
        commandline = ''.join(commandline)
        print(f"{commandline}", file=inp) 
        print("", file=inp) 
        print(f"Title Card", file=inp) 
        print("", file=inp) 
        print(f"{gmol.totcharge} {comp.spin_config.multiplicity}", file=inp) 

        ####################################################
        ### Coordinates, which are taken from gmol object ###
        ####################################################
        if hasattr(comp.qc_data,"coord_tag"):  coord_tag = comp.qc_data.coord_tag.lower()
        else:                                  coord_tag = "coord"
        assert hasattr(gmol,coord_tag)
        for a in gmol.atoms:
            ta = getattr(a,coord_tag)
            if jobtype.lower() != "opth": print("%s  %.6f  %.6f  %.6f" % (a.label, ta[0], ta[1], ta[2]), file=inp)
            else: 
                if a.label == 'H': print("%s %s %.6f  %.6f  %.6f" % (a.label, " 0", ta[0], ta[1], ta[2]), file=inp)
                else:              print("%s $s %.6f  %.6f  %.6f" % (a.label, "-1", ta[0], ta[1], ta[2]), file=inp)
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
