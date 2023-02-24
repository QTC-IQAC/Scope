#!/usr/bin/env python3
import sys
import os
import numpy as np
import pickle
from collections import Counter
from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal
from cell2mol.elementdata import ElementData 
elemdatabase = ElementData()

from Scope.Control_Jobs import set_user, set_cluster

#######################
def gen_G16_iso_input(mol: object, path: str, name: str, suffix: str="", extension: str="", coord_tag: str="coord", jobtype: str="scf", functional: str="B3LYP*", basis: str='def2SVP', spin: str='gmol', isGrimme: bool=True, loose_opt: bool=False, tight_opt: bool=False, useconnec: bool=False, append: bool=False):
## useconnec is not implemented. It is meant to call the generation of the connectivity section for G16
## append is not implemented. Will be used to append a computation to an existing input file

    assert hasattr(mol,coord_tag)  ## Asserts that the coordinates exist

    ### Tries to find spin
    if spin.lower() == "gmol": 
        if hasattr(mol, 'spin'):
            spinval = mol.spin
        else:
            if ((mol.eleccount + mol.totcharge) % 2 == 0): spinval = 1         #Singlet
            else:  spinval = 2                                                 #Doublet
    elif spin.lower() == "hs" or spin.lower() == "ls":
        spinval = 0
        for met in mol.metalist:
            ## I Should add Oxidation-State-dependent Rules
            if met.label == "Fe" and spin.lower() == "ls": spinval = spinval + 1
            elif met.label == "Fe" and spin.lower() == "hs": spinval = spinval + 5
    
    ### Starts printing input 
    with open(path+name+suffix+extension, 'w') as inp:
        #print(f"%nprocs={nproc}", file=inp)
        #if nproc*float(1.5) >= float(8): print(f"%mem={int(nproc*1.5)}gb", file=inp)
        #else: print(f"%mem=4gb", file=inp)
        print(f"%chk={name+suffix}.chk", file=inp)
        
        ## General keywords
        commandline = []
        commandline.append("#p")
        commandline.append(" nosymm")
        commandline.append(" scf=(maxconventionalcycles=200,xqc)")

        ## Functional
        if functional.lower() == "pbe": commandline.append(" UPBEPBE")
        elif functional.lower() == "b3lyp": commandline.append(" UB3LYP")
        elif functional.lower() == "b3lyp*": commandline.append(" UB3LYP IOp(3/76=1000002000) IOp(3/77=0720008000) IOp(3/78=0810010000)")
        elif functional.lower() == "b3lyp**": commandline.append(" UB3LYP IOp(3/76=1000001500) IOp(3/77=0720008500) IOp(3/78=0810010000)")

        ## Basis
        if basis.lower() == "def2sv": commandline.append(" def2SV")
        elif basis.lower() == "def2svp": commandline.append(" def2SVP")
        elif basis.lower() == "def2tzv": commandline.append(" def2TZV")
        elif basis.lower() == "def2tzvp": commandline.append(" def2TZVP")
        elif basis.lower() == "def2tzvpp": commandline.append(" def2TZVPP")

        ## Jobtype
        if jobtype.lower() == "opt" or jobtype.lower() == "opth" or jobtype.lower() == "opt&freq": 
            if   loose_opt: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral,loose)")
            elif tight_opt: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral,VeryTight)")
            else: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral)")
        if jobtype.lower() == "opt&freq" or jobtype.lower() == "freq": commandline.append(" freq")
        if isGrimme: commandline.append(" EmpiricalDispersion=GD3BJ")
        if tight_opt: commandline.append(" Int=UltraFine")

        commandline = ''.join(commandline)
        print(f"{commandline}", file=inp) 
        print("", file=inp) 
        print(f"Title Card", file=inp) 
        print("", file=inp) 
        print(f"{mol.totcharge} {spinval}", file=inp) 

        ### Retrieves the desired coordinates from mol object
        for a in mol.atoms:
            ta = getattr(a,coord_tag)
            if jobtype.lower() != "opth": print("%s  %12.6f  %12.6f  %12.6f" % (a.label, ta[0], ta[1], ta[2]), file=inp)
            else: 
                if a.label == 'H': print("%s %s %12.6f  %12.6f  %12.6f" % (a.label, " 0", ta[0], ta[1], ta[2]), file=inp)
                else:              print("%s $s %12.6f  %12.6f  %12.6f" % (a.label, "-1", ta[0], ta[1], ta[2]), file=inp)
        print("", file=inp) 

###################################################
def gen_G16_subfile(path, name, suffix, extension="", procs: int=1, queue: str="iqtc09", cluster: str=set_cluster(), user: str=set_user(), savechk: bool=False):
    ## savechk must be tested
   
    if    procs*float(1.5) >= float(8): mem=int(procs*1.5)
    else:                               mem=int(4)

    if 'login' in cluster or 'csuc' in cluster:
        with open(path+name+suffix+extension, 'w+') as sub:
            print(f"#!/bin/bash", file=sub)
            print(f"#SBATCH -J {name}{suffix}", file=sub)
            print(f"#SBATCH -e /scratch/{user}/std_files/{name}{suffix}.err", file=sub)
            print(f"#SBATCH -o /scratch/{user}/std_files/{name}{suffix}.out", file=sub)
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
            print(f"cp $JOBDIR/{name}{suffix}.com .", file=sub)
            print(f"echo '%nprocs={procs}' >  tmp1", file=sub)
            print(f"echo '%mem={mem}gb'    >> tmp1", file=sub)
            print(f"cat tmp1 {name}{suffix}.com > tmp2 | mv -f tmp2 {name}{suffix}.com", file=sub)
            #print(f"sed -i '1s/^/%mem={mem}gb\n' {name}{suffix}.com", file=sub)
            #print(f"sed -i '1s/^/%nprocs={procs}\n' {name}{suffix}.com", file=sub)
            print(f"g16 < {name}{suffix}.com > {name}{suffix}.log", file=sub)
            print(f"cp -pr *.log $JOBDIR/", file=sub)
            if savechk: print(f"cp -pr *.chk $JOBDIR", file=sub)
            os.chmod(path+name+suffix+extension, 0o777)

    elif 'portal' in cluster:
        with open(path+name+suffix+extension, 'w+') as sub:
            print(f"#!/bin/bash", file=sub)
            print(f"#$ -N {name}{suffix}", file=sub) 
            print(f"#$ -pe smp {procs}", file=sub)
            print(f"#$ -cwd", file=sub)
            print(f"#$ -o /home/{user}/x-stds/{name}{suffix}.stdout", file=sub)
            print(f"#$ -e /home/{user}/x-stds/{name}{suffix}.stderr", file=sub)
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
            print(f"cp $WORKDIR/{name}{suffix}.com .", file=sub)
            print(f"echo '%nprocs={procs}' >  tmp1", file=sub)
            print(f"echo '%mem={mem}gb'    >> tmp1", file=sub)
            print(f"cat tmp1 {name}{suffix}.com > tmp2 | mv -f tmp2 {name}{suffix}.com", file=sub)
#            print(f"sed -i '1s/^/%mem={mem}gb\n' {name}{suffix}.com", file=sub)
#            print(f"sed -i '1s/^/%nprocs={procs}\n' {name}{suffix}.com", file=sub)
            print(f"g16 < {name}{suffix}.com > {name}{suffix}.log", file=sub)
            print(f"cp -pr *.log $WORKDIR", file=sub)
            if savechk: print(f"cp -pr *.chk $WORKDIR", file=sub)
            os.chmod(path+name+suffix+extension, 0o777)
        
###################################################
def g16_from_inpt(inpt: object, mol: object, path: str, name: str, suffix: str="", extension: str="", spin: str='gmol'):

    ### Tries to find spin
    if spin.lower() == "gmol": 
        if hasattr(mol, 'spin'):
            spinval = mol.spin
        else:
            if ((mol.eleccount + mol.totcharge) % 2 == 0): spinval = 1         #Singlet
            else:  spinval = 2                                                 #Doublet
    elif spin.lower() == "hs" or spin.lower() == "ls":
        spinval = 0
        for met in mol.metalist:
            ## I Should add Oxidation-State-dependent Rules
            if met.label == "Fe" and spin.lower() == "ls": spinval = spinval + 1
            elif met.label == "Fe" and spin.lower() == "hs": spinval = spinval + 5
    
    ### Starts printing input 
    with open(path+name+suffix+extension, 'w') as inp:
        print(f"%chk={name+suffix}.chk", file=inp)
        
        ## General keywords
        commandline = []
        commandline.append("#p")
        commandline.append(" nosymm")
        commandline.append(" scf=(maxconventionalcycles=200,xqc)")

        ## Functional
        if "functional" in inpt.keys(): functional = inpt['functional'].lower()
        else:                           functional = "B3LYP**"
        if functional == "pbe": commandline.append(" UPBEPBE")
        elif functional == "b3lyp": commandline.append(" UB3LYP")
        elif functional == "b3lyp*": commandline.append(" UB3LYP IOp(3/76=1000002000) IOp(3/77=0720008000) IOp(3/78=0810010000)")
        elif functional == "b3lyp**": commandline.append(" UB3LYP IOp(3/76=1000001500) IOp(3/77=0720008500) IOp(3/78=0810010000)")
        else: print("G16_INPUT: functional", functional, "not recognized")

        ## Basis
        if "basis" in inpt.keys(): basis = inpt['basis'].lower()
        else:                      basis = "def2SVP"
        if basis.lower() == "def2sv": commandline.append(" def2SV")
        elif basis.lower() == "def2svp": commandline.append(" def2SVP")
        elif basis.lower() == "def2tzv": commandline.append(" def2TZV")
        elif basis.lower() == "def2tzvp": commandline.append(" def2TZVP")
        elif basis.lower() == "def2tzvpp": commandline.append(" def2TZVPP")
        else: print("G16_INPUT: basis", basis, "not recognized")

        ## Jobtype
        if "jobtype" in inpt.keys():   jobtype = inpt['jobtype'].lower()
        else:                          jobtype = "scf"
        if "loose_opt" in inpt.keys(): loose_opt = inpt['loose_opt']
        else:                          loose_opt = False
        if "tight_opt" in inpt.keys(): tight_opt = inpt['tight_opt']
        else:                          tight_opt = False
        if jobtype.lower() == "opt" or jobtype.lower() == "opth" or jobtype.lower() == "opt&freq": 
            if   loose_opt: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral,loose)")
            elif tight_opt: commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral,VeryTight)")
            else:           commandline.append(" opt=(RecalcFC=30,cartesian,skipdihedral)")
        elif jobtype.lower() == "opt&freq" or jobtype.lower() == "freq": commandline.append(" freq")
        else: print("G16_INPUT: jobtype", jobtype, "not recognized")

        ## Grimme
        if "isGrimme" in inpt.keys(): isGrimme = inpt['isGrimme']
        else:                         isGrimme = False
        if isGrimme: commandline.append(" EmpiricalDispersion=GD3BJ")

        ## Integral Grid
        if tight_opt: commandline.append(" Int=UltraFine")

        ## Commandline is put together
        commandline = ''.join(commandline)
        print(f"{commandline}", file=inp) 
        print("", file=inp) 
        print(f"Title Card", file=inp) 
        print("", file=inp) 
        print(f"{mol.totcharge} {spinval}", file=inp) 

        ####################################################
        ### Coordinates, which are taken from mol object ###
        ####################################################
        if "coord_tag" in inpt.keys(): coord_tag = inpt['coord_tag'].lower()
        else:                          coord_tag = "coord"
        assert hasattr(mol,coord_tag)
        for a in mol.atoms:
            ta = getattr(a,coord_tag)
            if jobtype.lower() != "opth": print("%s  %.6f  %.6f  %.6f" % (a.label, ta[0], ta[1], ta[2]), file=inp)
            else: 
                if a.label == 'H': print("%s %s %.6f  %.6f  %.6f" % (a.label, " 0", ta[0], ta[1], ta[2]), file=inp)
                else:              print("%s $s %.6f  %.6f  %.6f" % (a.label, "-1", ta[0], ta[1], ta[2]), file=inp)
        print("", file=inp) 

###################################################

def identify_metals(labels: list):
    elems = list(set(labels))
    thereismetal = []
    for idx, l in enumerate(elems):
        if (elemdatabase.elementblock[l] == 'd' or elemdatabase.elementblock[l] == 'f'): thereismetal.append(True)
        else: thereismetal.append(False)
    if any(thereismetal):
        nummetal = sum(thereismetal)
        trues = [i for i, x in enumerate(thereismetal) if x]
    return nummetal, trues
