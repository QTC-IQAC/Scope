#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sklearn as sk
import ase
import ase.io as aio
import qml
import time
import gc
import scipy as sp
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split
from sklearn.metrics.pairwise import pairwise_distances
#import ase.calculators.dftb as adftb
import sklearn.manifold 
import sklearn.manifold

print("Sklearn version:", sklearn.__version__)

###########################################
pwd = os.getcwd()
pwd = pwd.replace("\\","/")

#moldir = "/home/velallau/Postdoc/Marvel_TM_Database/Scripts/FPSampling/Test/"
moldir = pwd+"/3_reps/"
#moldir = "/home/velallau/Postdoc/Marvel_TM_Database/7-Nickel/3-Inputs/1-Mono/2-V15d/0-All/Complexes/XYZs/"
metals = ["Fe", "Mn", "Ru", "Cu", "Ni", "Re", "Co", "Cr"]

xyz_files = []
names = []
for f in sorted(os.listdir(moldir)):
    if f.endswith('.xyz'):
        xyz_files.append(f)
        splitname = f.split(".")
        corename = ''.join(str(i) for i in splitname[:-1])
        names.append(corename)

mols = [aio.read(moldir + mol_file) for mol_file in xyz_files]
#nuclear_charges = mols[0].get_atomic_numbers()
nsys=len(xyz_files)

natoms = []
atnums = []
for f in mols:
    mol_names = list(f.info.keys())#[0]
    atnums.append(f.get_atomic_numbers())
    natoms.append(len(f.get_atomic_numbers()))
    sizes = len(f.get_chemical_symbols())

mbtypes = qml.representations.get_slatm_mbtypes(atnums)
maxsize = np.max(natoms)

#np.save("mbtypes",mbtypes)
mbtypes = np.load("mbtypes.npy",allow_pickle=True)
grid_density=0.4 #was 1.5
cutoff=4.0       #was 3.5

#slatms = []
localslatms = []
#cms = []
for i, mol in enumerate(mols):
    try:
        np.load(names[i]+"_local.npy")
        print("loaded for", names[i])
    except:
        local = qml.representations.generate_slatm(mol.positions, atnums[i], mbtypes=mbtypes, local=True, dgrids=[grid_density,grid_density], rcut=cutoff)
        #full = qml.representations.generate_slatm(mol.positions, atnums[i], mbtypes=mbtypes, local=False, dgrids=[grid_density,grid_density], rcut=cutoff)
        if i == 0 : print("length of slatm:", len(local[0]))
     
        metal_index = np.where([x in metals for x in mol.get_chemical_symbols()])[0]
        metal_slatm = local[metal_index[0]]
        if len(metal_index) > 1: 
            print("More than one metal atom in", names[i])
        #else:
        #    try:
        #        cm = qml.representations.generate_atomic_coulomb_matrix(atnums[i], mol.positions, maxsize, sorting="distance", indices=metal_index)
        #    except:
        #        print(metal_index, atnums[i])
    
            #if (np.isnan(full).any() or np.isnan(metal_slatm).any()):
            #    print("molecule with NaN in slatm",mol)
            #else:
            #    full = np.array(full)
        marray = np.array(metal_slatm)
        #np.save(names[i]+"_full", full)
        np.save(names[i]+"_local", marray)
        #np.save(names[i]+"_cm", cm)
        print(names[i], "done")
        #slatms.append(full)
        localslatms.append(marray)
        #cms.append(cm)

#slatmarray = np.array(slatms)
#np.save("full_array",slatmarray)

locarray = np.array(localslatms)
np.save(f"4_local/{names[i]}_local.npy",locarray)

#cmsarray = np.array(cms)
#np.save("cm_array",cmsarray)
