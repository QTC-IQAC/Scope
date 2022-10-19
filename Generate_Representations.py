#!/usr/bin/env python3
import os
import sys
import numpy as np
import qml
import time

import pickle
from collections import Counter
from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal
from cell2mol.elementdata import ElementData

###########################################

def get_aSlatm(mol, mbtypes, dens=0.4, cutoff=4.0, debug=0):

    if hasattr(mol, 'atnums') and hasattr(mol, 'coord') and hasattr(mol, 'atoms'):
        local = qml.representations.generate_slatm(mol.coord, mol.atnums, mbtypes=mbtypes, local=True, dgrids=[dens, dens], rcut=cutoff)
        if np.isnan(local).any(): 
            print("molecule with NaN in slatm", mol.refcode)

        metal_index = [] 
        for idx, a in enumerate(mol.atoms):
            if a.block == 'd': metal_index.append(idx) 
        metal_slatm = local[metal_index[0]]
        if len(metal_index) > 1: 
            print("More than one metal atom in", names[i])
        return np.array(metal_slatm)
    else: 
        if hasattr(mol, 'refcode'):
            print(f"mol with refcode {mol.refcode} is missing one or more of necessary properties")
        else:
            print(f"mol is missing one or more of necessary properties")

def get_aSlatm_ligands(mol, mbtypes, dens=0.4, cutoff=4.0, debug=0):

    if hasattr(mol, 'atnums') and hasattr(mol, 'coord') and hasattr(mol, 'atoms'):
        local = qml.representations.generate_slatm(mol.coord, mol.atnums, mbtypes=mbtypes, local=True, dgrids=[dens, dens], rcut=cutoff)
        if np.isnan(local).any():
            print("molecule with NaN in slatm", mol.refcode)
        connected_slatm = [] 
        for idx, a in enumerate(mol.atoms):
            if a.mconnec > 0: connected_slatm.append(local[idx])
        connected_slatm = np.sum(connected_slatm, axis=0)
        return np.array(connected_slatm)
    else:
        if hasattr(mol, 'refcode'):
            print(f"mol with refcode {mol.refcode} is missing one or more of necessary properties")
        else:
            print(f"mol is missing one or more of necessary properties")

