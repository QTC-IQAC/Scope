#!/usr/bin/env python3
import os
import sys
import numpy as np
import time

import pickle
from collections import Counter
from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal
#from Test_V3.Elementdata import ElementData
#elemdatabase = ElementData()

###########################################
def atoms_for_representation(mol: object, target: list=['metal'], debug=0):
    target = list(bla.lower() for bla in target) 
    indices = []
    types = []
    if "metal" in target or "all" in target:
        for idx, a in enumerate(mol.atoms):
            if a.block == 'd': indices.append(idx)
            if a.block == 'd': types.append('metal')
    elif "ligands" in target or "all" in target:
        for idx, a in enumerate(mol.atoms):
            if a.mconnec > 0: indices.append(idx)
            if a.mconnec > 0: types.append('ligands')
    else:
        print("Target not identified. Options are METAL, LIGAND, or ALL")
    return indices, types 

def get_aSlatm(mol, indices, mbtypes, targ_coord: str="coord", dens=0.4, cutoff=4.0, debug=0):
    import qml

    if hasattr(mol, 'atnums') and hasattr(mol, targ_coord) and hasattr(mol, 'atoms'):
        coord = getattr(mol,targ_coord)  ## Gets the desired coordinates
        all_local = qml.representations.generate_slatm(coord, mol.atnums, mbtypes=mbtypes, local=True, dgrids=[dens, dens], rcut=cutoff)
        if np.isnan(all_local).any(): print("molecule with NaN in slatm", mol.refcode)

        aSlatm = []
        for idx, local in enumerate(all_local):
            if idx in indices: aSlatm.append(local)
        return np.array(aSlatm)
    else: 
        if hasattr(mol, 'refcode'):
            print(f"mol with refcode {mol.refcode} is missing one or more of necessary properties")
        else:
            print(f"mol is missing one or more of necessary properties")
        return None

#def get_SOAP(mol, indices, targ_coord: str="coord", debug=0):
#    import dscribe
