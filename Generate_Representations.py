#!/usr/bin/env python3
import os
import sys
import numpy as np
import time
try: 
    import qml
except: 
    pass

import pickle
from collections import Counter
from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal
#from Scope.Elementdata import ElementData
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

def get_aSlatm_ligands(mol, mbtypes, dens: float=0.4, cutoff: float=4.0, with_metal: bool=False, debug: int=0):

    if hasattr(mol, 'atnums') and hasattr(mol, 'coord') and hasattr(mol, 'atoms'):

        ## Move coords and atnums to independent lists in case we need to add the metal_atoms
        #print(mol.coord)

        coords = []
        atnums = []
        for idx, a in enumerate(mol.atoms):
            coords.append(a.coord)
            atnums.append(a.atnum)
        if with_metal == True:
            if len(mol.metalatoms) > 0:
                for m in mol.metalatoms:
                    atnums.append(m.atnum)
                    coords.append(m.coord)

        #print(mol.refcode, coords)
        #local = qml.representations.generate_slatm(mol.coord, mol.atnums, mbtypes=mbtypes, local=True, dgrids=[dens, dens], rcut=cutoff)
        local = qml.representations.generate_slatm(coords, atnums, mbtypes=mbtypes, local=True, dgrids=[dens, dens], rcut=cutoff)
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
