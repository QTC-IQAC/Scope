from cell2mol.tmcharge_common import *

def gmol_update_geom(mol: object, new_coord: list):
    mol.old_coord = mol.coord
    mol.coord = new_coord
    for idx, a in enumerate(mol.atoms):
        a.old_coord = a.coord
        a.coord = new_coord[idx]
