from cell2mol.tmcharge_common import *
from copy import deepcopy

def gmol_update_geom(mol: object, new_coord: list, debug: int=0):
    mol.old_coord = mol.coord
    mol.coord = new_coord
    for idx, a in enumerate(mol.atoms):
        a.old_coord = a.coord
        a.coord = new_coord[idx]

def cell_update_geom(cell: object, moleclist: list, new_coord: list, debug: int=0):
    assert hasattr(cell,"cellparam")
    cell.old_moleclist = cell.moleclist
    cell.moleclist = moleclist
    cell.old_pos = cell.pos
    cell.pos = new_coord

