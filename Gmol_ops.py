from cell2mol.tmcharge_common import *

def gmol_update_geom(mol: object, new_coord: list, debug: int=0):
    mol.old_coord = mol.coord
    mol.coord = new_coord
    for idx, a in enumerate(mol.atoms):
        a.old_coord = a.coord
        a.coord = new_coord[idx]

def cell_update_geom(cell: object, new_coord: object, debug: int=0):
    assert hasattr(cell,"cellparam")
    cell.old_moleclist = cell.moleclist.deepcopy()
    cell.moleclist = moleclist.deepcopy()
    cell.old_pos = pos.deepcopy()
    cell.pos = pos

