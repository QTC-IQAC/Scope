from cell2mol.tmcharge_common import *
from Scope.VNM_tools import vnm_displacement 
#from copy import deepcopy

def cell_postocoord(cell: object) -> None:
    if not hasattr(cell,"coord") and hasattr(cell,"pos"): setattr(cell,"coord",cell.pos)

def gmol_update_geom(mol: object, new_coord: list, name: str="coord", debug: int=0) -> None:
    assert hasattr(mol,name)
    old_coord = getattr(mol,name)
    setattr(mol,"old_"+name,old_coord)
    setattr(mol,name,new_coord)
    for idx, at in enumerate(mol.atoms):
        old_coord = getattr(at,name)
        setattr(at,"old_"+name,old_coord)
        setattr(at,name,new_coord[idx])

def gmol_create_geom(mol: object, new_coord: list, name: str="new_coord", debug: int=0) -> None:
    setattr(mol,name,new_coord)
    for idx, at in enumerate(mol.atoms):
        setattr(at,name,new_coord[idx]

def cell_update_geom(cell: object, moleclist: list, new_coord: list, name: str="coord", debug: int=0) -> None:
    assert hasattr(cell,"cellparam")
    assert hasattr(cell,name)
    cell.old_moleclist = cell.moleclist
    cell.moleclist = moleclist
    cell.old_pos = cell.pos
    cell.pos = new_coord
    for idx, mol in enumerate(cell.moleclist):
        for jdx, at in enumerate(mol.atoms):
            setattr(at,name,new_coord[mol.atlist[jdx]])

    old_coord = getattr(cell,name)
    setattr(cell,"old_"+name,old_coord)
    setattr(cell,name,new_coord)
    for idx, mol in enumerate(cell.moleclist):
        for jdx, at in enumerate(mol.atoms):
            old_coord = getattr(at,name)
            setattr(at,"old_"+name,old_coord)
            setattr(at,name,new_coord[idx])

def cell_create_geom(cell: object, moleclist: list, new_coord: list, name: str="new_coord", debug: int=0) -> None:
    assert hasattr(cell,"moleclist")        ## Makes sure it is a cell
    setattr(cell,name,new_coord)            ## Creates new attribute to object cell, with name "name" and value "new_coord"
    for idx, mol in enumerate(cell.moleclist):
        for jdx, at in enumerate(mol.atoms):
            setattr(at,name,new_coord[mol.atlist[jdx]])

def displace_neg_freqs(gmol: object, debug: int=0) -> list:
    disp_coord = vnm_displacement(gmol.VNMs, gmol.coord)
    return disp_coord

