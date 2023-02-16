from cell2mol.tmcharge_common import *
from Scope.VNM_tools import vnm_displacement 
#from copy import deepcopy

def gmol_update_geom(mol: object, new_coord: list, tag: str="coord", debug: int=0) -> None:
    assert hasattr(mol,tag)
    old_coord = getattr(mol,tag)
    setattr(mol,"old_"+tag,old_coord)
    setattr(mol,tag,new_coord)
    for idx, at in enumerate(mol.atoms):
        old_coord = getattr(at,tag)
        setattr(at,"old_"+tag,old_coord)
        setattr(at,tag,new_coord[idx])

def gmol_create_geom(mol: object, new_coord: list, tag: str="new_coord", debug: int=0) -> None:
    setattr(mol,tag,new_coord)
    for idx, at in enumerate(mol.atoms):
        setattr(at,tag,new_coord[idx])

def cell_update_geom(cell: object, moleclist: list, new_coord: list, tag: str="coord", debug: int=0) -> None:
    assert hasattr(cell,"cellparam")
    assert hasattr(cell,tag)
    cell.old_moleclist = cell.moleclist
    cell.moleclist = moleclist
    cell.old_pos = cell.pos
    cell.pos = new_coord
    for idx, mol in enumerate(cell.moleclist):
        for jdx, at in enumerate(mol.atoms):
            setattr(at,tag,new_coord[mol.atlist[jdx]])

    old_coord = getattr(cell,tag)
    setattr(cell,"old_"+tag,old_coord)
    setattr(cell,tag,new_coord)
    for idx, mol in enumerate(cell.moleclist):
        for jdx, at in enumerate(mol.atoms):
            old_coord = getattr(at,tag)
            setattr(at,"old_"+tag,old_coord)
            setattr(at,tag,new_coord[idx])

def cell_create_geom(cell: object, moleclist: list, new_coord: list, tag: str="new_coord", debug: int=0) -> None:
    assert hasattr(cell,"moleclist")        ## Makes sure it is a cell
    setattr(cell,tag,new_coord)            ## Creates new attribute to object cell, with tag "tag" and value "new_coord"
    for idx, mol in enumerate(cell.moleclist):
        for jdx, at in enumerate(mol.atoms):
            setattr(at,tag,new_coord[mol.atlist[jdx]])

def displace_neg_freqs(gmol: object, ini_coord_tag: str="coord", debug: int=0) -> list:
    ini_coord = getattr(gmol,ini_coord_tag)
    disp_coord = vnm_displacement(gmol.VNMs, ini_coord)
    return disp_coord

