#from cell2mol.tmcharge_common import *
from Scope.VNM_tools import vnm_displacement 
from Scope.Workflow import Branch
from Scope.Workflow.Branch import *
from Scope.Classes_State import *
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

def gmol_remove_geom(mol: object, tag: str="new_coord", debug: int=0) -> None:
    if hasattr(mol,tag):
        delattr(mol,tag)
        for idx, at in enumerate(mol.atoms):
            if hasattr(at,tag):
                delattr(at,tag)

#def cell_update_geom(cell: object, moleclist: list, new_coord: list, tag: str="coord", debug: int=0) -> None:
#    assert hasattr(cell,"cellparam")
#    assert hasattr(cell,tag)
#   
#    ## Updates Moleclist
#    setattr(cell,"old_moleclist",cell.moleclist)
#    cell.moleclist = moleclist
#
#    ## Updates Coordinates
#    old_coord = getattr(cell,tag)
#    setattr(cell,"old_"+tag,old_coord)
#    setattr(cell,tag,new_coord)
#    for idx, mol in enumerate(cell.moleclist):
#        print("doing molecule", idx)
#        for jdx, at in enumerate(mol.atoms):
#            old_coord = getattr(at,tag)
#            setattr(at,"old_"+tag,old_coord)
#            setattr(at,tag,new_coord[idx])
#
#def cell_create_geom(cell: object, moleclist: list, new_coord: list, tag: str="new_coord", debug: int=0) -> None:
#    assert hasattr(cell,"moleclist")        ## Makes sure it is a cell
#    setattr(cell,tag,new_coord)            ## Creates new attribute to object cell, with tag "tag" and value "new_coord"
#    for idx, mol in enumerate(cell.moleclist):
#        for jdx, at in enumerate(mol.atoms):
#            setattr(at,tag,new_coord[mol.atlist[jdx]])

def displace_neg_freqs(ini_coord, VNMs: object, debug: int=0) -> list:
    neg_VNMs = list([vnm for vnm in initial_state.VNMs if vnm.freq_cm < 0.0])
    disp_coord = vnm_displacement(neg_VNMs, ini_coord)
    return disp_coord

def cell_update_geom(cell: object, new_coord: list, tag: str="coord", debug: int=0) -> None:
    assert hasattr(cell,"cellparam")
    assert hasattr(cell,tag)
   
    ## Updates Coordinates
    old_coord = getattr(cell,tag)
    setattr(cell,"old_"+tag,old_coord)
    setattr(cell,tag,new_coord)

#    for idx, mol in enumerate(cell.moleclist):
#        print("doing molecule", idx)
#        for jdx, at in enumerate(mol.atoms):
#            old_coord = getattr(at,tag)
#            setattr(at,"old_"+tag,old_coord)
#            setattr(at,tag,new_coord[mol.atlist[idx]])
#            assert at.label == cell.labels[mol.atlist[idx]]

def cell_create_geom(cell: object, new_coord: list, tag: str="new_coord", debug: int=0) -> None:
    #assert hasattr(cell,"moleclist")        ## Makes sure it is a cell
    setattr(cell,tag,new_coord)            ## Creates new attribute to object cell, with tag "tag" and value "new_coord"
    #for idx, mol in enumerate(cell.moleclist):
    #    for jdx, at in enumerate(mol.atoms):
    #        setattr(at,tag,new_coord[mol.atlist[jdx]])


#######################
### Molecules as Sys ##
#######################
def gmol2sys(gmol: object, spin: str, name: str):
    gmol.name        = name
    gmol.refcode     = name
    gmol.branches    = []
    gmol._sys        = gmol
    gmol.spin        = 'LS'
    return gmol

def find_branch_gmol(gmol, keyword: str, debug: int=0):
    if not hasattr(gmol,"branches"): setattr(gmol,"branches",[])
    if debug > 1: print("finding branch with keyword:", keyword)
    if debug > 1: print("there are", len(gmol.branches), "branches in system")
    if len(gmol.branches) == 0: return False, None
    for idx, br in enumerate(gmol.branches):
        if debug > 1: print("evaluating branch with keyword:", br.keyword, "and path:", br.path)
        if br.keyword == keyword:
            if not os.path.isdir(br.path): print(f"WARNING: branch path does not exist. Loading the branch anyway")
            return True, br
    return False, None

def add_branch_gmol(gmol, keyword: str, calcs_path: str, debug: int=0):
    if calcs_path[-1] != '/': calcs_path += '/'
    new_branch = Branch.branch(calcs_path+keyword, keyword, gmol, debug=debug)
    if not os.path.isdir(calcs_path+keyword):
        try: os.makedirs(calcs_path+keyword)
        except Exception as exc:
             print(f"Error creating branch folder in {calcs_path+keyword}")
             print(exc)

    ## Creates recipes for the branch. One for each object.
    new_recipe = new_branch.add_recipe(gmol)
    gmol.branches.append(new_branch)

    new_state = state(gmol,"initial") 
    new_state.set_geometry(gmol.labels, gmol.coord)

    return new_branch

##########

##########################
### Periodic xyz as Sys ## (see Alternative_Sys_Classes.py for definition of a periodic xyz)
##########################
def perxyz2sys(perxyz: object, name: str):
    perxyz.name        = name
    perxyz.refcode     = name
    perxyz.branches    = []
    perxyz._sys        = perxyz
    perxyz.type        = "perxyz"

    if   perxyz.phase == "HT": perxyz.spin = "HS"
    elif perxyz.phase == "LT": perxyz.spin = "LS"
    else:                      perxyz.spin = "LS"
    return perxyz

def find_branch_perxyz(perxyz, keyword: str, debug: int=0):
    if not hasattr(perxyz,"branches"): setattr(perxyz,"branches",[])
    if debug > 1: print("finding branch with keyword:", keyword)
    if debug > 1: print("there are", len(perxyz.branches), "branches in system")
    if len(perxyz.branches) == 0: return False, None
    for idx, br in enumerate(perxyz.branches):
        if debug > 1: print("evaluating branch with keyword:", br.keyword, "and path:", br.path)
        if br.keyword == keyword:
            if not os.path.isdir(br.path): print(f"WARNING: branch path does not exist. Loading the branch anyway")
            return True, br
    return False, None

def add_branch_perxyz(perxyz, keyword: str, calcs_path: str, debug: int=0):
    if calcs_path[-1] != '/': calcs_path += '/'
    new_branch = Branch.branch(calcs_path+keyword, keyword, perxyz, debug=debug)
    if not os.path.isdir(calcs_path+keyword):
        try: os.makedirs(calcs_path+keyword)
        except Exception as exc:
             print(f"Error creating branch folder in {calcs_path+keyword}")
             print(exc)

    ## Creates recipes for the branch. One for each object.
    new_recipe = new_branch.add_recipe(perxyz)
    perxyz.branches.append(new_branch)

    new_state = state(perxyz,"initial")
    new_state.set_geometry(perxyz.labels, perxyz.coord)
    new_state.set_cell(perxyz.cellvec, perxyz.cellparam)

    return new_branch

#############################
### Simple Molecule as Sys ## (see Alternative_Sys_Classes.py for definition of a simple molecule)
#############################
def smol2sys(smol: object, name: str):
    smol.type        = "smol"
    smol.name        = name
    smol.refcode     = name
    smol.branches    = []
    smol._sys        = smol
    return smol

def find_branch_smol(smol, keyword: str, debug: int=0):
    if not hasattr(smol,"branches"): setattr(smol,"branches",[])
    if debug > 1: print("finding branch with keyword:", keyword)
    if debug > 1: print("there are", len(smol.branches), "branches in system")
    if len(smol.branches) == 0: return False, None
    for idx, br in enumerate(smol.branches):
        if debug > 1: print("evaluating branch with keyword:", br.keyword, "and path:", br.path)
        if br.keyword == keyword:
            if not os.path.isdir(br.path): print(f"WARNING: branch path does not exist. Loading the branch anyway")
            return True, br
    return False, None

def add_branch_smol(smol, keyword: str, calcs_path: str, debug: int=0):
    if calcs_path[-1] != '/': calcs_path += '/'
    new_branch = Branch.branch(calcs_path+keyword, keyword, smol, debug=debug)
    if not os.path.isdir(calcs_path+keyword):
        try: os.makedirs(calcs_path+keyword)
        except Exception as exc:
             print(f"Error creating branch folder in {calcs_path+keyword}")
             print(exc)

    ## Creates recipes for the branch. One for each object.
    new_recipe = new_branch.add_recipe(smol)
    smol.branches.append(new_branch)

    new_state = state(smol,"initial")
    new_state.set_geometry(smol.labels, smol.coord)

    return new_branch
