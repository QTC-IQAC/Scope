
from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal
from cell2mol.tmcharge_common import labels2formula
from cell2mol.elementdata import ElementData
Elemdatabase = ElementData()

#################################
def compare_species(mol1, mol2):
    # a pair of species is compared on the basis of:
    # 1) the total number of atoms
    # 2) the number of atoms of each type
    # 3) the number of adjacencies between each pair of element types
    issame = False
    if (mol1.natoms == mol2.natoms):
        if (mol1.eleccount == mol2.eleccount):
            # Checks element count
            issame = True
            for kdx, elem in enumerate(mol1.elemcountvec):
                if elem != mol2.elemcountvec[kdx]:
                    issame = False
            # Checks adjacency type count
            if hasattr(mol1,"adjtypes") and hasattr(mol2,"adjtypes") and issame:
                for kdx, elem in enumerate(mol1.adjtypes):
                    for ldx, elem2 in enumerate(elem):
                        if elem2 != mol2.adjtypes[kdx, ldx]:
                            issame = False
    return issame
#################################
