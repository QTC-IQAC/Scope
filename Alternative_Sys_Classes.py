from Scope.Adapted_from_cell2mol import labels2formula

from Scope import Software
from Scope.Software import Quantum_Espresso
from Scope.Software.Quantum_Espresso import Parse_QE_outputs

#########################
class simple_molecule(object):
    def __init__(self, atom_idx: list, labels: list, coord: list, radii: list) -> None:
        self.atom_idx             = atom_idx
        self.labels               = labels
        self.coord                = coord
        self.radii                = radii

class periodic_xyz(object):
    def __init__(self, name: str, labels: list, coord: list, path: str) -> None:
        self.name                 = name
        self.labels               = labels
        self.coord                = coord
        self.path                 = path
        self.moleclist            = []
        self.formula              = labels2formula(labels)
        self.phase                = str(name.split("_")[0])
        self.pressure             = float(name.split("_")[1])

    def add_cell_info(self, cellvec, celldim, cellparam) -> None:
        self.cellvec              = cellvec
        self.celldim              = celldim
        self.cellparam            = cellparam
        self.volume               = get_unit_cell_volume(*cellparam)

    def add_molecule(self, mol: object) -> None:
        if mol not in self.moleclist: self.moleclist.append(mol)
#########################
