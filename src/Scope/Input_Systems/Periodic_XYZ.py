from Scope.Adapted_from_cell2mol import labels2formula, get_radii
from Scope import Software
from Scope.Software import Quantum_Espresso
from Scope.Software.Quantum_Espresso import Parse_QE_outputs

class periodic_xyz(object):
    def __init__(self, name: str, labels: list, coord: list, path: str) -> None:
        self.version              = "0.2"
        self.type                 = "perxyz"
        self.subtype              = "perxyz"
        self.name                 = name
        self.labels               = labels
        self.coord                = coord
        self.path                 = path
        self.moleclist            = []
        self.formula              = labels2formula(labels)
        self.phase                = str(name.split("_")[0])
        self.pressure             = float(name.split("_")[1])
        self.natoms               = len(labels)

    def add_cell_info(self, cellvec, celldim, cellparam) -> None:
        from Scope.Geometry  import get_unit_cell_volume
        self.cellvec              = cellvec
        self.celldim              = celldim
        self.cellparam            = cellparam
        self.volume               = get_unit_cell_volume(*cellparam)

    def add_molecule(self, mol: object) -> None:
        if mol not in self.moleclist: self.moleclist.append(mol)

    def sum_charges(self):
        self.totcharge = 0
        for mol in self.moleclist: self.totcharge += mol.totcharge
        return self.totcharge

    def __repr__(self):
        to_print = ""
        to_print  += f'---------- SCOPE Periodic xyz Object -----------\n'
        to_print += f' Version               = {self.version}\n'
        to_print += f' Type                  = {self.type}\n'
        to_print += f' Number of Atoms       = {self.natoms}\n'
        to_print += f' Formula               = {self.formula}\n'
        if hasattr(self,"phase"):    to_print += f' Phase                 = {self.phase}\n'
        if hasattr(self,"pressure"): to_print += f' Pressure              = {self.pressure}\n'
        to_print += '------------------------------------------------\n'
        return to_print

#########################
