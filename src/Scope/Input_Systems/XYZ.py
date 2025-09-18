from Scope.Adapted_from_cell2mol import labels2formula, get_radii
from Scope import Software
from Scope.Software import Quantum_Espresso
from Scope.Software.Quantum_Espresso import Parse_QE_outputs

#########################
class simple_molecule(object):
    def __init__(self, atom_idx: list, labels: list, coord: list, totcharge: int=0, spin: str='LS') -> None:
        self.version              = "0.2"
        self.type                 = "smol"
        self.subtype              = "smol"
        self.atom_idx             = atom_idx
        self.labels               = labels
        self.coord                = coord
        self.radii                = get_radii(labels)
        self.totcharge            = totcharge
        self.spin                 = spin  

