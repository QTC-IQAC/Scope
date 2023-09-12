from Scope.Adapted_from_cell2mol import get_molecules
from Scope.Elementdata import ElementData
elemdatabase = ElementData()

from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal

##################
### COLLECTION ###
##################
class state(object):
    def __init__(self, _subject: object, name: str):
        self.type         = "state"
        self._subject     = _subject
        self.name         = name

        if not hasattr(self._subject,"states"): self._subject.states = []
        for idx, st in enumerate(self._subject.states):
            if st.name == name: self._subject.states[idx] = self  ## Replaces (updates) state?

    def set_geometry(self, labels, pos):
        self.labels      = labels
        self.pos         = pos 
         
    def set_cell(self, cellvec, cellparam):
        self.cellvec     = cellvec
        self.cellparam   = cellparam
        
    def set_moleclist(self):
        if hasattr(self,"labels") and hasattr(self,"pos"):
            if hasattr(self._subject,"factor"): factor = self._subject.factor
            else: factor = 1.3
            self.moleclist = get_molecules(self.labels, self.pos, factor=factor)
        else: 
            self.moleclist = []
        return self.moleclist
        
