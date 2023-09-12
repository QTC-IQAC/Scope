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
        self.results      = dict()

        ## Creates the variable "states" to the _subject, which should be a "cell" or "molecule"-class object
        if not hasattr(self._subject,"states"): self._subject.states = []
        updated = False
        for idx, st in enumerate(self._subject.states):
            if st.name == name: self._subject.states[idx] = self; updated = True; print("UPDATED")  ## Replaces (updates) state?
        if not updated: self._subject.states.append(self); print("ADDED")

    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  self.results[result.key] = result

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

    def link_to_computation(self, _computation: object):
        self._computation = _computation
        
