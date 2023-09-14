from Scope.Adapted_from_cell2mol import get_molecules
from Scope.Elementdata import ElementData
elemdatabase = ElementData()

from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal

##############
### STATES ###
##############
class state(object):
    def __init__(self, _subject: object, name: str):
        self.type         = "state"
        self._subject     = _subject
        self.name         = name
        self.results      = dict()
        self.computations = []

        ## Creates the variable "states" to the _subject, which should be a "cell" or "molecule"-class object
        if not hasattr(self._subject,"states"): self._subject.states = []
        updated = False
        for idx, st in enumerate(self._subject.states):
            if st.name == name: self._subject.states[idx] = self; updated = True #; print("UPDATED")  ## Replaces (updates) state?
        if not updated: self._subject.states.append(self) #; print("ADDED")

    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  self.results[result.key] = result

    def set_geometry(self, labels, pos):
        self.labels      = labels
        self.pos         = pos 
        self.coord       = pos 
        for idx, at in enumerate(_subject.atoms):
            setattr(at,self.name,    
         
    def set_cell(self, cellvec, cellparam):
        self.cellvec     = cellvec
        self.cellparam   = cellparam

    def set_VNMs(self, VNMs):
        self.VNMs = VNMs
        self.freqs_cm = [vnm.freq_cm for vnm in VNMs]
        if all(vnm.freq >= 0.0 for vnm in self.VNMs): self.isminimum = True
        else:                                         self.isminimum = False
        
    def set_moleclist(self):
        if hasattr(self,"labels") and hasattr(self,"pos"):
            if hasattr(self._subject,"factor"): factor = self._subject.factor
            else: factor = 1.3
            self.moleclist = get_molecules(self.labels, self.pos, factor=factor)
        else: 
            self.moleclist = []
        return self.moleclist

    def set_energy(self, energy):
        self.energy     = energy

    def add_computation(self, computation: object):
        self.computations.append(computation)
        
#########################################################################
## Tools associated with states. Normally, these would be class functions...
## However, in this case the classes are defined in cell2mol, and I don't want to change them 
#########################################################################

def find_state(subject: object, search_name: str):
    if hasattr(subject,"states"):
        found = False
        for idx, sta in subject.states:
            if sta.name == search_name: found = True; return sta 
        if not found: 
            sta = state(subject, search_name)
            subject.states.append(sta)
            return sta
    else: 
        setattr(subject,"states",[])
        sta = state(subject, search_name)
        subject.states.append(sta)
        return sta 
