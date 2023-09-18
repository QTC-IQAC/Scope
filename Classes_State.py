from Scope.Adapted_from_cell2mol import get_molecules, printxyz
from Scope.Classes_Data import collection, data
from Scope.Elementdata import ElementData
elemdatabase = ElementData()

#from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal

##############
### STATES ###
##############
class state(object):
    def __init__(self, _subject: object, name: str, debug: int=0):
        self.type         = "state"
        self._subject     = _subject
        self.name         = name
        self.results      = dict()
        self.computations = []

#        ## Creates the variable "states" to the _subject, which should be a "cell" or "molecule"-class object
        if not hasattr(self._subject,"states"): self._subject.states = []
        found = False
        for idx, st in enumerate(self._subject.states):
            if st.name == name: found = True; return st
        if not found: self._subject.states.append(self)

#        updated = False
#        for idx, st in enumerate(self._subject.states):
#            if st.name == name: 
#                st._update(self, debug=debug)
#                if debug > 0: print("UPDATED SELF", dir(self))
#                updated = True
#                if debug > 0: print("UPDATED state", name) 
#        if not updated: 
#            self._subject.states.append(self)
#            if debug > 0: print("ADDED state", name)

    def _update(self, other, debug: int=0):
        if debug > 0: print("Updating State", self.name)
        if debug > 0: print("SELF", dir(self))
        if debug > 0: print("OTHER", dir(other))
        if not isinstance(other, type(self)): 
            if debug > 0: print("_update STATE: other is not of the same type than self")
            return self
        for d in dir(other):
            if debug > 0: print("_update STATE: checking instance", d)
            if d in dir(other):
                if '_' not in d and not callable(getattr(other,d)) and d != "type" and d != "name" and d != "results" and d != "computations":
                    if debug > 0: print("_update STATE: updating instance", d)
                    at1 = getattr(other,d)
                    try:      attr = literal_eval(at1)
                    except:   attr = value
                    setattr(self, d, attr)
                elif d == "results":
                    for r in other.results:
                        if debug > 0: print("_update STATE: updating result", r)
                        self.add_result(r, overwrite=True)
                elif d == "computations": 
                    for c in other.computations:
                        if c not in self.computations: 
                            if debug > 0: print("_update STATE: updating computations, adding", c)
                            self.computations.append(c)
        return self

    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  self.results[result.key] = result

    def set_geometry(self, labels, pos):
        self.labels      = labels
        self.pos         = pos 
        self.coord       = pos 
        self.natoms      = len(labels)
        assert len(self.labels) == len(self.pos)

    def set_geometry_from_moleclist(self):
        self.labels      = []
        self.pos         = []
        self.coord       = []
        indices     = []
        for mol in self.moleclist:
            for idx, at in enumerate(mol.atoms):
                self.labels.append(at.label)
                self.pos.append(at.coord)
                self.coord.append(at.coord)
                indices.append(mol.atlist[idx])
        ## Below is to order the atoms as in the original cell, using the indices stored in the molecule object
        self.labels = [x for _, x in sorted(zip(indices, self.labels), key=lambda pair: pair[0])]
        self.pos = [x for _, x in sorted(zip(indices, self.pos), key=lambda pair: pair[0])]
        self.coord = [x for _, x in sorted(zip(indices, self.coord), key=lambda pair: pair[0])]
        self.natoms = len(self.labels)
        assert len(self.labels) == len(self.pos)
         
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
            if len(self.labels) > 0 and len(self.pos) > 0:
                if hasattr(self._subject,"factor"): factor = self._subject.factor
                else: factor = 1.3
                try: 
                    self.moleclist = get_molecules(self.labels, self.pos, factor=factor)
                except: 
                    print("ERROR setting moleclist for state:", self.name)
                    printxyz(self.labels, self.pos)
        else: 
            self.moleclist = []
        return self.moleclist

    def set_spin_config(self, spin_config):
        self.spin_config = spin_config

    def set_atoms(self):
        if not hasattr(self,"moleclist"): self.set_moleclist()
        self.atoms = []
        indices = []
        for mol in self.moleclist:
            for idx, at in enumerate(mol.atoms):
                self.atoms.append(at)
                indices.append(mol.atlist[idx])
        self.atoms = [x for _, x in sorted(zip(indices, self.atoms), key=lambda pair: pair[0])]

    def set_energy(self, energy, units):
        self.add_result(data("Energy",energy,units,"state.set_energy"))
        #self.energy     = energy

    def add_computation(self, computation: object):
        self.computations.append(computation)

    def reconstruct(self, debug: int=0):
        if hasattr(self._subject,"type") and hasattr(self,"cellvec"):
            if not hasattr(self,"moleclist"): self.set_moleclist() 
            if self._subject.type.lower() == "cell":
                from cell2mol import cell_reconstruct
                from cell2mol.cell_reconstruct import identify_frag_molec_H, getmolecs, tmatgenerator, additem, assigntype, fragments_reconstruct
                from cell2mol.cellconversions import frac2cart_fromparam, cart2frac, translate
                import itertools
                blocklist = self.moleclist.copy()
                moleclist = []
                refmoleclist = self._subject.refmoleclist.copy()
                covalent_factor = refmoleclist[0].factor
                metal_factor = refmoleclist[0].metal_factor
                moleclist, fraglist, Hlist, init_natoms = identify_frag_molec_H(blocklist, moleclist, refmoleclist, self.cellvec) 
                if len(fraglist) > 0 or len(Hlist) > 0: 
                    moleclist, finalmols, Warning = fragments_reconstruct(moleclist,fraglist,Hlist,refmoleclist,self.cellvec,covalent_factor,metal_factor)
                    moleclist.extend(finalmols)
                    self.moleclist = moleclist
                    self.set_geometry_from_moleclist()
 
            else: print("WARNING: reconstruct state, _subject is not a cell. I will not reconstruct")
        else: print("WARNING: reconstruct state, _subject does not have 'type' or 'cellvec' variables")
        return self.moleclist


    def __repr__(self) -> None:
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   STATE                                           \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Name                  = {self.name}\n'
        if hasattr(self,"labels"): to_print += f' Labels                = {self.labels}\n'
        if hasattr(self,"coord"):  to_print += f' Coord                 = {self.coord}\n'
        return to_print

#########################################################################
## Tools associated with states. Normally, these would be class functions...
## However, in this case the classes are defined in cell2mol, and I don't want to change them 
#########################################################################

def find_state(subject: object, search_name: str, debug: int=0):
    if debug >= 1: print("FIND_STATE: enters",search_name," with", len(subject.states),"states in subject")
    if not hasattr(subject,"states"): return False, None
    else: 
        found = False
        for idx, sta in enumerate(subject.states):
            if sta.name == search_name: 
                found = True
                if debug >= 1: print("FIND STATE: state found")
                return True, sta
        if not found: 
            return False, None
