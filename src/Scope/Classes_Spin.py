import numpy as np
from Scope.Other import get_metal_idxs

#################
### SPIN INFO ###
#################
class spin_config(object):
    def __init__(self, _source: object):
        self.type           = "spin_state"
        self._source        = _source
        self.atomic_spins   = []

    def add_atomic_spin(self, label, index, spin):
        new_atomic_spin = atomic_spin(label, index, spin, _spin_config=self)
        self.atomic_spins.append(new_atomic_spin)
        return new_atomic_spin 

    def get_total_magnetization(self):
        self.total_magnetization = 0
        self.ismagnetic = False
        for idx, atom_spin in enumerate(self.atomic_spins):
            self.total_magnetization += atom_spin.magnetization
        if self.total_magnetization != 0: self.ismagnetic = True
        return self.total_magnetization
   
    def get_multiplicity(self):
        self.multiplicity = 0
        for idx, atom_spin in enumerate(self.atomic_spins):
            if   atom_spin.orientation == "up":   self.multiplicity += atom_spin.multiplicity 
            elif atom_spin.orientation == "down": self.multiplicity -= atom_spin.multiplicity
        if len(self.atomic_spins) == 0 and self.multiplicity == 0: self.multiplicity = 1   ## For organic molecules
        return self.multiplicity    

    def get_QE_data(self):
        self.elems =  list(set(self._source.labels))
        self.nelems = len(self.elems)
        if len(self.atomic_spins) == 0:
            self.magn_pairs = []
            self.magn_uniques = []
        else:
            self.magn_pairs = list(set(tuple([atsp.label, atsp.magnetization]) for atsp in self.atomic_spins)) 
            self.magn_uniques = list(set(self.magn_pairs))
 
    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   SCOPEs Spin Configuration Class                 \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Source Name                  = {self._source.name}\n'
        to_print += f' Source Type                  = {self._source.type}\n'
        to_print += f'---------------------------------------------------\n'
        if hasattr(self,"elems"):               to_print += f' Elements                     = {self.elems}\n'
        if hasattr(self,"ismagnetic"):          to_print += f' Is Magnetic?                 = {self.ismagnetic}\n'
        if hasattr(self,"multiplicity"):        to_print += f' Multiplicity                 = {self.multiplicity}\n'
        if hasattr(self,"total_magnetization"): to_print += f' Total Magnetization          = {self.total_magnetization}\n'
        to_print += '----------------------------------------------------\n'
        return to_print

###################
### ATOMIC SPIN ###
###################
class atomic_spin(object):
    def __init__(self, label, index, spin, _spin_config):
        self.type           = "atomic_spin"
        self.label          = label
        self.index          = index
        self.spin           = spin
        self._spin_config   = _spin_config
        self.unpaired_elec  = get_unpaired_elec(label, spin)
        self.orientation    = "up"
        
        self.ms             = float(self.unpaired_elec/2)
        self.magnetization  = int(self.unpaired_elec)
        self.multiplicity   = int(2*self.ms + 1)

    def get_mod_label(self):
        if hasattr(self, "spin"):
            if   self.spin == "HS":       suffix = str("4")
            elif self.spin == "LS":       suffix = str("0")
            elif self.spin == "IS":       suffix = str("2")
            return self.label + suffix 
        elif hasattr(self, "magnetization"):
            if   self.magnetization == 4: suffix = str("4")
            elif self.magnetization == 0: suffix = str("0")
            elif self.magnetization == 2: suffix = str("2")
            return self.label + suffix 
        else: print("ATOMIC SPIN: I do not have information of either spin or magnetization")

    def set_orientation(self, orientation: str="up"):
        if   orientation.lower() == "up":   self.orientation = "up" 
        elif orientation.lower() == "down": self.orientation = "down"
        else: print("SET_ORIENTATION: I do not understand orientation"); return None 
        return self.orientation

#################
def get_unpaired_elec(label, spin):
    if   label == "Fe" and spin == "HS": return int(4)
    elif label == "Fe" and spin == "IS": return int(2)
    elif label == "Fe" and spin == "LS": return int(0)
    else: print("get_unpaired_elec: label and/or spin not in library"); return None

#################
def get_spin_config(source: object, metal_spins: list, debug: int=0):

    if debug > 0: print(f"GET_SPIN_CONFIG: Preparing Spin Configuration for New Computation Involving source: {source.name}")
    if debug > 0: print(f"GET_SPIN_CONFIG: Received metal_spins", metal_spins)
    assert hasattr(source,"labels"), f"GET_SPIN_CONFIG got object without labels"

    #########################
    ### IDENTIFIES METALS ###
    #########################
    metal_indices = get_metal_idxs(source.labels)

    ## if the user provides an abbreviated list of spin states. For instance, metal_spins="HS"
    if type(metal_spins) == list:
        if len(metal_spins) == 1:   # user sends 'LS'/'HS'
            tmp = metal_spins[0]
            is_abbr = True
        if len(metal_spins) == 0:   # user sends 'LS'/'HS' but there is no metal. Then assume 'LS'
            tmp = 'LS'
            is_abbr = True
    elif type(metal_spins) == str:
        is_abbr = True
        tmp = metal_spins

    if is_abbr:
        metal_spins = []
        for i in range(len(metal_indices)):
            metal_spins.append(tmp)

    ## Create spin_config-class object and fill it with spins
    pointer = 0
    new_spcf = spin_config(source) 
    for idx, l in enumerate(source.labels):
        if idx in metal_indices:
            desired_spin = metal_spins[pointer]
            new_atomic_spin = new_spcf.add_atomic_spin(l, idx, desired_spin)
            new_atomic_spin.get_mod_label
            pointer += 1 
    new_spcf.get_total_magnetization()
    new_spcf.get_multiplicity()

    return new_spcf
