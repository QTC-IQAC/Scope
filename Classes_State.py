import numpy as np
from Scope.Adapted_from_cell2mol import *
from Scope.Classes_Data import collection, data
from Scope.Classes_Molecule import *
from Scope.Reconstruct    import *
from Scope.Elementdata import ElementData
elemdatabase = ElementData()

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

        ## Creates the variable "states" to the _subject, which should be a "cell" or "molecule"-class object
        if not hasattr(self._subject,"states"): self._subject.states = []
        found = False
        
        ## Verifies that a state with the same name does not exist already. If not, appends it
        for idx, st in enumerate(self._subject.states):
            if st.name == name: 
                found = True
                print("WARNING from CLASS STATE: you're trying to create a state that already exists in _subject.")
                print("WARNING from CLASS STATE: use function called 'find_state' instead to retrieve existing state")
        if not found: self._subject.states.append(self)

    def add_result(self, result: object, overwrite: bool=False):
        result._object = self
        if overwrite or result.key not in self.results.keys():  
            self.results[result.key] = result

    def set_geometry(self, labels, pos):
        self.labels      = labels
        self.pos         = pos 
        self.coord       = pos 
        self.natoms      = len(labels)
        self.formula     = labels2formula(self.labels)
        self.radii       = get_radii(labels)
        assert len(self.labels) == len(self.pos)

    def set_geometry_from_moleclist(self):
        self.labels      = []
        self.pos         = []
        self.coord       = []
        indices     = []
        for mol in self.moleclist:
            if not hasattr(mol,"atoms"): mol.set_atoms()
            for idx, at in enumerate(mol.atoms):
                self.labels.append(at.label)
                self.pos.append(at.coord)
                self.coord.append(at.coord)
                indices.append(mol.indices[idx])
        ## Below is to order the atoms as in the original cell, using the indices stored in the molecule object
        self.labels = [x for _, x in sorted(zip(indices, self.labels), key=lambda pair: pair[0])]
        self.pos = [x for _, x in sorted(zip(indices, self.pos), key=lambda pair: pair[0])]
        self.coord = [x for _, x in sorted(zip(indices, self.coord), key=lambda pair: pair[0])]
        self.natoms = len(self.labels)
        self.formula = labels2formula(self.labels)
        assert len(self.labels) == len(self.pos)
         
    def set_cell(self, cellvec, cellparam):
        self.cellvec     = cellvec
        self.cellparam   = cellparam

    def set_forces(self, forces):
        self.forces      = forces

    def set_VNMs(self, VNMs):
        self.VNMs = VNMs
        self.freqs_cm = [vnm.freq_cm for vnm in VNMs]
        if all(vnm.freq_cm >= 0.0 for vnm in self.VNMs): self.isminimum = True
        else:                                            self.isminimum = False

        ## If it is not a minimum, evaluates if, at least, is close
        if not self.isminimum:
            self.num_neg_freqs = 0
            for vnm in self.VNMs: 
                if vnm.freq_cm < 0.0: self.num_neg_freqs += 1
            if self.num_neg_freqs <= 3 and VNMs[0].freq_cm > -50: self.almost_minimum = True
            else:                                                 self.almost_minimum = False
        
    def get_moleclist(self):
        from Scope.Classes_Molecule import specie, molecule
        if not hasattr(self,"labels") or not hasattr(self,"pos"): return None
        if len(self.labels) == 0 or len(self.pos) == 0: return None
        if hasattr(self._subject,"factor"):       cov_factor = self._subject.factor
        elif hasattr(self._subject,"cov_factor"): cov_factor = self._subject.cov_factor
        else: cov_factor = 1.3

        blocklist = split_species(self.labels, self.pos, cov_factor=cov_factor)
        self.moleclist = [] 
        for b in blocklist:
            mol_labels  = extract_from_list(b, self.labels, dimension=1)
            mol_coords  = extract_from_list(b, self.coord, dimension=1)
            #mol_radii   = extract_from_list(b, self.radii, dimension=1)
            newmolec    = molecule(mol_labels, mol_coords)
            if newmolec.iscomplex: newmolec.split_complex()
            if self._subject.type == "cell": newmolec.get_fractional_coord(cell_vector = self._subject.cellvec)
            self.moleclist.append(newmolec)
        return self.moleclist

    def get_ncomplex(self, debug: int=0):
        if not hasattr(self,"fragmented"): self.check_fragmentation(reconstruct=True, debug=debug)
        assert not self.fragmented, f"Found Fragmented molecules in the geometry of state: {self.name}"
         
        if not hasattr(self,"moleclist"): self.get_moleclist()
        self.ncomplex = 0
        for mol in self.moleclist:
            if mol.iscomplex: self.ncomplex += 1
        if debug > 0: print(f"State.get_ncomplex {self.ncomplex} complexes found in state: {self.name}")
        return self.ncomplex

#######
#    Not sure it is correct to extract Z like that
#######
#    def get_Z(self):
#        if not hasattr(self,"labels"): return None
#        import numpy as np
#        ratio = labels2ratio(self.labels)
#        self.Z = np.gcd.reduce(ratio)
#        return self.Z

    def check_fragmentation(self, reconstruct: bool = False, debug: int=0):
        if self._subject.type == "cell": 
            assert hasattr(self,"cellvec")
            assert hasattr(self._subject,"refmoleclist")
        else:
            self.fragmented = False
            return self.fragmented

        if not hasattr(self,"moleclist"): self.get_moleclist()
        self.fragmented = False

        # First comparison with current moleclist
        for mol in self.moleclist:
            found = False
            for rmol in self._subject.refmoleclist:
                issame = compare_species(mol, rmol)  
                if issame: found = True
            if not found: self.fragmented = True
        # If there are fragments and user wants reconstruction, it tries to reconstruct and checks the new moleclist
        if self.fragmented and reconstruct:
            new_moleclist = self.reconstruct(debug=debug)
            self.fragmented = False
            for mol in new_moleclist:
                found = False
                for rmol in self._subject.refmoleclist:
                    issame = compare_species(mol, rmol)  
                    if issame: found = True
                if not found: self.fragmented = True
            if not self.fragmented: 
                self.moleclist = new_moleclist
                self.set_geometry_from_moleclist()

        return self.fragmented

    def set_spin_config(self, spin_config):
        self.spin_config = spin_config

    def set_atoms(self):
        if not hasattr(self,"moleclist"): self.get_moleclist()
        self.atoms = []
        indices = []
        for mol in self.moleclist:
            for idx, at in enumerate(mol.atoms):
                self.atoms.append(at)
                indices.append(mol.atlist[idx])
        self.atoms = [x for _, x in sorted(zip(indices, self.atoms), key=lambda pair: pair[0])]

    def set_energy(self, energy, units, overwrite: bool=True):
        self.add_result(data("energy",energy,units,"state.set_energy"), overwrite=overwrite)
        #self.energy     = energy

    def find_computation(self, keyword: str='', step: int=1, run_number: int=1, debug: int=0):
        for idx, comp in enumerate(self.computations):
            if comp.keyword == keyword and comp.step == step and comp.run_number == run_number: this_comp = comp; return True, this_comp
        return False, None

    def add_computation(self, computation: object):
        found, comp = self.find_computation(computation.keyword, computation.step, computation.run_number)
        if not found: 
            self.computations.append(computation)
            computation.add_state(self)

    def reconstruct(self, debug: int=0):
        assert hasattr(self,"cellvec")
        if not hasattr(self._subject,"refmoleclist"): print("CLASS STATE.RECONSTRUCT: _subject does not have refmoleclist"); return None
        from Scope.Other import HiddenPrints
        if debug > 0: print("CLASS_STATE.RECONSTRUCT: reconstructing cell of state", self)
        with HiddenPrints():
            finished = False
            if hasattr(self._subject,"type") and hasattr(self,"cellvec"):
                if not hasattr(self,"moleclist"): self.get_moleclist() 
                if self._subject.type.lower() == "cell":
                    import itertools
                    blocklist = self.moleclist.copy()
                    refmoleclist = self._subject.refmoleclist.copy()
                    cov_factor = refmoleclist[0].cov_factor
                    metal_factor = refmoleclist[0].metal_factor
                    moleclist, fraglist, Hlist = classify_fragments(blocklist, refmoleclist, debug=debug) 
                    if len(fraglist) > 0 or len(Hlist) > 0: 
                        moleclist, finalmols, Warning = fragments_reconstruct(moleclist,fraglist,Hlist,refmoleclist,self.cellvec,cov_factor,metal_factor, debug=debug)
                        moleclist.extend(finalmols)
                        self.moleclist = moleclist
                        self.set_geometry_from_moleclist()
                        finished = True
     
                else: print("WARNING: reconstruct state, _subject is not a cell. I will not reconstruct"); return None
            else: print("WARNING: reconstruct state, _subject does not have 'type' or 'cellvec' variables"); return None
        if debug > 0 and finished: print("CLASS STATE.RECONSTRUCT: state reconstructed to", self)
        return self.moleclist

####################
    def get_thermal_data(self, Trange: range=range(10,501,1), Helec=None, Selec=None, Hvib=None, Svib=None, Gtot=None, overwrite: bool=False, debug: int=0):
        from Scope.Thermal_Corrections import get_Selec, get_Hvib, get_Svib, get_Gibbs

        if not hasattr(self,"ncomplex"): self.get_ncomplex(debug=debug)

        if Hvib is None and Svib is None:
            assert hasattr(self,"isminimum"), f"I can't compute thermal data on this state. Missing VNMs"
    
        ############## Helec ##############
        if Helec is None:   ### One can provide specific values for Helec, Selec, Hvib, Svib and Gtot 
            if overwrite or not "Helec" in self.results.keys():
                self.add_result(data("Helec",self.results["energy"].value/self.ncomplex,self.results["energy"].units,"state.get_thermal_data"), overwrite=overwrite)
        else: 
            if isinstance(Helec, data):
                if overwrite or not "Helec" in self.results.keys():
                    self.add_result(data("Helec",Helec.value,Helec.units,"enforced in get_thermal_data"), overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Helec. It must be a DATA-class object")

        ############## Selec ##############
        if Selec is None:
            if overwrite or not "Selec" in self.results.keys():
                self.add_result(get_Selec(self._subject.spin, outunits='au', nmol=self.ncomplex), overwrite=overwrite)
        else: 
            if isinstance(Selec, data):
                if overwrite or not "Selec" in self.results.keys():
                    self.add_result(data("Selec",Selec.value,Selec.units,"enforced in get_thermal_data"), overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Selec. It must be a DATA-class object")

        ############## Hvib ##############
        if Hvib is None:
            if overwrite or not "Hvib" in self.results.keys():
                Hvib = collection("Hvib")
                for temp in Trange:
                    Hvib.add_data(get_Hvib(np.abs(self.freqs_cm), temp, freq_units='cm', outunits='au', nmol=self.ncomplex))
                self.add_result(Hvib, overwrite=overwrite)
        else: 
            if isinstance(Hvib, collection):
                if overwrite or not "Hvib" in self.results.keys():
                    self.add_result(Hvib, overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Hvib. It must be a COLLECTION-class object")

        ############## Svib ##############
        if Svib is None:
            if overwrite or not "Svib" in self.results.keys():
                Svib = collection("Svib")
                for temp in Trange:
                    Svib.add_data(get_Svib(np.abs(self.freqs_cm), temp, freq_units='cm', outunits='au', nmol=self.ncomplex))
                self.add_result(Svib, overwrite=overwrite)
        else: 
            if isinstance(Svib, collection):
                if overwrite or not "Svib" in self.results.keys():
                    self.add_result(Svib, overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Svib. It must be a COLLECTION-class object")

        ############## Gtot ##############
        if Gtot is None:
            if overwrite or not "Gtot" in self.results.keys():
                Gtot = collection("Gtot")
                for temp in Trange:
                    # Retrieve data (not value)
                    Helec = self.results["Helec"]
                    Selec = self.results["Selec"]
                    Hvib_i = Hvib.find_value_with_property("temp", temp)
                    Svib_i = Svib.find_value_with_property("temp", temp)
                    assert Helec.units == Selec.units == Hvib_i.units == Svib_i.units, f"{Helec.units=}, {Selec.units=}, {Hvib_i.units=}, {Svib_i.units=}"
                    key = "Gtot"
                    value = get_Gibbs(Helec.value, Hvib_i.value, Selec.value, Svib_i.value, temp)
                    units = Helec.units
                    function = "state.get_thermal_data"
                    new_data = data(key, value, units, function)
                    new_data.add_property("temp", temp, overwrite=overwrite)
                    Gtot.add_data(new_data)
                self.add_result(Gtot, overwrite=overwrite)
        else: 
            if isinstance(Gtot, collection):
                if overwrite or not "Gtot" in self.results.keys():
                    self.add_result(Gtot, overwrite=overwrite)
            else:
                print("Get_Thermal_Data: wrong type of data provided when enforcing Gtot. It must be a COLLECTION")

    def __repr__(self) -> None:
        to_print  = f'---------------------------------------------------\n'
        to_print +=  '   STATE                                           \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Name                  = {self.name}\n'
        if hasattr(self,"labels"):         to_print += f' Labels                = {self.labels[0]}...\n'
        if hasattr(self,"coord"):          to_print += f' Coord                 = {self.coord[0]}...\n'
        if hasattr(self,"ncomplex"):       to_print += f' Number of Complexes   = {self.ncomplex}\n' 
        if hasattr(self,"isminimum"):      to_print += f' Is Minimum            = {self.isminimum}\n'
        if hasattr(self,"almost_minimum"): to_print += f' Almost a Minimum      = {self.almost_minimum}\n'
        if hasattr(self,"freq_cm"):        to_print += f' First Frequency (cm-1)= {self.freq_cm[0]}...\n'
        if hasattr(self,"moleclist"):  
            to_print += f' # Molecules:          = {len(self.moleclist)}\n'
            to_print += f' With Formulae:                               \n'
            for idx, m in enumerate(self.moleclist):
                to_print += f'    {idx}: {m.formula} \n'
        return to_print

#########################################################################
## Tools associated with states. Normally, these would be class (i.e. molecule) functions...
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
                if debug >= 1: print(f"FIND STATE: state {search_name} found")
                return True, sta
        if not found: 
            if debug >= 1: print(f"FIND STATE: state {search_name} not found")
            return False, None




## OLD ATTEMPT AT UPDATE

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

#    def _update(self, other, debug: int=0):
#        if debug > 0: print("Updating State", self.name)
#        if debug > 0: print("SELF", dir(self))
#        if debug > 0: print("OTHER", dir(other))
#        if not isinstance(other, type(self)): 
#            if debug > 0: print("_update STATE: other is not of the same type than self")
#            return self
#        for d in dir(other):
#            if debug > 0: print("_update STATE: checking instance", d)
#            if d in dir(other):
#                if '_' not in d and not callable(getattr(other,d)) and d != "type" and d != "name" and d != "results" and d != "computations":
#                    if debug > 0: print("_update STATE: updating instance", d)
#                    at1 = getattr(other,d)
#                    try:      attr = literal_eval(at1)
#                    except:   attr = value
#                    setattr(self, d, attr)
#                elif d == "results":
#                    for r in other.results:
#                        if debug > 0: print("_update STATE: updating result", r)
#                        self.add_result(r, overwrite=True)
#                elif d == "computations": 
#                    for c in other.computations:
#                        if c not in self.computations: 
#                            if debug > 0: print("_update STATE: updating computations, adding", c)
#                            self.computations.append(c)
#        return self

