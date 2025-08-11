import sys
import os
import numpy as np

from Scope.Parse_General import search_string, read_lines_file 
from Scope.Unit_cell_tools import get_unit_cell_volume
from Scope import Constants

from Scope.Adapted_from_cell2mol import labels2formula
from Scope.Elementdata import ElementData
elemdatabase = ElementData()

######################
##### QC Objects #####
######################
class VNM(object):
    def __init__(self, index: int, freq: float, red_mass: float=1.0, force_cnt: float=0.0, IR_int: float=0.0, sym: str='A'):
        self.index = index 
        self.freq_cm = freq                            ## In cm-1
        self.freq = freq*Constants.cm2har              ## In atomic units 
        self.red_mass = red_mass                       ## In AMU     as in Gaussian
        self.force_cnt = force_cnt                     ## In mDyne/A as in Gaussian
        self.IR_int = IR_int                           ## In KM/Mole as in Gaussian
        self.sym = sym
        self.haseigenvec = False

    def eigenvec(self, atomidxs: list, atnums: list, xs: list, ys: list, zs: list):
        self.atomidxs = atomidxs
        self.atnums = atnums
        self.labels = [elemdatabase.elementsym[atnum] for atnum in atnums]
        self.masses = [elemdatabase.elementweight[l] for l in self.labels]
        self.xs = xs
        self.ys = ys
        self.zs = zs
        self.eigenvec_format1 = np.column_stack([xs, ys, zs])
        self.eigenvec_format2 = np.column_stack([xs, ys, zs]).reshape(-1)
        self.haseigenvec = True

    def __repr__(self) -> None:
        to_print  = f'-----------------------------\n'
        to_print +=  '   Vibrational Normal Mode   \n'
        to_print += f'-----------------------------\n'
        to_print += f' Index                  = {self.index}\n'
        to_print += f' Freq (cm-1)            = {self.freq_cm}\n'
        to_print += f' IR Intensity (KM/Mole) = {self.IR_int}\n'
        to_print += f' Reduced Mass (AMU)     = {self.red_mass}\n'
        to_print += f' Has eigenvector        = {self.haseigenvec}\n'
        return to_print

class orbital_set(object):
    def __init__(self, name: str) -> None:
        self.name = name
        self.orbitals = []
    
    def get_eigenvalues(self, eigen_search_string: str, lines: list):
        self.eigen_search_string = eigen_search_string
        self.eigen_search_line_numbers, self.eigen_search_line_found = search_string(eigen_search_string, lines, typ='all')
        self.eigen_search_lines = []
        for l in self.eigen_search_line_numbers:
            self.eigen_search_lines.append(lines[l])        
        self.eigenvalues = []
        
    def add_orbital(self, orb):
        self.orbitals.append(orb)

class orbital(object):
    def __init__(self, index: int, occupation: int, energy: float) -> None:
        self.name = index
        self.occupation = occupation
        self.energy = energy

#class periodic_xyz(object):
#    def __init__(self, name: str, index: int, labels: list, pos: list, cellparam: list) -> None:
#        self.name = name
#        self.index = index
#        self.labels = labels  
#        self.pos = pos               
#        self.cellparam = cellparam  # in Bohr
#        self.volume = get_unit_cell_volume(cellparam[0], cellparam[1], cellparam[2], cellparam[3], cellparam[4], cellparam[5]) # in Bohr3
#        self.natoms = len(labels)
#        self.formula = labels2formula(labels) 
#        self.moleclist = []
#
#    def add_cellvec(self, cellvec, celldim):
#        self.cellvec = cellvec
#        self.celldim = celldim
