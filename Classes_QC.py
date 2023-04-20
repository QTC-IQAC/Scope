import sys
import os

from Test_V3.Parse_General import search_string, read_lines_file 
from Test_V3.Unit_cell_tools import get_unit_cell_volume
from Test_V3 import Constants

from Test_V3.Adapted_from_cell2mol import labels2formula
from Test_V3.Elementdata import ElementData
elemdatabase = ElementData()

######################
##### QC Objects #####
######################
class VNM(object):
    def __init__(self, index: int, freq: float, red_mass: float, force_cnt: float, IR_int: float, sym: str='A'):
        self.index = index 
        self.freq_cm = freq                            ## It must be in cm-1
        self.freq = freq*Constants.cm2har              ## It must be in atomic units 
        self.red_mass = red_mass
        self.force_cnt = force_cnt
        self.IR_int = IR_int
        self.sym = sym

    def eigenvec(self, atomidxs: list, atnums: list, xs: list, ys: list, zs: list):
        self.atomidxs = atomidxs
        self.atnums = atnums
        self.labels = [elemdatabase.elementsym[atnum] for atnum in atnums]
        self.xs = xs
        self.ys = ys
        self.zs = zs

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

class periodic_xyz(object):
    def __init__(self, name: str, index: int, labels: list, pos: list, cellparam: list) -> None:
        self.name = name
        self.index = index
        self.labels = labels  
        self.pos = pos               
        self.cellparam = cellparam  # in Bohr
        self.volume = get_unit_cell_volume(cellparam[0], cellparam[1], cellparam[2], cellparam[3], cellparam[4], cellparam[5]) # in Bohr3
        self.natoms = len(labels)
        self.formula = labels2formula(labels) 
        self.moleclist = []

    def add_cellvec(self, cellvec, celldim):
        self.cellvec = cellvec
        self.celldim = celldim
