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
        self.index        = index 
        self.freq_cm      = freq                     ## In cm-1
        self.freq         = freq*Constants.cm2har    ## In atomic units 
        self.red_mass     = red_mass                 ## In AMU     as in Gaussian
        self.force_cnt    = force_cnt                ## In mDyne/A as in Gaussian
        self.IR_int       = IR_int                   ## In KM/Mole as in Gaussian
        self.sym          = sym
        self.has_mode     = False

    def set_mode(self, atomidxs: list, atnums: list, xs: list, ys: list, zs: list):
        self.atomidxs     = atomidxs
        self.atnums       = atnums
        self.labels       = [elemdatabase.elementsym[atnum] for atnum in atnums]
        self.masses       = [elemdatabase.elementweight[l] for l in self.labels]
        self.mode         = np.column_stack([xs, ys, zs])
        self.mode_format2 = np.column_stack([xs, ys, zs]).reshape(-1)
        self.has_mode     = True

    def mass_weight_mode(self):
        if not self.has_mode: return None
        mw = np.repeat(np.sqrt(self.masses), 3)
        self.mode_mw      = self.mode * mw[np.newaxis, :]         ## Mass_Weighted Version of the Mode
        return self.mode_mw
    
    def write_dyn(self, initial_coord: list, amplitude: int=10, outfolder: str='./', labels: None=list):
        ## Writes a file with a trajectory representing the displacement of the VNM
        from Scope.Read_Write import write_xyz
        filename: str="dyn_vnm_"+str(self.index)+".xyz"
        if outfolder[-1] != '/': outfolder += '/'
        initial_coord = np.array(initial_coord)
        for f in range(-amplitude,amplitude,1):
            labels = []
            coords = []
            for idx in range(natoms):
                vector = self.mode[idx]
                coord  = initial_coord[idx] + vector*f*0.1
                if labels is None: label = elemdatabase.elementsym[vnm.atnums[idx]]
                else:              label = self.labels[idx]
                labels.append(label)
                coords.append(coord)
            write_xyz(outfolder+filename, labels, coords, append=True) 

    def __repr__(self) -> None:
        to_print  = f'-----------------------------\n'
        to_print +=  '   Vibrational Normal Mode   \n'
        to_print += f'-----------------------------\n'
        to_print += f' Index                  = {self.index}\n'
        to_print += f' Freq (cm-1)            = {self.freq_cm}\n'
        to_print += f' IR Intensity (KM/Mole) = {self.IR_int}\n'
        to_print += f' Reduced Mass (AMU)     = {self.red_mass}\n'
        to_print += f' Has Mode               = {self.has_mode}\n'
        return to_print

######
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
