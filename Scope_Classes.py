import sys
import copy
from copy import deepcopy

scopepath = '/home/g4vela/SCOPE/Database_SCO/Scripts/Scope'
sys.path.append(scopepath)

from cell2mol.tmcharge_common import labels2formula
from cell2mol.elementdata import ElementData
elemdatabase = ElementData()

from Parse_General import search_string, read_lines_file
from Parse_Cif import *
from unit_cell_tools import get_unit_cell_volume #, cellvec_2_cellparam 
from Geom_SCO_V1 import *

import Constants

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

############################
##### SYSTEM & CRYSTAL #####
############################
class sco_system(object):
    def __init__(self, refcode: str, cell2mol_path: str) -> None:
        self.type = "sco_system" 
        self.version = "0.1" 
        self.refcode = refcode
        self.refcode_wo_digits = ''.join([i for i in refcode if not i.isdigit()])
        self.cell2mol_path = cell2mol_path
        self.list_of_crystals = []

    def add_iso_calcs_path(self, iso_calcs_path: str='Unk'):
        self.iso_calcs_path = iso_calcs_path

    ####
    def set_reference_molecs(self, debug: int=0):

        self.list_of_refmolecs = []
        self.hasLS = False
        self.hasHS = False
        for crys in self.list_of_crystals:
            done = False
            while not done:
                for idx, mol in enumerate(crys.cell.moleclist):
                    if debug > 0: print(len(crys.cell.moleclist), "molecules in crystal", crys.name, "of", self.refcode)
                    if mol.type == "Complex":
    
                        # Checks it is an Fe-N6
                        keepit = False
                        for met in mol.metalist:
                            if hasattr(met, "coord_sphere"):
                                formula = labels2formula(met.coord_sphere)
                                if formula == "N6": keepit = True
                            else: print("No coord_sphere variable in metal object")
    
                        # If so...
                        if keepit:
                            ox_state = mol.metalist[0].totcharge
                            dist, angle = geom_sco_from_xyz(mol.labels, mol.coord)
                            guess_spin = guess_spin_state(int(ox_state), dist[0])
                            crys.add_spin_data(guess_spin)
    
                            mol.scope_FeNdist = dist
                            mol.scope_FeNangle = angle
                            mol.scope_guess_spin = guess_spin
                            #mol.scope_cell2mol_path = crys.cell2mol_path+'/'+gmolname+'.sco'
                            #mol.scope_job_story = job_story(mol.scope_gmol_path)
    
                            self.list_of_refmolecs.append(mol)
                            if debug > 0: print(f"appended reference molecule")
                            done = True

        pool = self.list_of_refmolecs.copy()

        if debug > 0: print(f"{len(pool)} tmcs in pool")
        if len(pool) > 0:
            for idx, tmc in enumerate(pool):
                if tmc.scope_guess_spin == 'LS' and not self.hasLS:
                    self.hasLS = True
                    self.LSid = idx
                    self.LSref = deepcopy(tmc)
                elif tmc.scope_guess_spin == 'HS' and not self.hasHS:
                    self.hasHS = True
                    self.HSid = idx
                    self.HSref = deepcopy(tmc)
            ### If it hasn't found any molecule that can be classified as HS and LS... then takes anything
            if not self.hasLS:
                self.LSid = 0 
                self.LSref = deepcopy(pool[0])
                self.LSref.scope_guess_spin == 'LS'
            if not self.hasHS:
                self.HSid = 0 
                self.HSref = deepcopy(pool[0])
                self.HSref.scope_guess_spin == 'HS'
        else: print("Empty pool of reference molecules")

class crystal(object):
    def __init__(self, refcode: str, name: str, cell2mol_path: str, cell: object) -> None:
        self.type = "crystal" 
        self.version = "0.1" 
        self.refcode = refcode
        self.name = name
        self.cell = cell 
        self.cell2mol_path = cell2mol_path

    def read_cif_data(self, cifpath: str) -> None:
        self.diff_temp = get_cif_diffraction_data(cifpath) 
        self.authors = get_cif_authors(cifpath)
        self.journal_year =   get_cif_journal(cifpath)[0]
        self.journal_name =   get_cif_journal(cifpath)[1] 
        self.journal_volume = get_cif_journal(cifpath)[2]
        self.journal_page =   get_cif_journal(cifpath)[3]

    def add_spin_data(self, guess_spin_state: str='Unk'):
        self.guess_spin_state = guess_spin_state

    def add_iso_calcs_path(self, iso_calcs_path: str='Unk'):
        self.iso_calcs_path = iso_calcs_path

################
##### JOBS #####
################
class job_story(object):
    def __init__(self, gmol_path: str) -> None:
        self.gmol_path = gmol_path
        self.job_list = []

    def add_iso_calcs_path(self, iso_calcs_path: str='Unk'):
        self.iso_calcs_path = iso_calcs_path

class job_class(object):
    def __init__(self, input_path: str, output_path: str, subfile_path: str, code: str, summary: str='') -> None:
        self.input_path = input_path
        self.output_path = output_path
        self.subfile_path = subfile_path
        self.code = code
        self.summary = summary

    def add_submission_init(self, nprocs: str='Unk', queue: str='Unk', issubmitted: bool=False) -> None:
        self.nprocs = nprocs
        self.queue = queue
        self.issubmitted = issubmitted

    def add_submission_end(self, isfinished: bool=False, isgood: bool=False, isregistered: bool=False, elapsed_time: float=float(0)) -> None:
        self.isfinished = isfinished 
        self.isgood = isgood
        self.isregistered = isregistered
        self.elapsed_time = elapsed_time 
