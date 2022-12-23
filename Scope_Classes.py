import sys
import copy
from copy import deepcopy
import os

from Scope.Parse_Cif import get_cif_diffraction_data, get_cif_authors, get_cif_journal 
from Scope.Parse_General import search_string, read_lines_file 
from Scope.Geom_SCO_V1 import geom_sco_from_xyz, guess_spin_state
#from Scope.Parse_G16_outputs.py import G16_time_to_sec

from cell2mol.tmcharge_common import labels2formula
from cell2mol.elementdata import ElementData
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

################
##### JOBS #####
################
class recipe(object):
    def __init__(self, obj: object, path: str, keyword: str) -> None:
        self.object = obj
        self.path = path
        self.keyword = keyword
        self.jobs = []

class job(object):
    def __init__(self, input_path: str, output_path: str, subfile_path: str, software: str, code: str='') -> None:
        self.input_path = input_path
        self.output_path = output_path
        self.subfile_path = subfile_path
        self.software = software
        self.code = code
        self.requisites = []

    def add_requisite(self, requisite) -> None:
        self.requisites.append(requisite)

    def add_submission_init(self, nprocs: str='Unk', queue: str='Unk', issubmitted: bool=False) -> None:
        self.nprocs = nprocs
        self.queue = queue
        self.issubmitted = issubmitted

############################
##### SYSTEM & CRYSTAL #####
############################
class sco_system(object):
    def __init__(self, refcode: str, cell2mol_path: str, scope_path: str) -> None:
        self.type = "sco_system" 
        self.version = "0.2" 
        self.refcode = refcode
        self.refcode_wo_digits = ''.join([i for i in refcode if not i.isdigit()])
        self.cell2mol_path = cell2mol_path
        self.scope_path = scope_path
        self.list_of_crystals = []
        self.list_of_recipes = []
    
    ##########
    def set_reference_molecs(self, debug: int=0):
        self.hasLS = False
        self.hasHS = False

        pool = []
        ##### All these below could go somewhere else #####
        ### Colects reference molecules from all crystals
        for crys in self.list_of_crystals:
            if len(crys.list_of_molecules) == 0: crys.get_FeN6_molecules(debug=debug)
            for mol in crys.list_of_molecules:
                pool.append(mol) 
            
        ### Selects reference molecules
        if debug > 0: print(f"{len(pool)} tmcs in pool")
        if len(pool) > 0:
            for idx, mol in enumerate(pool):
                if mol.scope_guess_spin == 'LS' and not self.hasLS:
                    self.hasLS = True
                    self.LSid = idx
                    self.LSref = deepcopy(mol)
                elif mol.scope_guess_spin == 'HS' and not self.hasHS:
                    self.hasHS = True
                    self.HSid = idx
                    self.HSref = deepcopy(mol)
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

    ##########
    def set_reference_crystals(self, debug: int=0):
        self.has_HS_ref_crys = False
        self.has_LS_ref_crys = False

        #### Selects reference crystals:
        HS_ref_crys_temp = 1000
        LS_ref_crys_temp = 1000
        for idx, crys in enumerate(self.list_of_crystals):
  
            if crys.phase == '': crys.get_spin_and_phase_data(debug=debug)

            ### This try/except block is to skip cases in which the diffraction temperature is not known
            try: crys.diff_temp = float(crys.diff_temp)
            except: continue

            if crys.phase == "HS" and crys.diff_temp < HS_ref_crys_temp: 
                self.HS_ref_crys = deepcopy(crys)
                self.HS_ref_crys_id = idx
                self.HS_ref_crys_temp = crys.diff_temp                
                self.hasHScrys = True
            elif crys.phase == "LS" and crys.diff_temp < LS_ref_crys_temp: 
                self.LS_ref_crys = deepcopy(crys)
                self.LS_ref_crys_id = idx
                self.LS_ref_crys_temp = crys.diff_temp                
                self.has_LS_ref_crys = True
            elif crys.phase == "IS":
                pass

        if not self.has_HS_ref_crys:
            self.HS_ref_crys = deepcopy(self.list_of_crystals[0])
            self.HS_ref_crys_id = 0
            self.HS_ref_crys.phase == 'HS'
        if not self.has_LS_ref_crys:
            self.LS_ref_crys = deepcopy(self.list_of_crystals[0])
            self.LS_ref_crys_id = 0
            self.LS_ref_crys.phase == 'LS'

    ##########
    def get_recipe(self, obj, keyword, debug: int=0):
        ## obj is the actual object that will be computed. For instance, a molecule, or the cell.
        recipe_exists = False
        for idx, rec in enumerate(self.list_of_recipes):
            if rec.keyword == keyword and rec.object == obj and os.path.isdir(rec.path):
                recipe_exists = True
                current_recipe = rec

        ## If such entry doesn't exist. Then we create it
        if not recipe_exists:
            current_recipe = recipe(obj, self.scope_path+'/'+keyword, keyword)
            self.list_of_recipes.append(current_recipe)
            if not os.path.isdir(self.scope_path+'/'+keyword): os.makedirs(self.scope_path+'/'+keyword)
        return recipe_exists, current_recipe

class crystal(object):
    def __init__(self, refcode: str, name: str, cell2mol_path: str, cell: object) -> None:
        self.type = "crystal" 
        self.version = "0.2" 
        self.refcode = refcode
        self.name = name
        self.cell = cell 
        self.cell2mol_path = cell2mol_path
        self.list_of_molecules = []

    def read_cif_data(self, cifpath: str) -> None:
        self.diff_temp = get_cif_diffraction_data(cifpath) 
        self.authors = get_cif_authors(cifpath)
        self.journal_year =   get_cif_journal(cifpath)[0]
        self.journal_name =   get_cif_journal(cifpath)[1] 
        self.journal_volume = get_cif_journal(cifpath)[2]
        self.journal_page =   get_cif_journal(cifpath)[3]

    def add_iso_calcs_path(self, iso_calcs_path: str='Unk'):
        self.iso_calcs_path = iso_calcs_path

    def get_spin_and_phase_data(self, debug: int=0):
        self.guess_spins = [] 
        self.phase = ''
        self.HSmolarfraction = 0.0
        for idx, mol in enumerate(self.list_of_molecules):
            self.guess_spins.append(mol.scope_guess_spin)

        if all(sp == "HS" for sp in self.guess_spins): 
            self.phase = "HS"
            self.HSmolarfraction = 1.0
        elif all(sp == "LS" for sp in self.guess_spins): 
            self.phase = "LS"
            self.HSmolarfraction = 0.0
        elif "HS" in self.guess_spins and "LS" in self.guess_spins: 
            self.phase = "IS"
            self.HSmolarfraction = float(self.guess_spins.count("HS")/len(self.guess_spins))

    def get_FeN6_molecules(self, debug: int=0):
        for idx, mol in enumerate(self.cell.moleclist):
            if mol.type == "Complex":
                keepit = False
                for met in mol.metalist:
                    if hasattr(met, "coord_sphere"):
                        formula = labels2formula(met.coord_sphere)
                        if formula == "N6": keepit = True
                    else: print("No coord_sphere variable in metal object")
                if keepit:
                    ox_state = mol.metalist[0].totcharge
                    mol.scope_FeNdist, mol.scope_FeNangle = geom_sco_from_xyz(mol.labels, mol.coord)
                    mol.scope_guess_spin = guess_spin_state(int(ox_state), mol.scope_FeNdist[0])
                    self.list_of_molecules.append(mol)
       

