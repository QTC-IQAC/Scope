import Scope.Constants
from Scope.Parse_General import search_string, read_lines_file
from Scope.Software.Gaussian.Parse_G16_outputs import *

class g16_output(object):
    def __init__(self, lines: list, computation: object=None):
        self._computation   = computation
        self.lines          = lines
        
    def clear_lines(self):
        self.lines          = []
        
    def read_lines(self):
        if hasattr(self,"_computation"): self.lines = read_lines_file(self._computation.out_path)
        
    def get_status_finished(self):
        keys=["Normal termination", "Error termination"]
        bools = []        
        ## Search Key Elements of Output
        for sx, string in enumerate(keys):
            search, found = search_string(string, self.lines, typ="last")
            bools.append(found)        
        ## Set Status based on Findings
        if bools[0] or bools[1]:             self.finished = True 
        else:                                self.finished = False
        return self.finished

    def get_status_good(self):
        keys=["SCF Done"]
        bools = []        
        ## Search Key Elements of Output
        for sx, string in enumerate(keys):
            search, found = search_string(string, self.lines, typ="last")
            bools.append(found)        
        ## Set Status based on Findings
        if bools[0]:                         self.good = True 
        else:                                self.good = False
        return self.good
    
###############
### BLOCKS ###
###############
    def get_scf_blocks(self, debug: int=0):
        self.scf_blocks = []        
        ## Strings to Search
        start_SCF, found = search_string("Cycle   1  Pass 1", self.lines)
        end_SCF, found   = search_string("SCF Done", self.lines)        
        ## Finds Blocks
        last_start = 0; last_end   = 0
        for s in start_SCF:
            found_end = False
            for e in end_SCF:
                if s > last_end and e > last_start and e > s: 
                    found_end = True; last_start = s; last_end = e; self.scf_blocks.append([s,e])
            if not found_end:
                if s == start_SCF[-1]: 
                    if debug > 0: print(f"SCF step starting at line {s} is not finished")
                    self.scf_blocks.append([s,len(self.lines)])                    
        return self.scf_blocks
    
    def get_opt_blocks(self, debug: int=0):
        if not hasattr(self,"scf_blocks"): self.get_scf_blocks()
        self.opt_blocks = []
        if len(self.scf_blocks) > 0: 
            for idx, b in enumerate(self.scf_blocks[:-1]):
                self.opt_blocks.append([self.scf_blocks[idx][0], self.scf_blocks[idx+1][0]])
            self.opt_blocks.append([self.scf_blocks[-1][0],len(self.lines)]) 
        return self.opt_blocks
    
##############
### GEOMS ###
##############
    def get_all_geometries(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_labels = []
        self.all_coords = []
        if len(self.opt_blocks) > 0:
            for idx, b in enumerate(self.opt_blocks):
                step_labels, step_coords       = parse_geometry_from_step(self.lines[b[0]+1:b[1]+1])
                self.all_labels.append(step_labels)
                self.all_coords.append(step_coords)
                if (step_labels is None or step_coords is None) and debug > 0: print("GET_ALL_GEOMETRIES: error parsing geometry between lines", b)
        else:
            if debug > 0: print("GET_ALL_GEOMETRIES: empty list of opt_blocks")
            self.all_labels = None
            self.all_coords = None
        return self.all_labels, self.all_coords

    def get_last_geometry(self, debug: int=0):
        if hasattr(self,"all_labels"):
            if len(self.all_labels) > 0 and len(self.all_coords) > 0:
                for idx in range(len(self.all_labels)-1,-1,-1):
                    if self.all_labels[idx] is not None and self.all_coords is not None:
                        self.last_labels = self.all_labels[idx]
                        self.last_coords = self.all_coords[idx]
                self.last_labels = None
                self.last_coords = None
                return self.last_labels, self.last_coords
            else:
                if debug > 0: print("GET_LAST_GEOMETRY: No geometries in list")
                self.last_labels = None
                self.last_coords = None
                return self.last_labels, self.last_coords
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0:
                if debug > 0: print("GET_LAST_GEOMETRY: No geometries in list")
                self.last_labels = None
                self.last_coords = None
                return self.last_labels, self.last_coords
            for idx in range(len(self.opt_blocks)-1,-1,-1):
                init_line = self.opt_blocks[idx][0]+1
                last_line = self.opt_blocks[idx][1]+1
                tmp1, tmp2 = parse_geometry_from_step(self.lines[init_line:last_line])
                if tmp1 is not None and tmp2 is not None:
                    self.last_labels = tmp1
                    self.last_coords = tmp2
                    return self.last_labels, self.last_coords
            if debug > 0: print("GET_LAST_GEOMETRY: all geometries are None")
            self.last_labels = None
            self.last_coords = None
            return self.last_labels, self.last_coords

##############
### FORCES ###
##############
    def get_all_forces(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_at_forces = []
        if len(self.opt_blocks) > 0:
            for idx, b in enumerate(self.opt_blocks):
                step_at_forces = parse_forces_from_step(self.lines[b[0]+1:b[1]+1])
                self.all_at_forces.append(step_at_forces)
                if step_at_forces is None and debug > 0: print("GET_ALL_FORCES: error parsing forces between lines", b)
        else:
            if debug > 0: print("GET_ALL_FORCES: empty list of opt_blocks")
            self.all_at_forces = None
        return self.all_at_forces

    def get_last_forces(self, debug: int=0):
        if hasattr(self,"all_at_forces"):
            if len(self.all_at_forces) > 0: 
                for idx in range(len(self.all_at_forces)-1,-1,-1):
                    if self.all_at_forces[idx] is not None:
                        self.last_at_forces  = self.all_at_forces[idx]
                self.last_at_forces = None
                return self.last_at_forces
            else:
                if debug > 0: print("GET_LAST_FORCES: No all_forces in list")
                self.last_at_forces = None
                return self.last_at_forces
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0:
                if debug > 0: print("GET_LAST_FORCES: No opt_blocks in list")
                self.last_at_forces = None
                return self.last_at_forces
            for idx in range(len(self.opt_blocks)-1,-1,-1):
                init_line = self.opt_blocks[idx][0]+1
                last_line = self.opt_blocks[idx][1]+1
                tmp1 = parse_forces_from_step(self.lines[init_line:last_line])
                if tmp1 is not None:
                    self.last_at_forces = tmp1
                    return self.last_at_forces
            if debug > 0: print("GET_LAST_FORCES: all forces are None")
            self.last_at_forces = None
            return self.last_at_forces

################
### ENERGIES ###
################
    def get_all_energies(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_energies = []
        if len(self.opt_blocks) > 0:
            for idx, b in enumerate(self.opt_blocks):
                step_energy  = parse_energy_from_step(self.lines[b[0]+1:b[1]+1])
                self.all_energies.append(step_energy)
                if step_energy is None and debug > 0: print("GET_ALL_ENERGIES: error parsing energy between lines", b)
        else:
            if debug > 0: print("GET_ALL_ENERGIES: empty list of opt_blocks")
            self.all_energies = None
        return self.all_energies

    def get_last_energy(self, debug: int=0):
        if hasattr(self,"all_energies"):
            if len(self.all_energies) > 0:
                for e in self.all_energies[-1::-1]:
                    if e is not None:
                        self.last_energy = e; return self.last_energy
                self.last_energy = None
                return self.last_energy
            else:
                if debug > 0: print("GET_LAST_ENERGY: No all_energies in list")
                self.last_energy = None
                return self.last_energy
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0:
                if debug > 0: print("GET_LAST_ENERGY: No opt_blocks in list")
                self.last_energy = None
                return self.last_energy
            for idx in range(len(self.opt_blocks)-1,-1,-1):
                init_line = self.opt_blocks[idx][0]+1
                last_line = self.opt_blocks[idx][1]+1
                tmp = parse_energy_from_step(self.lines[init_line:last_line])
                if tmp is not None:
                    self.last_energy = tmp; return self.last_energy
            if debug > 0: print("GET_LAST_ENERGY: all energies are None")
            self.last_energy = None
            return self.last_energy

#########################
### FREQUENCIES & VNM ###
#########################
    def get_frequencies(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        if len(self.opt_blocks) == 0: return None
        init_line = self.opt_blocks[-1][0]+1
        last_line = self.opt_blocks[-1][1]+1
        self.frequencies = parse_freqs_from_step(self.lines[init_line:last_line])
        return self.frequencies

    def get_vnms(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        if len(self.opt_blocks) == 0: return None
        init_line = self.opt_blocks[-1][0]+1
        last_line = self.opt_blocks[-1][1]+1
        self.vnms = parse_vnms_from_step(self.lines[init_line:last_line])
        return self.vnms

    
############
### TIME ###
############
    def get_elapsed_time(self, debug: int=0):
        line_elapsed, found_elapsed   = search_string("Elapsed time:", self.lines, typ='last')
        if found_elapsed:
            elapsed_time_list = self.lines[line_elapsed].split()[2:]
            self.elapsed_time = G16_time_to_sec(elapsed_time_list)
        else: self.elapsed_time = float(0)
        return self.elapsed_time        
    
    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print += f'   FORMATTED REPRESENTATION OF OUTPUT CLASS        \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f'#Lines            = {len(self.lines)}\n'
        if hasattr(self,"status"):       to_print += f'Status            = {self.status}\n'
        if hasattr(self,"scf_blocks"):   to_print += f'#SCF Blocks       = {len(self.scf_blocks)}\n'
        if hasattr(self,"opt_blocks"):   to_print += f'#OPT Blocks       = {len(self.opt_blocks)}\n'
        if hasattr(self,"all_coords"):   to_print += f'#Coordinates      = {len(self.all_coords)}\n'
        if hasattr(self,"all_energies"): to_print += f'#Energies         = {len(self.all_energies)}\n'
        if hasattr(self,"last_energy"):  to_print += f'Last Energy       = {self.last_energy}\n'
        if hasattr(self,"elapsed_time"): to_print += f'Elapsed Time      = {self.elapsed_time}\n'
        if hasattr(self,"frequencies"):  
            if self.frequencies is not None: to_print += f'1st Frequency     = {self.frequencies[0]}\n'
        to_print += f'---------------------------------------------------\n'
        return to_print
