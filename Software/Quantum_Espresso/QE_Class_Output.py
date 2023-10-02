import Scope.Constants
from Scope.Parse_General import search_string, read_lines_file
from Scope.Software.Quantum_Espresso.Parse_QE_outputs import *

class qe_output(object):
    def __init__(self, lines: list, computation: object=None):
        self._computation   = computation
        self.lines          = lines
        
    def clear_lines(self):
        self.lines          = []
        
    def read_lines(self):
        self.lines          = read_lines_file(self._computation.out_path)
        
    def get_status(self):
        keys=["JOB DONE", "End final coordinates", "Maximum CPU time exceeded", "convergence NOT achieved after"]
        linenum = []
        bools = []        
        ## Search Key Elements of Output
        for sx, string in enumerate(keys):
            search, found = search_string(string, self.lines, typ="last")
            linenum.append(search)
            bools.append(found)        
        ## Set Status based on Findings
        if   bools[0] and bools[1] and not bools[2]:                      self.status = "worked"           ### Everything worked
        elif bools[0] and bools[3] and not bools[1] and not bools[2]:     self.status = "not converged"    ### Job ended but did not converge
        elif bools[2]:                                                    self.status = "stopped"          ### Job Stopped due to time limits
        elif all(not f for f in bools):                                   self.status = "aborted"          ### Nothing, probably a job terminated abruptly
        else:                                                             self.status = "unknown"
        return self.status
    
###############
### BLOCKS ###
###############
    def get_scf_blocks(self, debug: int=0):
        self.scf_blocks = []        
        ## Strings to Search
        start_SCF, found = search_string("Self-consistent Calculation", self.lines)
        end_SCF, found   = search_string("End of self-consistent calculation", self.lines)        
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
        for idx, b in enumerate(self.scf_blocks[:-1]):
            self.opt_blocks.append([self.scf_blocks[idx][0], self.scf_blocks[idx+1][0]])
        return self.opt_blocks
    
##############
### GEOMS ###
##############
    def get_all_geometries(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_labels = []
        self.all_coords = []
        for idx, b in enumerate(self.opt_blocks):
            step_labels, step_coords       = parse_geometry_from_step(self.lines[b[0]+1:b[1]+1])    
            if idx > 0: assert step_labels[0] == self.all_labels[0][0], f"There is an issue with the parsed labels"
            self.all_labels.append(step_labels)
            self.all_coords.append(step_coords)
        return self.all_labels, self.all_coords

    def get_last_geometry(self, debug: int=0):
        if hasattr(self,"all_labels"): 
            self.last_labels = self.all_labels[-1]
            self.last_coords = self.all_coords[-1]
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0: return None, None
            init_line = self.opt_blocks[-1][0]+1
            last_line = self.opt_blocks[-1][1]+1
            self.last_labels, self.last_coords  = parse_geometry_from_step(self.lines[init_line:last_line])
        return self.last_labels, self.last_coords

###################
### CELL PARAMS ###
###################
    def get_all_cell_param(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_cell_vec = []
        for idx, b in enumerate(self.opt_blocks):
            step_cell_vec  = parse_cell_vectors(self.lines[b[0]+1:b[1]+1])
            self.all_cell_vec.append(step_cell_vec)
        return self.all_cell_vec

    def get_last_cell_param(self, debug: int=0):
        if hasattr(self,"all_cell_vec"): 
            self.last_cell_vec = self.all_cell_vec[-1] 
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0: return None
            init_line = self.opt_blocks[-1][0]+1
            last_line = self.opt_blocks[-1][1]+1
            self.last_cell_vec  = parse_cell_vectors(self.lines[init_line:last_line])
        return self.last_cell_vec

##############
### FORCES ###
##############
    def get_all_forces(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_at_forces = []
        self.all_tot_forces = []
        for idx, b in enumerate(self.opt_blocks):
            step_at_forces, step_tot_force = parse_forces_from_step(self.lines[b[0]+1:b[1]+1])
            assert len(step_at_forces) == len(step_labels), f"{len(step_forces)} {len(step_labels)}"
            self.all_at_forces.append(step_at_forces)
            self.all_tot_forces.append(step_tot_force)
        return self.all_at_forces, self.all_tot_forces

    def get_last_forces(self, debug: int=0):
        if hasattr(self,"all_at_forces"): 
            self.last_at_forces = self.all_at_forces[-1] 
            self.last_tot_force = self.all_tot_forces[-1] 
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0: return None, None
            init_line = self.opt_blocks[-1][0]+1
            last_line = self.opt_blocks[-1][1]+1
            self.last_at_forces, self.last_tot_force = parse_forces_from_step(self.lines[init_line:last_line])
        return self.last_at_forces, self.last_tot_force
            
################
### ENERGIES ###
################
    def get_all_energies(self, debug: int=0):
        from Scope.Classes_Data import data, collection
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_energies = [] 
        for b in self.opt_blocks:
            val      = parse_final_energy(self.lines[b[0]+1:b[1]+1])
            self.all_energies.append(val)
        return self.all_energies

    def get_last_energy(self, debug: int=0):
        if hasattr(self,"all_energies"): 
            self.last_energy = self.all_energies[-1]
        else:
            from Scope.Classes_Data import data, collection
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0: return None, None
            init_line = self.opt_blocks[-1][0]+1
            last_line = self.opt_blocks[-1][1]+1
            self.last_energy = parse_final_energy(self.lines[init_line:last_line])
        return self.last_energy
    
############
### TIME ###
############
    def get_elapsed_time(self, debug: int=0):
        line_elapsed, found_elapsed = search_string("PWSCF        :", self.lines, typ='last')
        if found_elapsed:
            eline = self.lines[line_elapsed]
            self.elapsed_time = QE_elapsed_time(eline)
        else: self.elapsed_time = float(0)
        return self.elapsed_time        
    
    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print += f'   FORMATTED REPRESENTATION OF OUTPUT CLASS        \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f'# Lines           = {len(self.lines)}\n'
        if hasattr(self,"status"):       to_print += f'Status            = {self.status}\n'
        if hasattr(self,"scf_blocks"):   to_print += f'#SCF Blocks       = {len(self.scf_blocks)}\n'
        if hasattr(self,"opt_blocks"):   to_print += f'#OPT Blocks       = {len(self.opt_blocks)}\n'
        if hasattr(self,"all_coords"):   to_print += f'#Coordinates      = {len(self.all_coords)}\n'
        if hasattr(self,"all_energies"): to_print += f'#Energies         = {len(self.all_energies)}\n'
        if hasattr(self,"all_cell_vec"): to_print += f'#Cell Vectors     = {len(self.all_cell_vec)}\n'
        if hasattr(self,"elapsed_time"): to_print += f'Elapsed Time      = {self.elapsed_time}\n'
        return to_print
