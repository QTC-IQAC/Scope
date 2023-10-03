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
        if hasattr(self,"_computation"): self.lines = read_lines_file(self._computation.out_path)
        
    def get_output_status(self):
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
        start_SCF, found1  = search_string("Self-consistent Calculation", self.lines)
        end_SCF, found2    = search_string("End of self-consistent calculation", self.lines)        
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

###################
### CELL PARAMS ###
###################
    def get_all_cell_vectors(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_cell_vec = []
        if len(self.opt_blocks) > 0: 
            for idx, b in enumerate(self.opt_blocks):
                step_cell_vec  = parse_cell_vectors(self.lines[b[0]+1:b[1]+1])
                self.all_cell_vec.append(step_cell_vec)
                if step_cell_vec is None and debug > 0: print("GET_ALL_CELL_PARAM: error parsing cell_parameters between lines", b)
        else:
            if debug > 0: print("GET_ALL_CELL_PARAM: empty list of opt_blocks")
            self.all_cell_vec = None
        return self.all_cell_vec

    def get_last_cell_vectors(self, debug: int=0):
        if hasattr(self,"all_cell_vec"):
            if len(self.all_cell_vec) > 0:
                for c in self.all_cell_vec[-1::-1]:
                    if c is not None:
                        self.last_cell_vec = c; return self.last_cell_vec
                self.last_cell_vec = None
                return self.last_cell_vec
            else:
                if debug > 0: print("GET_LAST_CELL_VECTORS: No all_volumes in list")
                self.last_cell_vec = None
                return self.last_cell_vec
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0:
                if debug > 0: print("GET_LAST_CELL_VECTORS: No opt_blocks in list")
                self.last_cell_vec = None
                return self.last_cell_vec
            for idx in range(len(self.opt_blocks)-1,-1,-1):
                init_line = self.opt_blocks[idx][0]+1
                last_line = self.opt_blocks[idx][1]+1
                tmp = parse_cell_vectors(self.lines[init_line:last_line])
                if tmp is not None:
                    self.last_cell_vec = tmp; return self.last_cell_vec
            if debug > 0: print("GET_LAST_CELL_VECTORS: all cell vectors are None")
            self.last_cell_vec = None
            return self.last_cell_vec

##############
### VOLUME ###
##############
    def get_all_volumes(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_volumes = []
        if len(self.opt_blocks) > 0: 
            for idx, b in enumerate(self.opt_blocks):
                step_volume  = parse_volume_from_step(self.lines[b[0]+1:b[1]+1])
                self.all_volumes.append(step_volume)
                if step_volume is None and debug > 0: print("GET_ALL_VOLUMES: error parsing volume between lines", b) 
        else:
            if debug > 0: print("GET_ALL_VOLUMES: empty list of opt_blocks")
            self.all_volumes = None
        return self.all_volumes

    def get_last_volume(self, debug: int=0):
        if hasattr(self,"all_volumes"): 
            if len(self.all_volumes) > 0: 
                for v in self.all_volumes[-1::-1]:
                    if v is not None:
                        self.last_volume = v
                        return self.last_volume
                self.last_volume = None
                return self.last_volume
            else: 
                if debug > 0: print("GET_LAST_VOLUME: No all_volumes in list")
                self.last_volume = None
                return self.last_volume
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0: 
                if debug > 0: print("GET_LAST_VOLUME: No opt_blocks in list")
                self.last_volume = None
                return self.last_volume
            for idx in range(len(self.opt_blocks)-1,-1,-1):
                init_line = self.opt_blocks[idx][0]+1
                last_line = self.opt_blocks[idx][1]+1
                tmp = parse_volume_from_step(self.lines[init_line:last_line])
                if tmp is not None:
                    self.last_volume = tmp; return self.last_volume
            if debug > 0: print("GET_LAST_VOLUME: all volumes are None")
            self.last_volume = None
            return self.last_volume

###############
### DENSITY ###
###############
    def get_all_densities(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_densities = []
        if len(self.opt_blocks) > 0: 
            for idx, b in enumerate(self.opt_blocks):
                step_density  = parse_density_from_step(self.lines[b[0]+1:b[1]+1])
                self.all_densities.append(step_density)
                if step_density is None and debug > 0: print("GET_ALL_VOLUMES: error parsing density between lines", b) 
        else:
            if debug > 0: print("GET_ALL_DENSITIES: empty list of opt_blocks")
            self.all_densities = None
        return self.all_densities

    def get_last_density(self, debug: int=0):
        if hasattr(self,"all_densities"):
            if len(self.all_densities) > 0:
                for d in self.all_densities[-1::-1]:
                    if d is not None:
                        self.last_density = d; return self.last_density
                self.last_density = None 
                return self.last_density
            else:
                if debug > 0: print("GET_LAST_DENSITY: No all_densities in list")
                self.last_density = None 
                return self.last_density
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0:
                if debug > 0: print("GET_LAST_DENSITY: No opt_blocks in list")
                self.last_density = None 
                return self.last_density
            for idx in range(len(self.opt_blocks)-1,-1,-1):
                init_line = self.opt_blocks[idx][0]+1
                last_line = self.opt_blocks[idx][1]+1
                tmp = parse_density_from_step(self.lines[init_line:last_line])
                if tmp is not None:
                    self.last_density = tmp; return self.last_density
            if debug > 0: print("GET_LAST_DENSITY: all densities are None")
            self.last_density = None 
            return self.last_density

##############
### FORCES ###
##############
    def get_all_forces(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_at_forces = []
        self.all_tot_forces = []
        if len(self.opt_blocks) > 0: 
            for idx, b in enumerate(self.opt_blocks):
                step_at_forces, step_tot_force = parse_forces_from_step(self.lines[b[0]+1:b[1]+1])
                self.all_at_forces.append(step_at_forces)
                self.all_tot_forces.append(step_tot_force)
                if (step_at_forces is None or step_tot_force is None) and debug > 0: print("GET_ALL_FORCES: error parsing forces between lines", b) 
        else:
            if debug > 0: print("GET_ALL_FORCES: empty list of opt_blocks")
            self.all_at_forces = None
            self.all_tot_forces = None
        return self.all_at_forces, self.all_tot_forces

    def get_last_forces(self, debug: int=0):
        if hasattr(self,"all_at_forces"): 
            if len(self.all_at_forces) > 0 and len(self.all_tot_forces) > 0:  
                for idx in range(len(self.all_at_forces)-1,-1,-1): 
                    if self.all_at_forces[idx] is not None and self.all_tot_forces is not None: 
                        self.last_at_forces  = self.all_at_forces[idx]
                        self.last_tot_forces = self.all_tot_forces[idx]
                self.last_at_forces = None 
                self.last_tot_forces = None 
                return self.last_at_forces, self.last_tot_forces 
            else: 
                if debug > 0: print("GET_LAST_FORCES: No all_forces in list")
                self.last_at_forces = None 
                self.last_tot_forces = None 
                return self.last_at_forces, self.last_tot_forces 
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0: 
                if debug > 0: print("GET_LAST_FORCES: No opt_blocks in list")
                self.last_at_forces = None 
                self.last_tot_forces = None 
                return self.last_at_forces, self.last_tot_forces 
            for idx in range(len(self.opt_blocks)-1,-1,-1):
                init_line = self.opt_blocks[idx][0]+1
                last_line = self.opt_blocks[idx][1]+1
                tmp1, tmp2 = parse_forces_from_step(self.lines[init_line:last_line])
                if tmp1 is not None and tmp2 is not None:
                    self.last_at_forces = tmp1 
                    self.last_tot_force = tmp2 
                    return self.last_at_forces, self.last_tot_force
            if debug > 0: print("GET_LAST_FORCES: all forces are None")
            self.last_at_forces = None 
            self.last_tot_forces = None 
            return self.last_at_forces, self.last_tot_forces 
            
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
        if hasattr(self,"last_energy"):  to_print += f'Last Energy       = {self.last_energy}\n'
        if hasattr(self,"last_volume"):  to_print += f'Last Volume       = {self.last_volume}\n'
        if hasattr(self,"last_density"): to_print += f'Last Density      = {self.last_density}\n'
        if hasattr(self,"elapsed_time"): to_print += f'Elapsed Time      = {self.elapsed_time}\n'
        to_print += f'---------------------------------------------------\n'
        return to_print

################
#### GENERAL ###
################
#    def execute_all_parsing(self, prop_name: str, function: str, debug: int=0):
#        if function not in dir(self): print("Function not in class"); return None
#        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
#        method = getattr(self,"function")
#        prop = []
#        for idx, b in enumerate(self.opt_blocks):
#            init_line = b[0]+1
#            last_line = b[1]+1
#            step_prop  = method(self.lines[init_line:last_lines])
#            if step_prop is not None: prop.append(step_prop)
#        setattr(self,prop_name,prop)
#        return getattr(self,prop_name)
#
#    def execute_last_parsing(self, prop_name: str, all_prop_name: str, function: str, debug: int=0):
#        if not function in dir(self): return None
#        method = getattr(self,"function")
#        #if not hasattr(self,all_prop_name): self.execute_all_parsing(all_prop_name, function)
#        if hasattr(self,all_prop_name): 
#            var = getattr(self,all_prop_name)
#            if len(var) > 0: 
#                for v in var[-1::-1]:
#                    if v is not None: prop = v; setattr(self,prop_name,prop); return getattr(self,prop_name)
#        else:
#            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
#            if len(self.opt_blocks) == 0: return None
#
#            for idx in range(len(self.opt_blocks),-1,-1):
#                init_line = self.opt_blocks[idx][0]+1
#                last_line = self.opt_blocks[idx][1]+1
#                prop = method(self.lines[init_line:last_line])
#                if prop is not None:
#                    setattr(self,prop_name,prop)
#                    return getattr(self,prop_name)

