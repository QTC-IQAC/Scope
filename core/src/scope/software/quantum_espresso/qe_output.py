from scope.parse_general     import read_lines_file
from scope.software.quantum_espresso.qe_parse  import *

class QE_output(object):
    def __init__(self, lines: list, computation: object=None):
        self._computation   = computation
        self.lines          = lines
        self.set_jobtype()
        self.get_requisites()

    def clear_lines(self):
        self.lines          = []
        
    def read_lines(self):
        if hasattr(self,"_computation"): self.lines = read_lines_file(self._computation.out_path)
        
    def set_jobtype(self, jobtype: str=None):
        if jobtype is not None:
            self.jobtype = jobtype
        else:
            if self._computation is not None:
                if hasattr(self._computation,"jobtype"):             self.jobtype        = self._computation.jobtype
                elif hasattr(self._computation,"qc_data"):
                    if hasattr(self._computation.qc_data,"comp_type"): self.jobtype      = self._computation.qc_data.comp_type
                else:                                          self.jobtype        = "unknown"
            else:                                              self.jobtype        = "unknown"
        return self.jobtype

##################
### REQUISITES ###
##################
    # Not very useful now, the idea is to define the blocks of data that are needed for every jobtype
    def get_requisites(self):
        if   self.jobtype == 'scf':       self.requisites = ['scf'] 
        elif self.jobtype == 'relax':     self.requisites = ['scf','opt'] 
        elif self.jobtype == 'vc-relax':  self.requisites = ['scf','opt'] 
        else:                             self.requisites = []
        return self.requisites

###############
### STATUS ###
###############
    def get_status_finished(self):
        self.status_finished = parse_status_finished(self.lines)
        return self.status_finished

    def get_optimization_finished(self, debug: int=0):
        status = parse_opt_status(self.lines)
        if   status == 'finished':      self.optimization_finished = True
        elif status == 'stucked' :      self.optimization_finished = True
        elif status == 'not finished' : self.optimization_finished = False
        return self.optimization_finished

    def get_scf_finished(self, debug: int=0):
        self.scf_finished = parse_scf_status(self.lines)
        return self.scf_finished

    def get_last_block_status(self, debug: int=0):
        lines_last_opt_block = self.get_lines_last_opt_block()
        if lines_last_opt_block is None: return "aborted"   # When file is empty or killed early
        scf_convergence      = parse_scf_status(lines_last_opt_block)
        coordinates          = parse_coord_status(lines_last_opt_block)
        time_limit           = parse_timelimit_status(lines_last_opt_block)

        if   time_limit:              self.last_block_status = "time_stopped" 
        if   scf_convergence is None: self.last_block_status = "aborted"
        elif not scf_convergence:     self.last_block_status = "no_scf_convergence"
        elif scf_convergence: 
            if not coordinates:       self.last_block_status = "no_coord"
            else:                     self.last_block_status = "worked"
        return self.last_block_status

###############
### BLOCKS ###
###############
    def get_last_complete_block(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks(debug=debug)
        found = False
        if debug > 0: print(f"GET_LAST_COMPLETE_BLOCK: evaluating with jobtype={self.jobtype}")
        for idx in range(len(self.opt_blocks)-1,-1,-1):
            if not found:
                init_line = self.opt_blocks[idx][0]+1
                last_line = self.opt_blocks[idx][1]+1
                if debug > 0: print(f"GET_LAST_COMPLETE_BLOCK: evaluating block {idx}. Between {init_line} and {last_line}")
                block_lines = self.lines[init_line:last_line]
                scf_convergence      = parse_scf_status(block_lines)
                coordinates          = parse_coord_status(block_lines)
                time_limit           = parse_timelimit_status(block_lines)
                if debug > 0: print(f"GET_LAST_COMPLETE_BLOCK: found scf={scf_convergence}, coord={coordinates}, timelimit={time_limit}" )

                # Evaluates conditions depending on the requisites of this type of job
                good = True
                if good and time_limit:                                                                     good = False
                if good and 'scf'  in self.requisites and (not scf_convergence or scf_convergence is None): good = False  
                if good and 'opt'  in self.requisites and not coordinates:                                  good = False
                if good: found = True; self.last_complete_block = block_lines 
 
                # OLD WAY
#                if 'scf' in self.jobtype:
#                    if not time_limit and scf_convergence:
#                        self.last_complete_block = block_lines
#                        found = True
#                else:
#                    if not time_limit and scf_convergence and coordinates:
#                        self.last_complete_block = block_lines
#                        found = True
        if not found: 
            if debug > 0: print("GET_LAST_COMPLETE_BLOCK: No complete block was found")
            self.last_complete_block = None
        return self.last_complete_block
            
    def get_scf_blocks(self, debug: int=0):
        start_SCF, found1  = parse_start_scf(self.lines)
        end_SCF, found2    = parse_end_scf(self.lines)
        last_start = 0
        last_end   = 0
        self.scf_blocks = []        
        ## Finds Blocks
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
        if not hasattr(self,"scf_blocks"): self.get_scf_blocks(debug=debug)
        self.opt_blocks = []
        if len(self.scf_blocks) > 0: 
            for idx, b in enumerate(self.scf_blocks[:-1]):
                if idx == 0: self.opt_blocks.append([0, self.scf_blocks[idx+1][0]])
                else:        self.opt_blocks.append([self.scf_blocks[idx][0], self.scf_blocks[idx+1][0]])
            self.opt_blocks.append([self.scf_blocks[-1][0],len(self.lines)]) 
        return self.opt_blocks
    
    def get_lines_last_scf_block(self, debug: int=0):
        if not hasattr(self,"scf_blocks"): self.get_scf_blocks()
        if len(self.scf_blocks) == 0: return None
        init_line = self.scf_blocks[-1][0]
        last_line = self.scf_blocks[-1][1]
        return self.lines[init_line:last_line]

    def get_lines_last_opt_block(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        if len(self.opt_blocks) == 0: return None
        init_line = self.opt_blocks[-1][0]
        last_line = self.opt_blocks[-1][1]
        return self.lines[init_line:last_line]

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

    def get_geometry_last_complete_block(self, debug: int=0):
        if not hasattr(self,"last_complete_block"): self.get_last_complete_block()
        if self.last_complete_block is None: return None, None
        tmp1, tmp2 = parse_geometry_from_step(self.last_complete_block)
        if tmp1 is not None and tmp2 is not None:
            self.last_labels = tmp1
            self.last_coords = tmp2
        else:
            if debug > 0: print("GET_GEOMETRY_LAST_COMPLETE_BLOCK: geometry is None")
            self.last_labels = None
            self.last_coords = None
        return self.last_labels, self.last_coords

###################
### CELL PARAMS ###
###################
    def get_all_cell_vectors(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_cellvec = []
        self.all_celldim = []
        self.all_cellparam = []
        if len(self.opt_blocks) > 0: 
            for idx, b in enumerate(self.opt_blocks):
                #step_cell_vec  = parse_cell_vectors(self.lines[b[0]+1:b[1]+1])
                step_cellvec, step_celldim, step_cellparam  = parse_cell_vectors(self.lines[b[0]+1:b[1]+1])
                if step_cellvec is None:
                    step_cellvec, step_cellparam = parse_cell_parameters(self.lines[b[0]+1:b[1]+1])
                    step_celldim = None     ## I'm not sure I can get celldim from cellparameters
                self.all_cellvec.append(step_cellvec)
                self.all_celldim.append(step_celldim)
                self.all_cellparam.append(step_cellparam)
                if step_cellvec is None and debug > 0: print("GET_ALL_CELL_PARAM: error parsing cell_parameters between lines", b)
        else:
            if debug > 0: print("GET_ALL_CELL_PARAM: empty list of opt_blocks")
            self.all_cellvec = None
            self.all_celldim = None
            self.all_cellparam = None
        return self.all_cellvec, self.all_celldim, self.all_cellparam

    def get_last_cell_vectors(self, debug: int=0):
        if hasattr(self,"all_cellvec"):
            if len(self.all_cellvec) > 0:
                for idx in range(len(self.all_cellvec)-1,-1,-1): 
                    if self.all_cellvec[idx] is not None:
                        self.last_cellvec = self.all_cellvec[idx]
                        self.last_celldim = self.all_celldim[idx]
                        self.last_cellparam = self.all_cellparam[idx]
                        return self.last_cellvec, self.last_celldim, self.last_cellparam
                self.last_cellvec = None
                self.last_celldim = None
                self.last_cellparam = None
                return self.last_cellvec, self.last_celldim, self.last_cellparam
            else:
                if debug > 0: print("GET_LAST_CELL_VECTORS: No cellvectors in list")
                self.last_cellvec = None
                self.last_celldim = None
                self.last_cellparam = None
                return self.last_cellvec, self.last_celldim, self.last_cellparam
        else:
            if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
            if len(self.opt_blocks) == 0:
                if debug > 0: print("GET_LAST_CELL_VECTORS: No opt_blocks in list")
                self.last_cellvec = None
                self.last_celldim = None
                self.last_cellparam = None
                return self.last_cellvec, self.last_celldim, self.last_cellparam
            for idx in range(len(self.opt_blocks)-1,-1,-1):
                init_line = self.opt_blocks[idx][0]+1
                last_line = self.opt_blocks[idx][1]+1
                tmp1, tmp2, tmp3 = parse_cell_vectors(self.lines[init_line:last_line])
                if tmp1 is None:
                    tmp1, tmp3 = parse_cell_parameters(self.lines[init_line:last_line])
                    tmp2 = None # I'm not sure I can get celldim from cellparameters
                if tmp1 is not None:
                    self.last_cellvec = tmp1 
                    self.last_celldim = tmp2 
                    self.last_cellparam = tmp3 
                    return self.last_cellvec, self.last_celldim, self.last_cellparam
            if debug > 0: print("GET_LAST_CELL_VECTORS: all cell vectors are None")
            self.last_cellvec = None
            self.last_celldim = None
            self.last_cellparam = None
            return self.last_cellvec, self.last_celldim, self.last_cellparam

    def get_cell_vectors_last_complete_block(self, debug: int=0):
        if not hasattr(self,"last_complete_block"): self.get_last_complete_block()
        if self.last_complete_block is None: return None, None, None
        tmp1, tmp2, tmp3 = parse_cell_vectors(self.last_complete_block)
        if tmp1 is None:
            tmp1, tmp3 = parse_cell_parameters(self.last_complete_block)
            tmp2 = None # I'm not sure I can get celldim from cellparameters
        if tmp1 is not None:
            self.last_cellvec = tmp1 
            self.last_celldim = tmp2 
            self.last_cellparam = tmp3 
        else:
            if debug > 0: print("get_cell_vectors_last_complete_block: cell vector is None")
            self.last_cellvec = None
            self.last_celldim = None
            self.last_cellparam = None
        return self.last_cellvec, self.last_celldim, self.last_cellparam

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

    def get_volume_last_complete_block(self, debug: int=0):
        if not hasattr(self,"last_complete_block"): self.get_last_complete_block()
        if self.last_complete_block is None: return None
        tmp = parse_volume_from_step(self.last_complete_block)
        if tmp is not None: self.last_volume = tmp
        else:
            if debug > 0: print("get_volume_last_complete_block: volume is None")
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

    def get_density_last_complete_block(self, debug: int=0):
        if not hasattr(self,"last_complete_block"): self.get_last_complete_block()
        if self.last_complete_block is None: return None
        tmp = parse_density_from_step(self.last_complete_block)
        if tmp is not None: self.last_density = tmp
        else:
            if debug > 0: print("get_density_last_complete_block: density is None")
            self.last_density = None
        return self.last_density

##############
### FORCES ###
##############
    def get_all_forces(self, with_total: bool=False, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_at_forces = []
        self.all_tot_forces = []
        if len(self.opt_blocks) > 0: 
            for idx, b in enumerate(self.opt_blocks):
                step_at_forces = parse_forces_from_step(self.lines[b[0]+1:b[1]+1])
                step_tot_force = parse_total_force_from_step(self.lines[b[0]+1:b[1]+1])
                self.all_at_forces.append(step_at_forces)
                self.all_tot_forces.append(step_tot_force)
                if (step_at_forces is None or step_tot_force is None) and debug > 0: print("GET_ALL_FORCES: error parsing forces between lines", b) 
        else:
            if debug > 0: print("GET_ALL_FORCES: empty list of opt_blocks")
            self.all_at_forces = None
            self.all_tot_forces = None
        if with_total: return self.all_at_forces, self.all_tot_forces
        else:          return self.all_at_forces

    def get_last_forces(self, with_total: bool=False, debug: int=0):
        if hasattr(self,"all_at_forces"): 
            if len(self.all_at_forces) > 0 and len(self.all_tot_forces) > 0:  
                for idx in range(len(self.all_at_forces)-1,-1,-1): 
                    if self.all_at_forces[idx] is not None and self.all_tot_forces is not None: 
                        self.last_at_forces  = self.all_at_forces[idx]
                        self.last_tot_forces = self.all_tot_forces[idx]
                self.last_at_forces = None 
                self.last_tot_forces = None 
                if with_total: return self.last_at_forces, self.last_tot_forces 
                else:          return self.last_at_forces
            else: 
                if debug > 0: print("GET_LAST_FORCES: No all_forces in list")
                self.last_at_forces = None 
                self.last_tot_forces = None 
                if with_total: return self.last_at_forces, self.last_tot_forces 
                else:          return self.last_at_forces
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
                tmp1 = parse_forces_from_step(self.lines[init_line:last_line])
                tmp2 = parse_total_force_from_step(self.lines[init_line:last_line])
                if tmp1 is not None and tmp2 is not None:
                    self.last_at_forces = tmp1 
                    self.last_tot_force = tmp2 
                    if with_total: return self.last_at_forces, self.last_tot_forces 
                    else:          return self.last_at_forces
            if debug > 0: print("GET_LAST_FORCES: all forces are None")
            self.last_at_forces = None 
            self.last_tot_forces = None 
            if with_total: return self.last_at_forces, self.last_tot_forces 
            else:          return self.last_at_forces

    def get_forces_last_complete_block(self, with_total: bool=False, debug: int=0):
        if not hasattr(self,"last_complete_block"): self.get_last_complete_block()
        if self.last_complete_block is None: return None, None
        tmp1 = parse_forces_from_step(self.last_complete_block)
        tmp2 = parse_total_force_from_step(self.last_complete_block)
        if tmp1 is not None and tmp2 is not None:
            self.last_at_forces = tmp1 
            self.last_tot_force = tmp2 
        else:
            if debug > 0: print("get_forces_last_complete_block: forces are None")
            self.last_at_forces = None
            self.last_tot_force = None 
        if with_total: return self.last_at_forces, self.last_tot_forces 
        else:          return self.last_at_forces
            
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

    def get_energy_last_complete_block(self, debug: int=0):
        if not hasattr(self,"last_complete_block"): self.get_last_complete_block()
        if self.last_complete_block is None: return None
        ## If Makov Payne Energy has been requested, it is taken
        tmp = parse_makov_from_step(self.last_complete_block)
        if tmp is None: 
            tmp = parse_energy_from_step(self.last_complete_block)
        if tmp is not None: 
            self.last_energy = tmp
        else:
            if debug > 0: print("get_energy_last_complete_block: energy is None")
            self.last_energy = None
        return self.last_energy

############
### TIME ###
############
    def get_elapsed_time(self, debug: int=0):
        self.elapsed_time = parse_elapsed_time(self.lines)
        return self.elapsed_time        

    def get_cpu_time(self, debug: int=0):
        self.cpu_time = parse_cpu_time(self.lines)
        return self.cpu_time        
    
#############
### PRINT ###
#############
    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print += f'          SCOPE QEspresso OUTPUT CLASS             \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' #Lines           = {len(self.lines)}\n'
        to_print += f' Job Type         = {self.jobtype}\n'
        to_print += f' Requisites       = {self.requisites}\n'
        if hasattr(self,"last_block_status"):  to_print += f' Last Block Status = {self.last_block_status}\n'
        if hasattr(self,"scf_blocks"):    to_print += f' #SCF Blocks       = {len(self.scf_blocks)}\n'
        if hasattr(self,"opt_blocks"):    to_print += f' #OPT Blocks       = {len(self.opt_blocks)}\n'
        if hasattr(self,"all_coords"):    to_print += f' #Coordinates      = {len(self.all_coords)}\n'
        if hasattr(self,"all_energies"):  to_print += f' #Energies         = {len(self.all_energies)}\n'
        if hasattr(self,"last_energy"):   to_print += f' Last Energy       = {self.last_energy}\n'
        if hasattr(self,"last_volume"):   to_print += f' Last Volume       = {self.last_volume}\n'
        if hasattr(self,"last_density"):  to_print += f' Last Density      = {self.last_density}\n'
        if hasattr(self,"elapsed_time"):  to_print += f' Elapsed Time      = {self.elapsed_time}\n'
        to_print += f'---------------------------------------------------\n'
        return to_print
