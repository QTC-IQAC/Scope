import Scope_New.Constants 
from   Scope_New.Parse_General import search_string, read_lines_file
from . import Parse_G16_outputs 

class g16_output(object):
    def __init__(self, lines: list, computation: object=None):
        self._computation     = computation
        self.lines            = lines
        if computation is not None:
            if hasattr(computation,"jobtype"):           self.jobtype        = computation.jobtype
            elif hasattr(computation.qc_data,"jobtype"): self.jobtype        = computation.qc_data.jobtype
            self.requisites       = self.get_requisites()
        else: 
            self.jobtype        = "unknown"
            self.requisites     = []
        
    def clear_lines(self):
        self.lines          = []
        
    def read_lines(self):
        if hasattr(self,"_computation"): self.lines = read_lines_file(self._computation.out_path)

##################
### REQUISITES ###
##################
    # Not very useful now, the idea is to define the blocks of data that are needed for every jobtype
    def get_requisites(self):
        if   self.jobtype == 'scf':   self.requisites = ['scf'] 
        elif self.jobtype == 'opt':   self.requisites = ['scf','opt'] 
        elif self.jobtype == 'freq':  self.requisites = ['scf','freq'] 
        return self.requisites

###############
### STATUS ###
###############
    def get_status_finished(self):
        self.status_finished = parse_status_finished(self.lines)
        return self.status_finished

    def get_optimization_finished(self, debug: int=0):
        self.optimization_finished = parse_opt_status(self.lines)
        return self.optimization_finished

    def get_scf_finished(self, debug: int=0):
        self.scf_finished = parse_scf_status(self.lines)
        return self.scf_finished

    def get_last_block_status(self, debug: int=0):
        lines_last_opt_block = self.get_lines_last_opt_block()
        if lines_last_opt_block is None: return "aborted"   # When file is empty or killed early
        scf_convergence      = parse_scf_status(lines_last_opt_block)
        coordinates          = parse_coord_status(lines_last_opt_block)
        frequencies          = parse_freq_status(lines_last_opt_block)
        time_limit           = parse_timelimit_status(lines_last_opt_block)

        if debug > 0: print(f"GET_LAST_BLOCK_STATUS: found scf={scf_convergence}, coord={coordinates}, frequencies={frequencies}, timelimit={time_limit}" )

        if   time_limit:                                  self.last_block_status = "time_stopped"
        if   scf_convergence is None:                     self.last_block_status = "aborted"
        elif not scf_convergence:                         self.last_block_status = "no_scf_convergence"
        elif scf_convergence:
            if not coordinates and not frequencies:       self.last_block_status = "no_coord"
            else:                                         self.last_block_status = "worked"
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
                init_line = self.opt_blocks[idx][0]#+1
                last_line = self.opt_blocks[idx][1]#+1
                if debug > 0: print(f"GET_LAST_COMPLETE_BLOCK: evaluating block {idx}. Between {init_line} and {last_line}")
                block_lines = self.lines[init_line:last_line]
                scf_convergence      = parse_scf_status(block_lines)
                coordinates          = parse_coord_status(block_lines)
                frequencies          = parse_freq_status(block_lines)
                time_limit           = parse_timelimit_status(block_lines)
                if debug > 0: print(f"GET_LAST_COMPLETE_BLOCK: found scf={scf_convergence}, coord={coordinates}, frequencies={frequencies}, timelimit={time_limit}" )

                # Evaluates conditions depending on the requisites of this type of job
                good = True
                if good and time_limit:                                                                     good = False
                if good and 'scf'  in self.requisites and (not scf_convergence or scf_convergence is None): good = False  
                if good and 'opt'  in self.requisites and not coordinates:                                  good = False
                if good and 'freq' in self.requisites and not frequencies:                                  good = False
                if good: found = True; self.last_complete_block = block_lines 

#                if 'scf' in self.jobtype:
#                    if not time_limit and scf_convergence:
#                        self.last_complete_block = block_lines
#                        found = True
#                elif 'opt' in self.jobtype:
#                    if not time_limit and scf_convergence and coordinates and frequencies:
#                        self.last_complete_block = block_lines
#                        found = True
#                elif 'freq' in self.jobtype:
#                    if not time_limit and scf_convergence and frequencies:
#                        self.last_complete_block = block_lines
#                        found = True
#                else: print(f"GET_LAST_COMPLETE_BLOCK: unexpected output.jobtype={self.jobtype}")
 
        if not found:
            print("GET_LAST_COMPLETE_BLOCK: No complete block was found")
            self.last_complete_block = None
        return self.last_complete_block

#    def get_last_complete_block(self, debug: int=0):
#        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
#        found = False
#        for idx in range(len(self.opt_blocks)-1,-1,-1):
#            if not found:
#                init_line = self.opt_blocks[idx][0]+1
#                last_line = self.opt_blocks[idx][1]+1
#                block_lines = self.lines[init_line:last_line]
#                scf_convergence      = parse_scf_status(block_lines)
#                coordinates          = parse_coord_status(block_lines)
#                time_limit           = parse_timelimit_status(block_lines)
#                if not time_limit and scf_convergence and coordinates:
#                    self.last_complete_block = block_lines
#                    found = True
#        if not found:
#            if debug > 0: print("GET_LAST_COMPLETE_BLOCK: No complete block was found")
#            self.last_complete_block = None
#        return self.last_complete_block

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
        if not hasattr(self,"scf_blocks"): self.get_scf_blocks()
        self.opt_blocks = []
        if len(self.scf_blocks) > 0:
            for idx, b in enumerate(self.scf_blocks[:-1]):
                if idx == 0: self.opt_blocks.append([0, self.scf_blocks[idx+1][0]-1])
                else:        self.opt_blocks.append([self.scf_blocks[idx][0], self.scf_blocks[idx+1][0]-1])
            self.opt_blocks.append([self.scf_blocks[-1][0],len(self.lines)-1])
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
                step_labels, step_coords       = parse_geometry_from_step(self.lines[b[0]:b[1]])
                #step_labels, step_coords       = parse_geometry_from_step(self.lines[b[0]+1:b[1]+1])
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
                init_line = self.opt_blocks[idx][0]#+1
                last_line = self.opt_blocks[idx][1]#+1
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

##############
### FORCES ###
##############
    def get_all_forces(self, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        self.all_at_forces = []
        self.all_tot_forces = []
        if len(self.opt_blocks) > 0:
            for idx, b in enumerate(self.opt_blocks):
                step_at_forces = parse_forces_from_step(self.lines[b[0]:b[1]])
                #step_at_forces = parse_forces_from_step(self.lines[b[0]+1:b[1]+1])
                self.all_at_forces.append(step_at_forces)
                if step_at_forces is None and debug > 0: print("GET_ALL_FORCES: error parsing forces between lines", b)
        else:
            if debug > 0: print("GET_ALL_FORCES: empty list of opt_blocks")
            self.all_at_forces = None
        return self.all_at_forces

    def get_last_forces(self, debug: int=0):
        if hasattr(self,"all_at_forces"):
            if len(self.all_at_forces) > 0 and len(self.all_tot_forces) > 0:
                for idx in range(len(self.all_at_forces)-1,-1,-1):
                    if self.all_at_forces[idx] is not None and self.all_tot_forces is not None:
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
                init_line = self.opt_blocks[idx][0]#+1
                last_line = self.opt_blocks[idx][1]#+1
                tmp1 = parse_forces_from_step(self.lines[init_line:last_line])
                if tmp1 is not None:
                    self.last_at_forces = tmp1
                    return self.last_at_forces
            if debug > 0: print("GET_LAST_FORCES: all forces are None")
            self.last_at_forces = None
            return self.last_at_forces

    def get_forces_last_complete_block(self, debug: int=0):
        if not hasattr(self,"last_complete_block"): self.get_last_complete_block()
        if self.last_complete_block is None: return None, None
        tmp = parse_forces_from_step(self.last_complete_block)
        if tmp is not None:
            self.last_at_forces = tmp
        else:
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
                step_energy  = parse_energy_from_step(self.lines[b[0]:b[1]])
                #step_energy  = parse_energy_from_step(self.lines[b[0]+1:b[1]+1])
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
                init_line = self.opt_blocks[idx][0]#+1
                last_line = self.opt_blocks[idx][1]#+1
                tmp = parse_energy_from_step(self.lines[init_line:last_line])
                if tmp is not None:
                    self.last_energy = tmp; return self.last_energy
            if debug > 0: print("GET_LAST_ENERGY: all energies are None")
            self.last_energy = None
            return self.last_energy

    def get_energy_last_complete_block(self, debug: int=0):
        if not hasattr(self,"last_complete_block"): self.get_last_complete_block()
        if self.last_complete_block is None: return None
        tmp = parse_energy_from_step(self.last_complete_block)
        if tmp is not None: self.last_energy = tmp
        else:
            if debug > 0: print("get_energy_last_complete_block: energy is None")
            self.last_energy = None
        return self.last_energy

#########################
### FREQUENCIES & VNM ###
#########################
    #def get_frequencies(self, debug: int=0):
    #    if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
    #    if len(self.opt_blocks) == 0: return None
    #    for idx in range(len(self.opt_blocks)-1,-1,-1):
    #        init_line = self.opt_blocks[idx][0]#+1
    #        last_line = self.opt_blocks[idx][1]#+1
    #        tmp = parse_freqs_from_step(self.lines[init_line:last_line])
    #        if tmp is not None:
    #            self.frequencies = tmp; return self.frequencies
    #    self.frequencies = None
    #    return self.frequencies

    def get_vnms(self, witheigen: bool=False, debug: int=0):
        if not hasattr(self,"opt_blocks"): self.get_opt_blocks()
        if len(self.opt_blocks) == 0: return None
        for idx in range(len(self.opt_blocks)-1,-1,-1):
            init_line = self.opt_blocks[idx][0]#+1
            last_line = self.opt_blocks[idx][1]#+1
            ## After June 2025, I force frequencies in Gaussian to be computed with freq=hpmodes. 
            ## As a result, two parsing routines are needed.
            tmp = parse_hp_vnms_from_step(self.lines[init_line:last_line], witheigen=witheigen, debug=debug)
            if tmp is None: tmp = parse_vnms_from_step(self.lines[init_line:last_line], witheigen=witheigen, debug=debug)
            if tmp is not None:    
                self.vnms = tmp; return self.vnms
            else:                  
                if debug > 0: print("get_vnms: vnms is None")
        self.vnms = None
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

    def get_cpu_time(self, debug: int=0):
        line_elapsed, found_elapsed   = search_string("Job cpu time:", self.lines, typ='last')
        if found_elapsed:
            elapsed_time_list = self.lines[line_elapsed].split()[3:]
            self.cpu_time = G16_time_to_sec(elapsed_time_list)
        else: self.cpu_time = float(0)
        return self.cpu_time        
    
#############
### PRINT ###
#############
    def __repr__(self):
        to_print  = f'---------------------------------------------------\n'
        to_print += f'   FORMATTED REPRESENTATION OF OUTPUT CLASS        \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f'#Lines            = {len(self.lines)}\n'
        if hasattr(self,"last_block_status"):  to_print += f'Last Block Status = {self.last_block_status}\n'
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
