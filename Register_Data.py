import sys
import os
import numpy as np
from copy import deepcopy

from . import Constants
from .Connectivity import *
from .Classes_Data import *
from .Classes_State import *
#from .Gmol_ops import gmol_update_geom, cell_update_geom, gmol_create_geom, cell_create_geom
#from .Parse_General import read_lines_file, search_string 

######################################################################
# 0) HERE WE GATHER THE RULES TO REGISTER THE DIFFERENT TYPES OF JOBS ##
# 1) in computations.register() we manage the connection between registration and the other parts of the code
# 2) in Execute_Job, we manage the relationship between the registration and the flow of the "Job", "Recipe" and "Branch"

# NOTE: we could consider merging the computations.register with this (0+1).
######################################################################

def reg_general(comp: object, debug: int=0):

    if not hasattr(comp,"output"): comp.create_output() 
    comp.isfinished              = comp.output.get_status_finished()
    comp.elapsed_time            = comp.output.get_elapsed_time() 
    comp.cpu_time                = comp.output.get_cpu_time() 
    comp.status                  = comp.output.get_last_block_status()
    if debug > 0: 
        print("REG_GENERAL: comp.isfinished:", comp.isfinished)
        print("REG_GENERAL: comp.cpu_time:", comp.cpu_time)
        print("REG_GENERAL: comp.status:", comp.status)

###########################################
def reg_optimization(comp: object, debug: int=0):

    ### For simplicity...
    gmol = comp._job._recipe.source

    ### 0-In Case Reg_General hasn't been run:
    if not hasattr(comp,"output"): reg_general(comp)

    ### 1-Retrieve the last geometry, and evaluate if it is the converged one (or if it needs to run more)
    labels, new_coord = comp.output.get_geometry_last_complete_block()
    comp.isgood = comp.output.get_optimization_finished()
    if labels is not None:
        assert len(labels) == len(new_coord)
        assert len(labels) > 0

    ### 2-If it is a periodic object, gets the cell vectors as well
    if gmol.type == "cell": 
        cellvec, celldim, cellparam = comp.output.get_cell_vectors_last_complete_block()
        if cellvec is None: cellvec, celldim, cellparam = comp.output.get_last_cell_vectors()
        if cellvec is None: print("Couldn't Extract cell vectors from output", comp.out_path)

    worked = False
    if new_coord is not None:
        try: 
            ### 2-Stores Results in the State Object
            exists, fstate = find_state(gmol, comp.qc_data.fstate)   ## If exists, it will be updated 
            if not exists: fstate = state(gmol, comp.qc_data.fstate)
            fstate.set_geometry(labels, new_coord)
            fstate.set_spin_config(comp.spin_config)
            if gmol.type == "cell": 
                fstate.set_cell(cellvec, cellparam)
                fstate.get_moleclist()
                fstate.check_fragmentation(reconstruct=True, debug=debug)
            #fstate.add_computation(comp)
            worked = True
        except Exception as exc:
            print("Exception storing results for registry optimization", exc)
    return worked

###########################################
def reg_frequencies(comp: object, witheigen: bool=False, debug: int=0):

    ### For simplicity...
    gmol = comp._job._recipe.source

    ### 0-In Case Reg_General hasn't been run:
    if not hasattr(comp,"output"): reg_general(comp)
    
    worked = False
    ### 1-Parsing ###
    VNMs = comp.output.get_vnms(witheigen=witheigen)
    if VNMs is None:
        print("REG_FREQUENCIES: could not parse frequencies")

    ### 2-Storage ###
    if VNMs is not None:
        exists, fstate = find_state(gmol, comp.qc_data.fstate)   ## If exists, it will be updated 
        if not exists: fstate = state(gmol, comp.qc_data.fstate)
        fstate.set_VNMs(VNMs)
        fstate.add_computation(comp)
        worked = True
        comp.isgood = True 

    ### Also forces, but they do not condition "worked" or "isgood"
    forces = comp.output.get_forces_last_complete_block()
    if forces is not None:
        exists, fstate = find_state(gmol, comp.qc_data.fstate)   ## If exists, it will be updated 
        fstate.set_forces(forces)

    return worked

###########################################
def reg_energy(comp: object, debug: int=0):

    ### For simplicity...
    gmol = comp._job._recipe.source

    ### 0-In Case Reg_General hasn't been run:
    if not hasattr(comp,"output"): reg_general(comp)

    ### 1-Parses Energy
    energy = comp.output.get_energy_last_complete_block()       ## last_complete_block requires convergence, not necessary energy. Careful
    comp.isgood = comp.output.get_scf_finished()
    print("REG_ENERGY:", energy, comp.isgood)
    if debug > 0: print(f"REG_ENERGY: energy is {energy} a.u.") ## Parsing routines already convert energy to a.u.

    ## 2-Try to parse Forces if they're available, typically for finite differences
    if hasattr(comp.qc_data,"print_forces"):
        if comp.qc_data.print_forces: forces = comp.output.get_forces_last_complete_block()
        else:                         forces = None
    else:                         forces = None

    ### Storage ###
    worked = False
    if energy is not None or forces is not None:
        exists, fstate = find_state(gmol, comp.qc_data.fstate)   ## If exists, it will be updated
        if not exists: fstate = state(gmol, comp.qc_data.fstate)
        if energy is not None: fstate.set_energy(energy, 'au')
        if forces is not None: fstate.set_forces(forces)
        fstate.add_computation(comp)
        worked = True

    return worked
