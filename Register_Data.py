import sys
import os
import numpy as np
from copy import deepcopy

from Scope.Adapted_from_cell2mol import labels2formula, get_adjmatrix, get_radii, get_blocks, get_molecules, inv
from Scope.Classes_Data import *
from Scope.Classes_State import *
from Scope.Gmol_ops import gmol_update_geom, cell_update_geom, gmol_create_geom, cell_create_geom
from Scope.Parse_General import read_lines_file, search_string 

#from Scope.Software.Gaussian.G16_Class_Output import *
#from Scope.Software.Quantum_Espresso.QE_Class_Output import * 

from Scope import Constants

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
    comp.status                  = comp.output.get_last_block_status()
    if debug > 0: 
        print("REG_GENERAL: comp.isfinished:", comp.isfinished)
        print("REG_GENERAL: comp.elapsed_time:", comp.elapsed_time)
        print("REG_GENERAL: comp.status:", comp.status)

###########################################
def reg_optimization(comp: object, debug: int=0):

    ### For simplicity...
    gmol = comp._job._recipe.subject

    ### 0-In Case Reg_General hasn't been run:
    if not hasattr(comp,"output"): reg_general(comp)

    ### 1-Retrieve the last geometry, and evaluate if it is the converged one (or if it needs to run more)
    labels, new_coord = comp.output.get_geometry_last_complete_block()
    comp.isgood = comp.output.get_optimization_finished()
    if labels is not None:
        assert len(labels) == len(new_coord)
        assert len(labels) > 0

    ### 2-If it is a periodic objecet, gets the cell vectors as well
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
            fstate.add_computation(comp)
            worked = True
        except Exception as exc:
            print("Exception storing results for registry optimization", exc)
    return worked

###########################################
def reg_frequencies(comp: object, debug: int=0):

    ### For simplicity...
    gmol = comp._job._recipe.subject

    ### 0-In Case Reg_General hasn't been run:
    if not hasattr(comp,"output"): reg_general(comp)
    
    worked = False
    ### 1-Parsing ###
    VNMs = comp.output.get_vnms()

    ### 2-Storage ###
    if VNMs is not None: 
        exists, fstate = find_state(gmol, comp.qc_data.fstate)   ## If exists, it will be updated 
        if not exists: fstate = state(gmol, comp.qc_data.fstate)
        fstate.set_VNMs(VNMs)
        fstate.add_computation(comp)
        worked = True
        comp.isgood = True 

    return worked

###########################################
def reg_energy(comp: object, debug: int=0):

    ### For simplicity...
    gmol = comp._job._recipe.subject

    ### 0-In Case Reg_General hasn't been run:
    if not hasattr(comp,"output"): reg_general(comp)

    ### 1-Parses Energy
    energy = comp.output.get_energy_last_complete_block()       ## last_complete_block requires convergence, not necessary energy. Careful
    if debug > 0: print(f"REG_ENERGY: energy is {energy} a.u.")

    ### Storage ###
    #worked = False
    if energy is not None:
        try:
            exists, fstate = find_state(gmol, comp.qc_data.fstate)   ## If exists, it will be updated
            if not exists: fstate = state(gmol, comp.qc_data.fstate)
            fstate.set_energy(energy, 'au')
            fstate.add_computation(comp)
            worked = True
        except Exception as exc:
            print(exc)
    #        worked = False

    #return worked
