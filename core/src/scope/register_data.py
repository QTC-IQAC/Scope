import sys
import os
import numpy as np
from copy import deepcopy

from scope import constants
from scope.connectivity import *
from scope.classes_data import *
from scope.classes_state import *

######################################################################
# 0) HERE WE GATHER THE RULES TO REGISTER THE DIFFERENT TYPES OF JOBS 
# 1) in computations.register() we manage the connection between registration and the other parts of the code
# 2) in Execute_Job, we manage the relationship between the registration, and the computational workflow (Branch, Workflow, Job, Computation)
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
    """Store optimization geometry, available forces and optimization success.

    Keep each Computation's forces before discarding its output object. Store
    forces on the final State when available, in the parser's force units.
    """
    ### 0-In Case Reg_General hasn't been run:
    if not hasattr(comp,"output"): reg_general(comp)

    ### 1-Retrieve the last geometry, and evaluate if it is the converged one (or if it needs to run more)
    labels, new_coord = comp.output.get_geometry_last_complete_block()
    if new_coord is None:
        print(f"REG_OPTIMIZATION: received None coordinates from comp.output. State won't be updated")
    comp.isgood = comp.output.get_optimization_finished()
    comp.forces = comp.output.get_forces_last_complete_block()
    if labels is not None:
        assert len(labels) == len(new_coord)
        assert len(labels) > 0

    ### 2-If it is a periodic object, gets the cell vectors as well
    if comp.source.object_type == "cell": 
        cellvec, celldim, cellparam = comp.output.get_cell_vectors_last_complete_block()
        if cellvec is None: cellvec, celldim, cellparam = comp.output.get_last_cell_vectors()
        if cellvec is None: print("Couldn't Extract cell vectors from output", comp.out_path)

    worked = False
    if new_coord is not None:
        ### 3a-Stores Results in the State Object
        exists, fstate        = comp.source.find_state(comp.qc_data.fstate, debug=debug)  ## If exists, it will be updated
        if not exists: fstate = comp.source.add_state(comp.qc_data.fstate, debug=debug)   ## Otherwise, it is created
        print(f"REG_OPTIMIZATION: Setting new geometry on fstate={fstate.name}")
        fstate.set_geometry(labels, new_coord)                                            ## New geometry is stored
        if comp.forces is not None: fstate.set_forces(comp.forces)                        ## New forces are stored, if available
        if comp.source.object_type == "cell": 
            print(f"REG_OPTIMIZATION: Setting new cell parameters on fstate={fstate.name}")
            fstate.set_cell(cellvec, cellparam)
            print(f"REG_OPTIMIZATION: Getting molecules for fstate={fstate.name} and checking fragmentation")
            fstate.get_molecules()
            fstate.check_fragmentation(reconstruct=True, debug=debug)
        ### 3b-Computation is linked to the fstate
        fstate.add_computation(comp)
        worked = True
    return worked

###########################################
def reg_frequencies(comp: object, witheigen: bool=False, debug: int=0):

    ### 0-In Case Reg_General hasn't been run:
    if not hasattr(comp,"output"): reg_general(comp)
    
    worked = False
    ### 1-Parsing VNMs and Forces ###
    VNMs   = comp.output.get_vnms(witheigen=witheigen)
    forces = comp.output.get_forces_last_complete_block()

    ### 2a-Saving Data in fstate ###
    exists, fstate        = comp.source.find_state(comp.qc_data.fstate, debug=debug)  ## If exists, it will be updated
    if not exists: fstate = comp.source.add_state(comp.qc_data.fstate, debug=debug)   ## Otherwise, it is created 
    ### 2b-Computation is linked to the fstate
    fstate.add_computation(comp)

    if VNMs is not None:
        fstate.set_VNMs(VNMs)
        worked = True
        comp.isgood = True 
        if debug > 0: print("REG_FREQUENCIES: VNMs were parsed")
    else: print("REG_FREQUENCIES: could not parse frequencies")

    ### Forces do not condition "worked" nor "comp.isgood"
    if forces is not None: 
        fstate.set_forces(forces)
        if debug > 0: print("REG_FREQUENCIES: atomic forces were parsed")
    else: print("REG_FREQUENCIES: could not parse forces")

    return worked

###########################################
def reg_energy(comp: object, debug: int=0):

    ### 0-In Case Reg_General hasn't been run:
    if not hasattr(comp,"output"): reg_general(comp)

    ### 1-Parses Energy
    energy      = comp.output.get_energy_last_complete_block()
    comp.energy = energy  # Hartree; keep this attempt independent of later State updates.
    comp.isgood = comp.output.get_scf_finished()
    if debug > 0: print(f"REG_ENERGY: {comp.isgood=}")          ## Stores isgood, which determines if the computation converged, essentially
    if debug > 0: print(f"REG_ENERGY: {energy=} a.u.")          ## Parsing routines already convert energy to a.u.

    ## 2-Try to parse Forces if they're available, typically for finite differences
    if hasattr(comp.qc_data,"print_forces"):
        if comp.qc_data.print_forces: forces = comp.output.get_forces_last_complete_block()
        else:                         forces = None
    else:                             forces = None

    if getattr(comp.qc_data,'print_forces',False): 
        comp.forces = forces  # Keep this attempt's forces independent of State updates.

    ### Storage ###
    worked = False
    if energy is not None or forces is not None:
        # Saves data
        exists, fstate        = comp.source.find_state(comp.qc_data.fstate, debug=debug)  ## If exists, it will be updated
        if not exists: fstate = comp.source.add_state(comp.qc_data.fstate, debug=debug)   ## Otherwise, it is created 
        if energy is not None: 
            print(f"REG_ENERGY: Setting energy on fstate={fstate.name}")
            fstate.set_energy(energy, 'au')
        if forces is not None: 
            print(f"REG_ENERGY: Setting forces on fstate={fstate.name}")
            fstate.set_forces(forces)
        # Computation is linked to the fstate
        fstate.add_computation(comp)
        worked = True
    return worked

###########################################
def reg_excited_states(comp: object, debug: int=0):

    ### 0-In Case Reg_General hasn't been run:
    if not hasattr(comp,"output"): reg_general(comp)
    
    worked = False
    ### 1-Parsing VNMs and Forces ###
    exc_states = comp.output.get_exc_states(nstates=comp.qc_data.td_nstates, debug=debug)

    ### 2a-Saving Data in fstate ###
    exists, fstate        = comp.source.find_state(comp.qc_data.fstate, debug=debug)  ## If exists, it will be updated
    if not exists: fstate = comp.source.add_state(comp.qc_data.fstate, debug=debug)   ## Otherwise, it is created 

    ### 2b-Computation is linked to the fstate
    fstate.add_computation(comp)

    if exc_states is not None:
        fstate.set_exc_states(exc_states)
        worked = True
        comp.isgood = True 
        if debug > 0: print("REG_EXCITED_STATES: excited states were parsed")
    else: print("REG_EXCITED_STATES: could not parse excited states")

    return worked


def reg_findiff(job: object, debug: int=0) -> bool:
    """Register a complete QE finite-difference set of computations on job.fstate.

    Requires an undisplaced 'ref' SCF and all displacement pairs. The final State
    receives the reference geometry/cell, its SCF energy and forces, 3N-3 VNMs
    (excluding translations), and links to the contributing Computations.
    The function Returns False for incomplete or inconsistent sets of computations. 
    Negative frequencies do not invalidate registration.
    """
    from scope.findiff import get_VNM_from_findiff
    found, reference = job.source.find_state(job.istate)
    if not found: return False
    current = sorted((comp for comp in job.computations if not comp.has_update), key=lambda comp: comp.keyword)
    ref_comp = next((comp for comp in current if comp.keyword == 'ref'),None)
    if ref_comp is None: return False
    if not ref_comp.isregistered or not getattr(ref_comp,'isgood',False) or not getattr(ref_comp,'isfinished',False): return False
    if getattr(ref_comp,'energy',None) is None or not np.isfinite(ref_comp.energy) or getattr(ref_comp,'forces',None) is None: return False
    try:
        attempts = tuple((comp.keyword, comp.step, comp.run_number, getattr(comp,'registration_time',None)) for comp in current)
        geometry  = (tuple(reference.labels), np.asarray(reference.coord).tobytes(), np.asarray(reference.cell_vector).tobytes())
        settings  = (job.fstate, getattr(job.job_data,'findiff_displacement',0.01), getattr(job.job_data,'findiff_units','bohr'))
        signature = (attempts, geometry, settings)
        exists, fstate = job.source.find_state(job.fstate)
        if signature == getattr(job,'_findiff_signature',None) and exists and hasattr(fstate,'VNMs') and len(fstate.VNMs) == 3*len(reference.labels)-3: return True
        VNMs = get_VNM_from_findiff(job, debug=debug)
    except (OSError, ValueError) as exc:
        print(f"REG_FINDIFF: {exc}")
        return False
    if VNMs is None: return False

    if not exists: fstate = job.source.add_state(job.fstate)
    fstate.set_geometry(list(reference.labels), np.array(reference.coord, dtype=float, copy=True))
    fstate.set_cell(deepcopy(reference.cell_vector), deepcopy(reference.cell_param))
    fstate.set_energy(ref_comp.energy,'au')
    fstate.set_forces(np.array(ref_comp.forces,dtype=float,copy=True))
    fstate.set_VNMs(VNMs)
    for comp in current: fstate.add_computation(comp)
    job._findiff_signature = signature
    job._workflow._branch._system.should_save = True
    if debug > 0: print(f"REG_FINDIFF: {len(VNMs)} modes registered on State {fstate.name}")
    return True
