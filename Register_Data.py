import sys
import os
import numpy as np
from copy import deepcopy

from Scope.Adapted_from_cell2mol import labels2formula, get_adjmatrix, get_radii, get_blocks, get_molecules, inv
from Scope.Classes_Data import *
from Scope.Classes_State import *
from Scope.Gmol_ops import gmol_update_geom, cell_update_geom, gmol_create_geom, cell_create_geom
from Scope.Parse_General import read_lines_file, search_string 
from Scope.Software.Gaussian.Parse_G16_outputs import *
from Scope.Software.Quantum_Espresso.Parse_QE_outputs import * 
from Scope.Software.Quantum_Espresso.QE_Class_Output import * 
from Scope import Constants

###########################################
def reg_general(comp: object, debug: int=0):

    if not hasattr(comp,'output_lines'): comp.read_lines()
    lines = comp.output_lines
   
    ###################
    ### Gaussian 16 ###
    ###################
    if comp.software == 'g16': 
        line_normal, found_normal = search_string("Normal termination", lines, typ='all')
        line_error, found_error   = search_string("Error termination", lines, typ='all')
        #line_route, found_route   = search_string("Route card not found", lines, typ='all')
        line_time, found_time     = search_string("Elapsed time:", lines, typ='last')
        if found_time:
            elapsed_time_list = lines[line_time].split()[2:]
            comp.elapsed_time = G16_time_to_sec(elapsed_time_list)
        else: comp.elapsed_time = float(0)

        if found_error or found_normal:  comp.isfinished = True
        else:                            comp.isfinished = False

        if found_normal:                 comp.isgood = True
        else:                            comp.isgood = False

    #########################
    ### Quantum Espresso ###
    #########################
    elif comp.software == "qe": 

        ## With Output Class
        new_output = qe_output(lines, comp)

        status = new_output.get_status()
        if status == "worked": 
            comp.isfinished = True
            comp.isgood     = True
        else:
            comp.isfinished = False
            comp.isgood     = False

        comp.elapsed_time = new_output.get_elapsed_time()

    #############
    ### Other ###
    #############
    else:
        comp.isfinished   = False
        comp.isgood       = False
        comp.elapsed_time = float(0)
        print(f"    REG_GENERAL: Registry of {comp.software} comps is not implemented.")

###########################################
def reg_optimization(comp: object, debug: int=0):

    ### For simplicity...
    gmol = comp._job._recipe.subject

    ### Reads output
    comp.read_lines()
    lines = comp.output_lines

    worked = False
    ###############
    ## G09 & G16 ##
    ###############
    if comp.software == "g16": 
        new_coord = G16_get_last_geom(lines, debug=debug)
        labels = gmol.labels.copy()

    ########
    ## QE ##
    ########
    elif comp.software == "qe":

        ## This means that the optimization hasn't finished. 
        ## Even in that case, It will still try to retrieve the last geometry 
        line_BFGS, found_BFGS = search_string("End of BFGS Geometry Optimization", lines, typ='last')
        if not found_BFGS: comp.isgood = False    

        ## Extract Coordinates 
        tmp_labels, new_coord = parse_final_geometry(lines, debug=debug)
        assert len(tmp_labels) == len(new_coord)
        assert len(tmp_labels) > 0

        ## Originally, labels include digits to follow the spin state of metal atoms. For instance 'Fe4' indicates a HS Fe atom
        ## These digits must be removed from labels when storing the data, since the digit is only for QE  
        labels = []
        for l in tmp_labels:
            labels.append(str(''.join([c for c in l if not c.isdigit()])))

        # Extracts Cell Parameters
        if gmol.type == "cell": cellvec, celldim, cellparam = get_cell_vectors(lines, debug=debug) 

    #############
    ### Other ###
    #############
    else: print(f"    REG_OPTIMIZATION: Registry of {comp.software} comps is not implemented.")

    ######################
    ### Stores Results ###
    ######################
    assert len(labels) > 0
    assert len(new_coord) > 0
    assert len(new_coord) == len(labels)

    if len(new_coord) > 0: 

        ## Stores data in the corresponding state-class object
        exists, fstate = find_state(gmol, comp.qc_data.fstate)   ## If exists, it will be updated 
        if not exists: fstate = state(gmol, comp.qc_data.fstate)
        fstate.set_geometry(labels, new_coord)
        fstate.set_spin_config(comp.spin_config)
        if gmol.type == "cell": fstate.set_cell(cellvec, cellparam)
        if gmol.type == "cell": fstate.get_moleclist()
        if gmol.type == "cell": fstate.check_fragmentation(reconstruct=True, debug=debug)

        fstate.add_computation(comp)
        worked = True
    else: print("    REG_OPT: empty labels and positions received. Job could not be registered") 

    return worked

###########################################
def reg_frequencies(comp: object, debug: int=0):
    gmol = comp._job._recipe.subject

    ### Reads output
    comp.read_lines()
    lines = comp.output_lines

    ### Parsing ###
    worked = False
    if comp.software == "g16": 
        VNMs = G16_get_VNM(lines, witheigen=True, debug=debug)
    else: print(f"    REG_FREQUENCIES: Registry of {comp.software} comps is not implemented.")

    ### Storage ###
    exists, fstate = find_state(gmol, comp.qc_data.fstate)   ## If exists, it will be updated 
    if not exists: fstate = state(gmol, comp.qc_data.fstate)
    fstate.set_VNMs(VNMs)
    fstate.add_computation(comp)
    worked = True

    return worked

###########################################
def reg_energy(comp: object, debug: int=0):

    ### For simplicity...
    gmol = comp._job._recipe.subject

    ### Reads output
    comp.read_lines()
    lines = comp.output_lines

    ### Parsing ###
    worked = False
    try: 
        if comp.software == "g16": energy = G16_get_last_energy(lines, debug=debug)
        if comp.software == "qe":  energy = parse_final_energy(lines, debug=debug)
        if energy is not None:     worked = True
    except Exception as exc:
        print(exc)
        worked = False

    ### Storage ###
    if worked: 
        try:
            exists, fstate = find_state(gmol, comp.qc_data.fstate)   ## If exists, it will be updated 
            if not exists: fstate = state(gmol, comp.qc_data.fstate)
            fstate.set_energy(energy, 'au')
            fstate.add_computation(comp)
            worked = True
        except Exception as exc:
            print(exc)
            worked = False

    return worked

###########################################
#def reg_forces(comp: object, debug: int=0):
#    gmol = comp._job._recipe.subject
#
#    ### If not yet available, it reads the output lines
#    if not hasattr(comp,'output_lines'): comp.read_lines()
#    lines = comp.output_lines
#
#    worked = False
#    if comp.software == "qe": forces = QE_get_forces(lines, debug=debug)
#    # Maybe add an assert for the size of the forces vector
#    if len(forces) == len(gmol.natoms): worked = True; comp.forces = np.array(forces)
#    return worked


#def reg_findiff(job: object, debug: int=0):
#
#    ## We collect all computations in a simpler variable
#    comps = job.computations.sort(key=lambda x: (x.index))
#
#    ## Checks that all computations are good and are registered, otherwise quits
#    for idx, c in enumerate(comps):
#        worked = reg_forces(c)
#        if not worked: 
#            print(f"REG FINDIFF: WARNING. Registration of Forces for Computation {idx} didn't work")
#            return None
#
#    worked = False
#    gmol = job._recipe.subject
#    VNMs = get_VNM_from_findiff(job, debug=debug)
#    gmol.VNMs = [vnm for vnm in VNMs]
#    gmol.freqs_cm = [vnm.freq_cm for vnm in VNMs]

##   !!! WE MUST DECIDE WHAT TO DO ABOUT NEGATIVE or NEAR-ZERO FREQUENCIES. IN G16 there is no problem, but here yes
#    if all(vnm.freq >= 0.0 for vnm in VNMs): gmol.isminimum = True
#    else:                                    gmol.isminimum = False
#
#    worked = True
#    if gmol.isminimum:
#        new_coord = G16_get_last_geom(lines, debug=debug)
#        fstate = "min_coord"
#        if hasattr(gmol,fstate): gmol_update_geom(gmol, new_coord, tag=fstate, debug=debug)
#        else:                     gmol_create_geom(gmol, new_coord, tag=fstate, debug=debug)
