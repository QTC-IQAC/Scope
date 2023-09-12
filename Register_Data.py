import sys
import os
import numpy as np
from copy import deepcopy

from Scope.Adapted_from_cell2mol import labels2formula, get_adjmatrix, get_radii, get_blocks, get_molecules, inv
from Scope.Gmol_ops import gmol_update_geom, cell_update_geom, gmol_create_geom, cell_create_geom
from Scope.Parse_QE_outputs import * 
from Scope.Parse_General import read_lines_file, search_string 
from Scope.Parse_G16_outputs import *
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
        line_done, found_done = search_string("JOB DONE", lines, typ='last')
        line_time_start, found_time_start = search_string("Program PWSCF", lines, typ='first')
        line_time_end, found_time_end = search_string("This run was terminated", lines, typ='last')
        line_elapsed, found_elapsed = search_string("PWSCF        :", lines, typ='last')

        if found_elapsed:
            eline = lines[line_elapsed]
            comp.elapsed_time = QE_elapsed_time(eline)
        else: comp.elapsed_time = float(0)

        if found_done and found_time_end: comp.isfinished = True
        else:                             comp.isfinished = False

        ## Other conditions for "isgood" might appear when registering other items
        if found_done and found_time_end: comp.isgood = True
        else:                             comp.isgood = False

        if debug > 1: print(    "REG_GEN: isgood:    ", comp.isgood)
        if debug > 1: print(    "REG_GEN: isfinished:", comp.isfinished)
        if debug > 1: print(    "REG_GEN: found_done:", found_done)
        if debug > 1: print(    "REG_GEN: found_time_start:", found_time_start)
        if debug > 1: print(    "REG_GEN: found_time_end:", found_time_end)
        if debug > 1: print(    "REG_GEN: found_elapsed:", found_elapsed)
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

    gmol = comp._job._recipe.subject

    ### Reads output
    #if not hasattr(comp,'output_lines'): comp.read_lines()
    comp.read_lines()
    lines = comp.output_lines

    ### Defines tag for the new geometry. We avoid blanks
    if ' ' in comp._job.keyword:      new_tag = comp._job.keyword.replace(' ','_')
    else:                             new_tag = comp._job.keyword
    if '-' in comp._job.keyword:      new_tag = comp._job.keyword.replace('-','_')
    else:                             new_tag = comp._job.keyword

    worked = False
    ###############
    ## G09 & G16 ##
    ###############
    if comp.software == "g16": 
        new_coord = G16_get_last_geom(lines, debug=debug)

        ### If the desired tag already exists, then it updates it
        if len(new_coord) > 0:
            if hasattr(gmol,new_tag): 
                gmol_update_geom(gmol, new_coord, tag=new_tag, debug=debug)
                if debug >= 1: print(f"    REG_OPT: geom updated with tag={new_tag}")
            else:                     
                gmol_create_geom(gmol, new_coord, tag=new_tag, debug=debug)
                if debug >= 1: print(f"    REG_OPT: geom created with tag={new_tag}")
            worked = True
        else: 
            print("    REG_OPT: could not extract last coordinates from", comp.out_path)
            if len(lines) > 0: print("    REG_OPT: last output line:", lines[-1])

    ########
    ## QE ##
    ########
    elif comp.software == "qe":

        ## This means that the optimization hasn't finished. 
        ## Even in that case, It will still try to retrieve the last geometry 
        line_BFGS, found_BFGS = search_string("End of BFGS Geometry Optimization", lines, typ='last')
        if not found_BFGS: comp.isgood = False    

        labels, new_coord = parse_final_geometry(lines, debug=debug)
        if len(labels) > 0 and len(new_coord) > 0: 
            try:    factor = gmol.moleclist[0].factor
            except: factor = 1.3

            ### Here we only receive basic information of the molecule
            ### Removed 'cos no molecule information is retained at this stage
            #///////////////////////////
            #///warning, basic_moleclist = get_molecules(labels, new_coord, factor=factor, debug=debug)
            #///assert len(basic_moleclist) == len(gmol.moleclist)
            #///
            #///#### Here we recover the molecule class, but without any data from cell2mol
            #///moleclist = []
            #///for idx, bm in enumerate(basic_moleclist):
            #///    new_molec = molecule(gmol.refcode+str(idx),bm[0], bm[1], bm[2], bm[3])
            #///    moleclist.append(new_molec) 
            #///////////////////////////
            
            ### Updates Molecule or cell, depending on the type of object
            if hasattr(gmol,"moleclist"): 
                if hasattr(gmol,new_tag): cell_update_geom(gmol, new_coord, tag=new_tag, debug=debug)
                else:                     cell_create_geom(gmol, new_coord, tag=new_tag, debug=debug)
            else:                       
                if hasattr(gmol,new_tag): gmol_update_geom(gmol,            new_coord, tag=new_tag, debug=debug)
                else:                     gmol_create_geom(gmol,            new_coord, tag=new_tag, debug=debug)
            worked = True

        else: print("    REG_OPT: empty labels and positions received. Job could not be registered") 
    else: print("    REG_OPT: Software", comp.software, " not considered")

    return worked

###########################################
def reg_frequencies(comp: object, debug: int=0):
    gmol = comp._job._recipe.subject

    ### If not yet available, it reads the output lines
    if not hasattr(comp,'output_lines'): comp.read_lines()
    lines = comp.output_lines

    worked = False
    if comp.software == "g16": 
        ###############
        ## G09 & G16 ##
        ###############
        VNMs = G16_get_VNM(lines, witheigen=True, debug=debug)
        gmol.VNMs = [vnm for vnm in VNMs]
        gmol.freqs_cm = [vnm.freq_cm for vnm in VNMs]
        gmol.Helec = float(G16_get_last_energy(lines))
        if all(vnm.freq >= 0.0 for vnm in VNMs): gmol.isminimum = True
        else:                                    gmol.isminimum = False
            
        worked = True
        if gmol.isminimum:
            new_coord = G16_get_last_geom(lines, debug=debug)
            new_tag = "min_coord"
            if hasattr(gmol,new_tag): gmol_update_geom(gmol, new_coord, tag=new_tag, debug=debug)
            else:                     gmol_create_geom(gmol, new_coord, tag=new_tag, debug=debug)
    return worked

###########################################
def reg_energy(comp: object, debug: int=0):
    gmol = comp._job._recipe.subject

    ### If not yet available, it reads the output lines
    if not hasattr(comp,'output_lines'): comp.read_lines()
    lines = comp.output_lines

    worked = False
    if comp.software == "g16": energy = float(G16_get_last_energy(lines))
    if comp.software == "qe":  energy = float(parse_final_energy(lines, debug=debug))

    if type(energy) == float: worked = True; gmol.Helec = energy
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
#        new_tag = "min_coord"
#        if hasattr(gmol,new_tag): gmol_update_geom(gmol, new_coord, tag=new_tag, debug=debug)
#        else:                     gmol_create_geom(gmol, new_coord, tag=new_tag, debug=debug)
