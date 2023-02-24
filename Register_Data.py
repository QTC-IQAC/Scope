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
def reg_general(job: object, debug: int=0):
   
    if hasattr(job,"output_lines"):  lines = job.output_lines
    else: lines = read_lines_file(job.output_path, flat=True) 

    ###################
    ### Gaussian 16 ###
    ###################
    if job.software.lower() == 'g16' or job.software.lower() == 'g09':
        line_normal, found_normal = search_string("Normal termination", lines, typ='all')
        line_error, found_error   = search_string("Error termination", lines, typ='all')
        #line_route, found_route   = search_string("Route card not found", lines, typ='all')
        line_time, found_time     = search_string("Elapsed time:", lines, typ='last')
        if found_time:
            elapsed_time_list = lines[line_time].split()[2:]
            job.elapsed_time = G16_time_to_sec(elapsed_time_list)
        else: job.elapsed_time = float(0)

        if found_error or found_normal:  job.isfinished = True
        else:                            job.isfinished = False

        if found_normal:                 job.isgood = True
        else:                            job.isgood = False

    #########################
    ### Quantum Espresso ###
    #########################
    elif job.software.lower() == "qe" or job.software.lower() == "quantum_espresso" or job.software.lower() == "quantum espresso":
        line_done, found_done = search_string("JOB DONE", lines, typ='last')
        line_time_start, found_time_start = search_string("Program PWSCF", lines, typ='first')
        line_time_end, found_time_end = search_string("This run was terminated", lines, typ='last')
        line_elapsed, found_elapsed = search_string("PWSCF        :", lines, typ='last')

        if found_elapsed:
            eline = lines[line_elapsed]
            job.elapsed_time = QE_elapsed_time(eline)
        else: job.elapsed_time = float(0)

        if found_done and found_time_end: job.isfinished = True
        else:                             job.isfinished = False

        if found_done and found_time_end: job.isgood = True
        else:                             job.isgood = False

        if debug >= 1: print(    "REG_GEN: isgood:    ", job.isgood)
        if debug >= 1: print(    "REG_GEN: isfinished:", job.isfinished)
        if debug >= 1: print(    "REG_GEN: found_done:", found_done)
        if debug >= 1: print(    "REG_GEN: found_time_start:", found_time_start)
        if debug >= 1: print(    "REG_GEN: found_time_end:", found_time_end)
        if debug >= 1: print(    "REG_GEN: found_elapsed:", found_elapsed)
    #############
    ### Other ###
    #############
    else:
        job.isfinished   = False
        job.isgood       = False
        job.elapsed_time = float(0)
        print(f"    REG_GENERAL: Registry of {job.software} jobs is not implemented.")

###########################################
def reg_optimization(gmol: object, job: object, debug: int=0):
    if hasattr(job,"output_lines"):  lines = job.output_lines
    else: lines = read_lines_file(job.output_path, flat=True) 

    ### Defines tag for the new geometry. We avoid blanks
    if ' ' in job.code: new_tag = job.code.replace(' ','_')
    else:               new_tag = job.code

    worked = False
    ###############
    ## G09 & G16 ##
    ###############
    if job.software.lower() == "g16" or job.software.lower() == "g09":
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
        else: print("    REG_OPT: could not extract last coordinates from", job.output_path)

    ########
    ## QE ##
    ########
    elif job.software.lower() == "qe" or job.software.lower() == "quantum_espresso" or job.software.lower() == "quantum espresso":
        labels, new_coord = parse_final_geometry(lines, debug=debug)
        if len(labels) > 0 and len(new_coord) > 0: 
            try:    factor = gmol.moleclist[0].factor
            except: factor = 1.3

            ### Here we only receive basic information of the molecule
            warning, basic_moleclist = get_molecules(labels, new_coord, factor=factor, debug=debug)
            assert len(basic_moleclist) == len(gmol.moleclist)
            
            #### Here we recover the molecule class, but without any "charge" data 
            #moleclist = []
            #for idx, bm in enumerate(basic_moleclist):
            #    new_molec = molecule(gmol.refcode+str(idx),bm[0], bm[1], bm[2], bm[3])
            #    moleclist.append(new_molec) 
            
            #if not warning and debug >= 1: print(f"UPDATE_COORDINATES: found {len(moleclist)} molecules in cell")

            ### Updates Molecule or cell, depending on the type of object
            if hasattr(gmol,"moleclist"): 
#                if hasattr(gmol,new_tag): cell_update_geom(gmol, moleclist, new_coord, tag=new_tag, debug=debug)
#                else:                     cell_create_geom(gmol, moleclist, new_coord, tag=new_tag, debug=debug)
                if hasattr(gmol,new_tag): cell_update_geom(gmol, new_coord, tag=new_tag, debug=debug)
                else:                     cell_create_geom(gmol, new_coord, tag=new_tag, debug=debug)
            else:                       
                if hasattr(gmol,new_tag): gmol_update_geom(gmol,            new_coord, tag=new_tag, debug=debug)
                else:                     gmol_create_geom(gmol,            new_coord, tag=new_tag, debug=debug)
            worked = True

        else: print("    REG_OPT: empty labels and positions received. Job could not be registered") 
    else: print("    REG_OPT: Software", job.software, " not considered")

    return worked

###########################################
def reg_frequencies(gmol: object, job: object, debug: int=0):
    if hasattr(job,"output_lines"):  lines = job.output_lines
    else: lines = read_lines_file(job.output_path, flat=True) 

    worked = False
    ###############
    ## G09 & G16 ##
    ###############
    if job.software.lower() == "g16" or "g09":
        VNMs = G16_get_VNM(lines, witheigen=True, debug=debug)
        gmol.VNMs = [vnm for vnm in VNMs]
        gmol.freqs_cm = [vnm.freq_cm for vnm in VNMs]
        gmol.Helec = G16_get_last_energy(lines)
        if all(vnm.freq >= 0.0 for vnm in VNMs): gmol.isminimum = True
        else:                                    gmol.isminimum = False
            
        worked = True
        if gmol.isminimum:
            new_coord = G16_get_last_geom(lines, debug=debug)
            ### Defines tag for the new geometry. We avoid blanks
            new_tag = "min_coord"
            ### If the desired tag already exists, then it updates it
            if hasattr(gmol,new_tag): 
                gmol_update_geom(gmol, new_coord, tag=new_tag, debug=debug)
            else:                     
                gmol_create_geom(gmol, new_coord, tag=new_tag, debug=debug)

    ########
    ## QE ##
    ########
    #elif job.software.lower() == "qe" or job.software.lower() == "quantum_espresso" or job.software.lower() == "quantum espresso":

    #else: print("    REG_OPTIMIZATIONS: Software", job.software, " not considered")

    return worked
