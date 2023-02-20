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

    ###############
    ## G09 & G16 ##
    ###############
    if job.software.lower() == "g16" or job.software.lower() == "g09":
        if debug >= 1: print("    REG_OPTIMIZATION received", job.software, "and identified it as G16 or G09")
        new_coord = G16_get_last_geom(lines, debug=debug)

        ### If the desired tag already exists, then it updates it
        if hasattr(gmol,new_tag): 
            gmol_update_geom(gmol, new_coord, tag=new_tag, debug=debug)
            if debug >= 1: print(f"    REG_OPT: geom updated with tag={new_tag}")
        else:                     
            gmol_create_geom(gmol, new_coord, tag=new_tag, debug=debug)
            if debug >= 1: print(f"    REG_OPT: geom created with tag={new_tag}")

    ########
    ## QE ##
    ########
    elif job.software.lower() == "qe" or job.software.lower() == "quantum_espresso" or job.software.lower() == "quantum espresso":
        if debug >= 1: print("    REG_OPTIMIZATION received", job.software, "and identified it as QUANTUM ESPRESSO")
        labels, pos = parse_final_geometry(lines, debug=debug)
        if len(labels) > 0 and len(pos) > 0: 
            try:    factor = gmol.moleclist[0].factor
            except: factor = 1.3
            warning, moleclist = get_molecules(labels, pos, factor=factor, debug=debug)
            if not warning and debug >= 1: print(f"UPDATE_COORDINATES: found {len(moleclist)} molecules in cell")

            ### Updates Molecule or cell, depending on the type of object
            if hasattr(gmol,"moleclist"): 
                if hasattr(gmol,new_tag): cell_update_geom(gmol, moleclist, new_coord, tag=new_tag, debug=debug)
                else:                     cell_create_geom(gmol, moleclist, new_coord, tag=new_tag, debug=debug)
            else:                       
                if hasattr(gmol,new_tag): gmol_update_geom(gmol,            new_coord, tag=new_tag, debug=debug)
                else:                     gmol_create_geom(gmol,            new_coord, tag=new_tag, debug=debug)

        else: print("    REG_OPTIMIZATIONS: empty labels and positions received. Job could not be registered") 
    else: print("    REG_OPTIMIZATIONS: Software", job.software, " not considered")

###########################################
def reg_frequencies(gmol: object, job: object, debug: int=0):
    if hasattr(job,"output_lines"):  lines = job.output_lines
    else: lines = read_lines_file(job.output_path, flat=True) 

    ###############
    ## G09 & G16 ##
    ###############
    if job.software.lower() == "g16" or "g09":
        VNMs = G16_get_VNM(lines, witheigen=True, debug=debug)
        gmol.VNMs = [vnm for vnm in VNMs]
        gmol.freqs_cm = [vnm.freq_cm for vnm in VNMs]
        gmol.Helec = G16_get_last_energy(lines)
        if all(vnm.freq >= 0.0 for vnm in VNMs):
            gmol.isminimum = True
        else:
            gmol.isminimum = False

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
