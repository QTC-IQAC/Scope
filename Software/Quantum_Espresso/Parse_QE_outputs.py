import os
import numpy as np
import sys

from Scope.Parse_General import search_string, read_lines_file
from Scope.Unit_cell_tools import cellvec_2_cellparam, get_unit_cell_volume
from Scope.Classes_QC import periodic_xyz
import Scope.Constants

bohr2angs = Scope.Constants.bohr2angs

def parse_geometry_from_step(lines, debug=0):
    
    bohr2angs = Scope.Constants.bohr2angs
    
    init_str = "ATOMIC_POSITIONS"
    last_str = "Writing output data"
    init_linenum, found_init = search_string(init_str, lines, typ="first")
    last_linenum, found_end  = search_string(last_str, lines, typ="first")
    if not found_init or not found_end: return None, None
        
    last_linenum = last_linenum - 3   ## Writing output data file is written 3 lines after
    #print(init_linenum, last_linenum)
    geo_lines = lines[init_linenum:last_linenum]

    ### Reads units of coordinates from the ATOMIC POSITIONS line
    at_pos_line=geo_lines[0].lower()
    if "angstrom" in at_pos_line: units = 'angstrom'
    elif "bohr"   in at_pos_line: units = 'bohr'
    else: print("PARSE_GEOMETRY: error reading units of coordinates:", at_pos_line)

    ### Reads Coordinates
    coord = []
    labels = []
    for idx, l in enumerate(geo_lines[1:]):
        line_data = l.split()
        if len(line_data) == 4:
            label, x, y, z = l.split()
            try:
                if units == "angstrom": coord.append([float(x), float(y), float(z)])
                elif units == "bohr":   coord.append([float(x*bohr2angs), float(y*bohr2angs), float(z*bohr2angs)])
                labels.append(label)
            except exception as exc:
                if debug >= 1: print(f"PARSE_GEOMETRY: exception {exc} reading line: {l}")
        else:
            if debug >= 1: print(f"PARSE_GEOMETRY: line discarded: {line_data}")
    return labels, coord

def parse_forces_from_step(lines, debug=0):
    init_str = "Forces acting on atoms"
    last_str = "Total force ="
    init_linenum, found_init = search_string(init_str, lines, typ="first")
    last_linenum, found_end  = search_string(last_str, lines, typ="first")
    if not found_init or not found_end: return None
    
    init_linenum = init_linenum + 2
    last_linenum = last_linenum - 1   
    force_lines = lines[init_linenum:last_linenum]

    ### Reads Atomic Forces
    forces = []
    for idx, l in enumerate(force_lines):
        line_data = l.split()
        if len(line_data) == 9:
            dummy, atnum, dummy, specnum, dummy, dummy, x, y, z = l.split()
            try: 
                forces.append([float(x), float(y), float(z)])
            except exception as exc:
                if debug >= 1: print(f"PARSE_FORCES: exception {exc} reading line: {l}")
        else:
            if debug >= 1: print(f"PARSE_FORCES: line discarded: {line_data}")
                
    ### Reads Total Force
    total_force = lines[last_linenum+1].split()[3]
        
    return forces, total_force



def check_job_requirements(lines: list, key_str: list=["JOB DONE", "Begin final coordinates", "End final coordinates"], debug: int=0):
    item = []
    bools = []
    for sx, string in enumerate(key_str):
        search, found = search_string(string, lines, typ="last")
        item.append(search)
        bools.append(found)
    return item, bools

def check_maxseconds(lines: list, string: str="Maximum CPU time exceeded", debug: int=0):
    linenum, found = search_string(string, lines, typ="last")
    return linenum, found

def check_convergence(lines: list, string: str="convergence NOT achieved after", debug: int=0):
    linenum, found = search_string(string, lines, typ="last")
    return linenum, found
    
def parse_cell_vectors(lines, debug: int=0):
    string = "CELL_PARAMETERS"
    linenum, found = search_string(string, lines, typ="first")
    if not found: return None
    cellvec = []
    cellvec.append(list(float(i) for i in lines[linenum+1].split()))
    cellvec.append(list(float(i) for i in lines[linenum+2].split()))
    cellvec.append(list(float(i) for i in lines[linenum+3].split()))
    return cellvec

def get_cell_vectors(lines, debug: int=0):
    cellvec_lines = []
    cellvec_strings = ["celldm(1)", "celldm(4)", "crystal axes", "unit-cell volume"]
    for sdx, s in enumerate(cellvec_strings):
        cellvec_lines.append(search_string(s, lines, typ="last")[0])
    celldim = []
    try:
        celldim.append(float(lines[cellvec_lines[0]].split()[1])*bohr2angs)
        celldim.append(float(lines[cellvec_lines[0]].split()[3])*bohr2angs)
        celldim.append(float(lines[cellvec_lines[0]].split()[5])*bohr2angs)
        celldim.append(float(lines[cellvec_lines[1]].split()[1])*bohr2angs)
        celldim.append(float(lines[cellvec_lines[1]].split()[3])*bohr2angs)
        celldim.append(float(lines[cellvec_lines[1]].split()[5])*bohr2angs)
    except Exception as exc:
        print("Error trying to parse cell dimensions from file")
        print("Exception is:", exc)
        print("cellvec_lines:", cellvec_lines)
    cellvec = []
    v1 = np.array(lines[cellvec_lines[2]+1].split("=")[1].replace(")", "").replace("(", "").split()).astype(float)
    v2 = np.array(lines[cellvec_lines[2]+2].split("=")[1].replace(")", "").replace("(", "").split()).astype(float)
    v3 = np.array(lines[cellvec_lines[2]+3].split("=")[1].replace(")", "").replace("(", "").split()).astype(float)
    cellvec.append(v1*celldim[0])
    cellvec.append(v2*celldim[0])
    cellvec.append(v3*celldim[0])
    cellparam = cellvec_2_cellparam(cellvec)
    return np.array(cellvec), celldim, cellparam

def parse_final_geoopt_step(lines, debug: int=0):
    last_step_init = "Self-consistent Calculation"
    last_step_last = "End final coordinates"
    last_step_init_line, found1 = search_string(last_step_init, lines, typ="last")
    last_step_last_line, found2 = search_string(last_step_last, lines, typ="last")
    max_line, max_found = check_maxseconds(lines, debug=debug)
 
    ## In some cases, "End final coordinates" is not written. Probably because the geometry has not fished the optimization 
    ## Then, as an alternative, one can search for NEW-OLD 
    if found1 and found2 and last_step_init_line > last_step_last_line: 
        if debug >= 1: print(f"    PARSE_FINAL_GEOOPT_STEP: 'End final coordinates' is before 'Self-consistent'. It must be a vc-relax")
        if debug >= 1: print(f"    PARSE_FINAL_GEOOPT_STEP: Going for the alternative. Limiting to uplim={last_step_last_line}")
        last_step_init_line, found1 = search_string(last_step_init, lines, typ="last", uplim=last_step_last_line)
        if debug >= 1: print(f"    PARSE_FINAL_GEOOPT_STEP: Results: {last_step_init_line}, {found1}") 

    if not found2: 
        if debug >= 1: print(f"    PARSE_FINAL_GEOOPT_STEP: 'End final coordinates' not found")
        if debug >= 1: print(f"    PARSE_FINAL_GEOOPT_STEP: Going for the alternative")
        last_step_last = "Writing output data"   
        if max_found: last_step_last_line, found2 = search_string(last_step_last, lines, typ="all"); last_step_last_line = last_step_last_line[-2]
        else:         last_step_last_line, found2 = search_string(last_step_last, lines, typ="last")

    if last_step_init_line > last_step_last_line:
        last_step_init_line, found1 = search_string(last_step_init, lines, typ="last", uplim=last_step_last_line)

    if found1 and found2: worked = True
    else: worked = False
    return worked, last_step_init_line, last_step_last_line 
    
def parse_final_geometry(lines, debug: int=0):
    key_line, key_found = check_job_requirements(lines, debug=debug)
    worked, last_step_init_line, last_step_last_line = parse_final_geoopt_step(lines, debug=debug)
    if not worked: 
        print("PARSE_FINAL_GEOMETRY: final_geoopt_step didn't work. Retrieving empty labels and coordinates. It will fail")
        return [], []
    else:
        last_step_lines = lines[last_step_init_line:last_step_last_line+1]
        if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: final step lines between {last_step_init_line+1} and {last_step_last_line+1}")
        if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: first last_step_line: {last_step_lines[0]}")
        if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: last  last_step_line: {last_step_lines[-1]}")
    
        warning = False
        if all(f for f in key_found):               ### Everything worked
            last_geo_init_string = "ATOMIC_POSITIONS"
            last_geo_init_line = search_string(last_geo_init_string, last_step_lines, typ="last")[0] 
            last_geo_last_line = last_step_last_line-1
            if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: last geo lines between {last_step_init_line+last_geo_init_line+1} and {last_geo_last_line+1}")
        elif key_found[0] and not key_found[2]:     ### Job ended but did not converge
            last_geo_init_string = "ATOMIC_POSITIONS"
            last_geo_last_string = "Writing output data"
            last_geo_init_line = search_string(last_geo_init_string, last_step_lines, typ="last")[0] 
            last_geo_last_line = search_string(last_geo_last_string, last_step_lines, typ="first")[0] - 3
            if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: last geo lines between {last_step_init_line+last_geo_init_line+1} and {last_step_init_line+last_geo_last_line+1}")
        elif all(not f for f in key_found):         ### Nothing, probably a job terminated due to time limits
            last_geo_init_string = "ATOMIC_POSITIONS"
            last_geo_last_string = "Writing output data"
            last_geo_init_line = search_string(last_geo_init_string, last_step_lines, typ="last", debug=debug)[0] 
            last_geo_last_line = search_string(last_geo_last_string, last_step_lines, typ="first", debug=debug)[0] - 3
            if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: last geo lines between {last_step_init_line+last_geo_init_line+1} and {last_step_init_line+last_geo_last_line+1}")
        else: 
            print("no key was found", key_found)
            warning = True
     
        last_geo_lines = last_step_lines[last_geo_init_line:last_geo_last_line]
        if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: first last_geo_line: {last_geo_lines[0]}")
        if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: last  last_geo_line: {last_geo_lines[-1]}")
        if not warning:
            ### Reads units of coordinates from the ATOMIC POSITIONS line
            at_pos_line=last_geo_lines[0].lower()
            if "angstrom" in at_pos_line: units = 'angstrom'
            elif "bohr" in at_pos_line: units = 'bohr'
            else: print("PARSE_FINAL_GEOMETRY: error reading units of coordinates:", at_pos_line) 
            if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: found coordinates in", units)

            coord = []
            labels = []
            for idx, l in enumerate(last_geo_lines[1:]):
                line_data = l.split()
                if len(line_data) == 4:
                    label, x, y, z = l.split()
                    try: 
                        if units == "angstrom": coord.append([float(x), float(y), float(z)])
                        elif units == "bohr": coord.append([float(x*bohr2angs), float(y*bohr2angs), float(z*bohr2angs)])
                        labels.append(label)
                    except exception as exc:
                        if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: exception {exc} reading line: {l}")   
                else: 
                    if debug >= 1: print(f"PARSE_FINAL_GEOMETRY: line discarded: {line_data}")
            #print("coord[0]=",  coord[0])
            #print("coord[-1]=", coord[-1])
        return labels, coord

def parse_final_energy(lines, debug: int=0):
    #line_err, found_conv_error = check_convergence(lines)
    #if found_conv_error: return None 
    string="Final energy   ="
    linenum, found = search_string(string, lines, typ="last")
    if found: val = float(lines[linenum].split()[3])
    else:
        string = "!    total energy              ="
        linenum, found  = search_string(string, lines, typ="last")
        if found: val = float(lines[linenum].split()[4])
        else: return None
    return val
#    return val*Scope.Constants.ry2har 

def parse_hubbard_energy(lines, debug: int=0):
    string = "Hubbard energy            ="
    linenum = search_string(string, lines, typ="last")[0]
    val = lines[linenum].split()[2]
    return val 

def parse_grimme_energy(lines, debug: int=0):
    string = "DFT-D3 Dispersion         ="
    linenum = search_string(string, lines, typ="last")[0]
    val = lines[linenum].split()[2]
    return val 

def parse_total_energy(lines, debug: int=0):
    string = "!    total energy"
    linenum = search_string(string, lines, typ="last")[0]
    val = lines[linenum].split()[2]
    return val

def QE_elapsed_time(eline: str, debug: int=0):
    time = eline.split("WALL")[-2]
    time2 = time.split("CPU")[1].rstrip().lstrip()

    rest = time2
    days = float(0)
    hours = float(0)
    minutes = float(0)
    seconds = float(0)
    total = float(0)
    if len(rest) > 0:
        if "d" in rest:
            spl = rest.split("d")
            days = float(spl[0])
            rest = str(spl[1])
        if "h" in rest:
            spl = rest.split("h")
            hours = float(spl[0])
            rest = str(spl[1])
        if "m" in rest:
            spl = rest.split("m")
            minutes = float(spl[0])
            rest = str(spl[1])
        if "s" in rest:
            spl = rest.split("s")
            seconds = float(spl[0])
            rest = str(spl[1])

    total = days*86400 + hours*3600 + minutes*60 + seconds
    return total 


#def QE_get_forces():
#    WARNING: forces are read in Ry/au in QE. Must convert to Hartre/au




#
#        new_periodic_xyz.add_energy(energy)
#
#
#
#
#        new_periodic_xyz = periodic_xyz(name, index, labels, pos, cellparam)
#        cellvec, celldim, cellparam = get_cell_vectors(lines)
#
#        #### Creates object of class "periodic_xyz"
#        new_periodic_xyz.add_cellvec(cellvec, celldim)
#        return new_periodic_xyz, warning 
#    else: 
#        new_periodic_xyz = periodic_xyz("Empty", 0, [], [], np.zeros((6)))
#        return new_periodic_xyz, warning
#       
