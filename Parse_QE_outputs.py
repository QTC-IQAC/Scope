import os
import numpy as np
import sys

from Scope.Parse_General import search_string, read_lines_file
from Scope.unit_cell_tools import cellvec_2_cellparam, get_unit_cell_volume
from Scope.Scope_Classes import periodic_xyz

bohr2angs = 0.529177

def get_cell_vectors(lines, debug: int=0):
    cellvec_lines = []
    cellvec_strings = ["celldm(1)", "celldm(4)", "crystal axes", "unit-cell volume"]
    for sdx, s in enumerate(cellvec_strings):
        cellvec_lines.append(search_string(s, lines, typ="last")[0])
    celldim = []
    try:
        celldim.append(float(lines[cellvec_lines[0]].split()[1]))
        celldim.append(float(lines[cellvec_lines[0]].split()[3]))
        celldim.append(float(lines[cellvec_lines[0]].split()[5]))
        celldim.append(float(lines[cellvec_lines[1]].split()[1]))
        celldim.append(float(lines[cellvec_lines[1]].split()[3]))
        celldim.append(float(lines[cellvec_lines[1]].split()[5]))
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
    return cellvec, celldim, cellparam
    
def parse_final_geoopt_step(lines: str):
    last_step_init = "Self-consistent Calculation"
    last_step_last = "End final coordinates"
    last_step_init_line, dummy = search_string(last_step_init, lines, typ="last")
    last_step_last_line, dummy = search_string(last_step_last, lines, typ="last")
    return last_step_init_line, last_step_last_line 

def check_job_requirements(lines: str, key_str: list=["JOB DONE", "Begin final coordinates", "End final coordinates"]):
    item = []
    found = []
    for sx, string in enumerate(str):
        search, found = search_string(string, lines, typ="last")
        line.append(search)
        found.append(found)
    return line, found
    
def parse_final_geometry(lines: str):
    key_line, key_found = check_job_requirements(lines)
    last_step_init_line, last_step_last_line = parse_final_geoopt_step(lines)
    last_step_lines = lines[last_step_init_line:last_step_last_line]

    if all(f for f in key_found):  ## Everything worked
        last_geo_init_string = "ATOMIC_POSITIONS"
        last_geo_init_line = search_string(last_geo_init_string, last_step_lines, typ="last")[0] + 1
        last_geo_last_line = last_step_end_line-1
    elif key_found[0] and not key_found[3]:     ### Job ended but did not converge
        last_geo_init_string = "ATOMIC_POSITIONS"
        last_geo_last_string = "Writing output data"
        last_geo_init_line = search_string(last_geo_init_string, last_step_lines, typ="last")[0] + 1
        last_geo_last_line = search_string(last_geo_last_string, last_step_lines, typ="first")[0] - 4
    else: 
        warning = True
 
    last_geo_lines = last_step_lines[last_geo_init_line:last_geo_last_line]
    if not warning:
        ### Reads units of coordinates from the ATOMIC POSITIONS line
        if "angstrom" in str(lines[last_geo_init_line - 1].lower()): units = 'angstrom'
        elif "bohr" in str(lines[last_geo_init_line - 1].lower()): units = 'bohr'
        pos = []
        labels = []
        for idx, l in enumerate(last_geo_lines):
            line_data = l.split()
            if len(line_data) == 4:
                label, x, y, z = l.split()
                if units == "angstrom": pos.append([float(x), float(y), float(z)])
                elif units == "bohr": pos.append([float(x*bohr2angs), float(y*bohr2angs), float(z*bohr2angs)])
                labels.append(label)
    return labels, pos

def parse_final_energy(lines: str):
    string="Final energy   ="
    linenum = search_string(string, lines, typ="last")[0]
    val = lines[linenum].split()[2]
    return val 

def parse_hubbard_energy(lines: str):
    string = "Hubbard energy            ="
    linenum = search_string(string, lines, typ="last")[0]
    val = lines[linenum].split()[2]
    return val 

def parse_grimme_energy(lines: str):
    string = "DFT-D3 Dispersion         ="
    linenum = search_string(string, lines, typ="last")[0]
    val = lines[linenum].split()[2]
    return val 

def parse_total_energy(lines: str):
    string = "!    total energy"
    linenum = search_string(string, lines, typ="last")[0]
    val = lines[linenum].split()[2]
    return val

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
