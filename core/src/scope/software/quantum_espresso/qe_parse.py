from   scope.parse_general   import search_string
from   scope.geometry        import cellvec_2_cellparam
import scope.constants

bohr2angs = scope.constants.bohr2angs
ry2har    = scope.constants.ry2har

###############
### PARSING ###
###############
def parse_findiff_input(lines: list, natoms: int, debug: int=0):
    """Reads a SCOPE-generated periodic QE displacement input.

    Returns elemental labels, Angstrom coordinates, atomic masses in amu and
    Angstrom cell vectors, in input atom order. Requires ATOMIC_SPECIES and
    Angstrom ATOMIC_POSITIONS/CELL_PARAMETERS cards; malformed inputs raise
    ValueError. Numeric spin suffixes on SCOPE species names are removed.
    """
    species_idx, found_species = search_string('ATOMIC_SPECIES', lines, typ='first')
    coord_idx, found_coord     = search_string('ATOMIC_POSITIONS', lines, typ='first')
    cell_idx, found_cell       = search_string('CELL_PARAMETERS', lines, typ='first')
    if not (found_species and found_coord and found_cell): raise ValueError("PARSE_FINDIFF_INPUT: Missing species, positions or cell card")
    if 'angstrom' not in lines[coord_idx].lower() or 'angstrom' not in lines[cell_idx].lower(): raise ValueError("PARSE_FINDIFF_INPUT: Angstrom positions and cell vectors are required")

    species = {}
    for line in lines[species_idx+1:]:
        fields = line.split()
        if len(fields) != 3: break
        try: mass = float(fields[1])
        except ValueError: break
        species[fields[0]] = mass
    try:
        atoms  = [line.split() for line in lines[coord_idx+1:coord_idx+1+natoms]]
        labels = [atom[0].rstrip('0123456789.+-') for atom in atoms]
        coord  = [[float(value) for value in atom[1:4]] for atom in atoms]
        masses = [species[atom[0]] for atom in atoms]
        cell   = [[float(value) for value in line.split()] for line in lines[cell_idx+1:cell_idx+4]]
    except (KeyError, IndexError, ValueError) as exc:
        raise ValueError("PARSE_FINDIFF_INPUT: Could not read displacement geometry or masses") from exc
    if len(atoms) != natoms or any(len(row) != 3 for row in coord) or len(cell) != 3 or any(len(row) != 3 for row in cell): raise ValueError("PARSE_FINDIFF_INPUT: Incomplete positions or cell vectors")
    if debug > 0: print(f"PARSE_FINDIFF_INPUT: {natoms} atoms, positions in Angstrom, masses in amu")
    return labels, coord, masses, cell


def parse_start_scf(lines, debug=0):
    start_SCF, found  = search_string("Self-consistent Calculation", lines)
    return start_SCF, found

def parse_end_scf(lines, debug=0):
    end_SCF, found  = search_string("End of self-consistent calculation", lines)
    return end_SCF, found

def parse_scf_accuracies_from_step(lines: list, debug: int=0) -> list:
    """Returns the estimated SCF accuracies reported within one QE step."""
    linenums, found = search_string("estimated scf accuracy", lines)
    if not found: return []

    accuracies = []
    for linenum in linenums:
        try:                  accuracies.append(float(lines[linenum].split()[-2]))
        except (ValueError, IndexError):
            if debug > 0: print(f"PARSE_SCF_ACCURACIES: could not read line: {lines[linenum].rstrip()}")
    return accuracies

def parse_geometry_from_step(lines, debug=0):
    init_str = "ATOMIC_POSITIONS"
    last_str = "output data"
    init_linenum, found_init = search_string(init_str, lines, typ="first")
    last_linenum, found_end  = search_string(last_str, lines, typ="first")
    if not found_init or not found_end: return None, None
        
    last_linenum = last_linenum - 3   ## output data file is written 3 lines after
    #print(init_linenum, last_linenum)
    geo_lines = lines[init_linenum:last_linenum]

    ### Reads units of coordinates from the ATOMIC POSITIONS line
    at_pos_line=geo_lines[0].lower()
    if "angstrom" in at_pos_line: units = 'angstrom'
    elif "bohr"   in at_pos_line: units = 'bohr'
    else: print("PARSE_GEOMETRY: error reading units of coordinates:", at_pos_line)

    ### Reads Coordinates
    coord = []
    tmp_labels = []
    for idx, l in enumerate(geo_lines[1:]):
        line_data = l.split()
        if len(line_data) == 4:
            label, x, y, z = l.split()
            try:
                if units == "angstrom": coord.append([float(x), float(y), float(z)])
                elif units == "bohr":   coord.append([float(x*bohr2angs), float(y*bohr2angs), float(z*bohr2angs)])
                tmp_labels.append(label)
            except exception as exc:
                if debug >= 1: print(f"PARSE_GEOMETRY: exception {exc} reading line: {l}")
        else:
            if debug >= 1: print(f"PARSE_GEOMETRY: line discarded: {line_data}")
 
    ## Originally, labels include digits to follow the spin state of metal atoms. For instance 'Fe4' indicates a HS Fe atom
    ## These digits must be removed from labels when storing the data, since the digit is only for QE
    labels = []
    for l in tmp_labels:
        labels.append(str(''.join([c for c in l if not c.isdigit()])))

    return labels, coord

def parse_forces_from_step(lines, debug: int=0):
    init_str = "Forces acting on atoms"
    last_str = "Total force ="
    init_linenum, found_init = search_string(init_str, lines, typ="first")
    last_linenum, found_end  = search_string(last_str, lines, typ="first")
    if debug >= 1: print(f"PARSE_FORCES: {found_init=} {found_end=}")
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
            except ValueError as exc:
                if debug >= 1: print(f"PARSE_FORCES: exception {exc} reading line: {l}")
                return None
        else:
            if debug >= 1: print(f"PARSE_FORCES: line discarded: {line_data}")
    return forces
                
def parse_total_force_from_step(lines, debug: int=0):
    init_str = "Forces acting on atoms"
    last_str = "Total force ="
    init_linenum, found_init = search_string(init_str, lines, typ="first")
    last_linenum, found_end  = search_string(last_str, lines, typ="first")
    if not found_init or not found_end: return None

    init_linenum = init_linenum + 2
    last_linenum = last_linenum - 1   
    force_lines = lines[init_linenum:last_linenum]

    ### Reads Total Force
    total_force = lines[last_linenum+1].split()[3]
    return total_force

#def parse_cell_vectors(lines, debug: int=0):
#    string = "CELL_PARAMETERS"
#    linenum, found = search_string(string, lines, typ="first")
#    if not found: return None
#    cellvec = []
#    cellvec.append(list(float(i) for i in lines[linenum+1].split()))
#    cellvec.append(list(float(i) for i in lines[linenum+2].split()))
#    cellvec.append(list(float(i) for i in lines[linenum+3].split()))
#    return cellvec

def parse_energy_from_step(lines, debug: int=0):
    string="Final energy   ="
    linenum, found = search_string(string, lines, typ="last")
    if found: val = float(lines[linenum].split()[3])
    else:
        string = "!    total energy              ="
        linenum, found  = search_string(string, lines, typ="last")
        if found: val = float(lines[linenum].split()[4])
        else: return None
    return val * ry2har

#def parse_hubbard_energy(lines, debug: int=0):
#    string = "Hubbard energy            ="
#    linenum = search_string(string, lines, typ="last")[0]
#    val = lines[linenum].split()[2]
#    return val * ry2har
#
#def parse_grimme_energy(lines, debug: int=0):
#    string = "DFT-D3 Dispersion         ="
#    linenum = search_string(string, lines, typ="last")[0]
#    val = lines[linenum].split()[2]
#    return val * ry2har
#
#def parse_total_energy(lines, debug: int=0):
#    string = "!    total energy"
#    linenum = search_string(string, lines, typ="last")[0]
#    val = lines[linenum].split()[2]
#    return val * ry2har

def parse_makov_from_step(lines, debug: int=0):
    string = "!    Total+Makov-Payne energy"
    linenum, found = search_string(string, lines, typ="last")
    if found: val = float(lines[linenum].split()[4]); 
    else: return None
    return val * ry2har

def parse_cpu_time(lines: str, debug: int=0):
    # To be implemented
    return float(0)

def parse_elapsed_time(lines: str, debug: int=0):
    line_elapsed, found_elapsed = search_string("PWSCF        :", lines, typ='last')
    if not found_elapsed: return float(0)

    eline = lines[line_elapsed]

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

##############
### STATUS ###
##############
def parse_status_finished(lines):
    search, found = search_string("JOB DONE", lines, typ="last")
    return found

def parse_opt_status(lines):
    """Returns the latest explicit BFGS outcome, or 'not finished' if none is printed."""
    # 1) Locate convergence, step exhaustion and optimizer failure messages
    converged_line, found_converged = search_string('bfgs converged in', lines, typ='last')
    exhausted_line, found_exhausted = search_string('The maximum number of steps has been reached', lines, typ='last')
    stuck_line,     found_stuck     = search_string('history already reset at previous step:', lines, typ='last')
    if not found_converged: converged_line = -1
    if not found_exhausted: exhausted_line = -1
    if not found_stuck:     stuck_line     = -1

    # 2) Return the outcome of the last message in the file
    if converged_line > max(exhausted_line, stuck_line): return 'finished'
    if exhausted_line > stuck_line:                     return 'exhausted'
    if found_stuck:                                     return 'stucked'
    return 'not finished'

def parse_scf_status(lines):
    linenum1, found1 = search_string("convergence has been achieved", lines, typ="last")
    linenum2, found2 = search_string("convergence NOT achieved",      lines, typ="last")
    if found1 and (not found2 or linenum1 > linenum2): return True
    elif found2:                                    return False
    else:                                 return None

def parse_coord_status(lines):
    linenum1, found1 = search_string("ATOMIC_POSITIONS"   , lines, typ="last")
    linenum2, found2 = search_string("output data",   lines, typ="last")
    if found1 and found2:                 return True
    else:                                 return False

def parse_timelimit_status(lines: list, debug: int=0):
    linenum, found = search_string("Maximum CPU time exceeded", lines, typ="last")
    if found:                             return True
    else:                                 return False

def parse_volume_from_step(lines, debug: int=0):
    string = "unit-cell volume"
    linenum, found = search_string(string, lines, typ="first")
    if found:
        if   len(lines[linenum].split()) == 10: val = float(lines[linenum].split()[4])
        elif len(lines[linenum].split()) == 5:  val = float(lines[linenum].split()[3])
        else: print(len(lines[linenum].split()), "as length for volume line")
    else: val = None 
    return val

def parse_density_from_step(lines, debug: int=0):
    string = "density ="
    linenum, found = search_string(string, lines, typ="first")
    if found: val = float(lines[linenum].split()[2])
    else:     val = None
    return val

def parse_cell_parameters(lines, debug: int=0):
    import numpy as np
    string = "CELL_PARAMETERS (angstrom)"
    linenum, found = search_string(string, lines, typ="last")
    if not found: return None, None
    cellvec = []
    cellvec.append(list(map(float,lines[linenum+1].split()[0:3])))
    cellvec.append(list(map(float,lines[linenum+2].split()[0:3])))
    cellvec.append(list(map(float,lines[linenum+3].split()[0:3])))
    cellparam = cellvec_2_cellparam(cellvec)
    return np.array(cellvec), cellparam

def parse_cell_vectors(lines, debug: int=0):
    import numpy as np
    cellvec_lines = []
    cellvec_strings = ["celldm(1)", "celldm(4)", "crystal axes"] #, "unit-cell volume"]

    for sdx, s in enumerate(cellvec_strings):
        linenum, found = search_string(s, lines, typ="last")
        if not found: return None, None, None
        cellvec_lines.append(linenum)
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
    try:
        v1 = np.array(lines[cellvec_lines[2]+1].split("=")[1].replace(")", "").replace("(", "").split()).astype(float)
        v2 = np.array(lines[cellvec_lines[2]+2].split("=")[1].replace(")", "").replace("(", "").split()).astype(float)
        v3 = np.array(lines[cellvec_lines[2]+3].split("=")[1].replace(")", "").replace("(", "").split()).astype(float)
    except Exception as exc:
        print("Error trying to parse cellvec_lines[2] from file")
        print("Exception is:", exc)
        print("cellvec_lines:", cellvec_lines[2])
    cellvec.append(v1*celldim[0])
    cellvec.append(v2*celldim[0])
    cellvec.append(v3*celldim[0])
    cellparam = cellvec_2_cellparam(cellvec)
    return np.array(cellvec), celldim, cellparam
