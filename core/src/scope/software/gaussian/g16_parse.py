import scope.constants 
from   scope.parse_general import search_string, read_lines_file

##############
### STATUS ###
##############
def parse_status_finished(lines):
    search1, found1 = search_string("Normal termination", lines, typ="last")
    search2, found2 = search_string("Error termination", lines, typ="last")
    if found1 or found2: return True
    else:                return False

def parse_opt_status(lines):
    linenum, found = search_string("Optimization completed", lines, typ="last")
    return found

def parse_freq_status(lines):
    linenum, found = search_string("- Thermochemistry -", lines, typ="last")
    return found

def parse_force_status(lines):
    linenum, found = search_string("Forces (Hartrees/Bohr)", lines, typ='last')
    return found

def parse_scf_status(lines):
    linenum1, found1 = search_string("SCF Done", lines, typ="last")
    linenum2, found2 = search_string("Convergence criterion not met",      lines, typ="last")
    if found1 and not found2:             return True
    elif found2:                          
        if linenum1 > linenum2:           return True ## Case of Quadratic Convergence
        else:                             return False
    else:                                 return None

def parse_coord_status(lines):
    linenum, found = search_string("Coordinates (Angstroms)"   , lines, typ="last")
    return found

def parse_timelimit_status(lines):
    # No way to set time limit in gaussian jobs. So a time limit is never printed on output
    return False

###############
### PARSING ###
###############
def parse_start_scf(lines, debug=0):
    start_SCF, found  = search_string("Cycle   1", lines)
    return start_SCF, found

def parse_end_scf(lines, debug=0):
    end_SCF, found  = search_string("SCF Done", lines)
    return end_SCF, found

def parse_geometry_from_step(lines, debug: int=0):
    from scope.elementdata import ElementData
    elemdatabase = ElementData()
    labels = []
    coords = []
    ldx, found1 = search_string("Coordinates (Angstroms)", lines, typ='first')
    if found1:
        init_geom_line = ldx + 3
        ldx, found2 = search_string("--------------------------------------",lines,typ='first',lowlim=init_geom_line)
        ldx += 1
        if found2:
            end_geom_line = ldx - 1
    if found1 and found2:
        for l in lines[init_geom_line:end_geom_line]:
            if len(l.split()) == 6:
                tmp = []
                atnum = l.split()[1]
                x     = l.split()[3]
                y     = l.split()[4]
                z     = l.split()[5]
                tmp.append(float(x))
                tmp.append(float(y))
                tmp.append(float(z))
                coords.append(tmp)
                lab = elemdatabase.elementsym[int(atnum)]
                labels.append(lab)
    return labels, coords

def parse_forces_from_step(lines, debug: int=0):
    forces = []
    ldx, found1 = search_string("Forces (Hartrees/Bohr)", lines, typ='first')
    if found1:
        init_line = ldx + 3
        ldx, found2 = search_string("--------------------------------------",lines,typ='first',lowlim=init_line)
        ldx += 1
        if found2:
            end_line = ldx - 1
    if found1 and found2:
        for l in lines[init_line:end_line]:
            if len(l.split()) == 5:
                tmp = []
                x     = l.split()[2]
                y     = l.split()[3]
                z     = l.split()[4]
                tmp.append(float(x))
                tmp.append(float(y))
                tmp.append(float(z))
                forces.append(tmp)
        return forces
    else: 
        return None

####################
## PARSING ENERGY ##
####################
def parse_energy_from_step(lines, debug: int=0):
    ldx, found = search_string("SCF Done", lines, typ='first')        
    if found:
        ener = float(lines[ldx].split()[4])
        return ener
    else: 
        print("G16_PARSE: Energy not found in lines.")
        return None

def parse_energy(lines, typ: str='last', debug: int=0):
    ldx, found = search_string("SCF Done", lines, typ=typ)        
    if found:
        ener = float(lines[ldx].split()[4])
        return ener
    else: 
        print("G16_PARSE: Energy not found in lines.")
        return None

#########################
## PARSING FREE ENERGY ##
#########################
def parse_free_energy_from_step(lines, debug: int=0):
    ldx, found = search_string(" Sum of electronic and thermal Free Energies=", lines, typ='first')
    if found:
        free_energy = float(lines[ldx].split()[7])
        return free_energy
    else:
        print("G16_PARSE: Free energy not found in the output file.")
        return None

def parse_free_energy(lines, typ: str='last', debug: int=0):
    ldx, found = search_string(" Sum of electronic and thermal Free Energies=", lines, typ=typ)
    if found:
        free_energy = float(lines[ldx].split()[7])
        return free_energy
    else:
        print("G16_PARSE: Free energy not found in the output file.")
        return None

######################################
## PARSING VIBRATIONAL NORMAL MODES ##
######################################
def parse_hp_vnms_from_step(lines: list, witheigen: bool=False, debug: int=0):
    ## New implementation when using freq(hpmodes) in Gaussian16 input line
    from scope.classes_qc import VNM 
    vnms = []

    ldx, found1 = search_string("Frequencies ---", lines, typ='all')
    kdx, found2 = search_string("Harmonic frequencies", lines, typ='last')
    if not found1 or not found2: return None

    index = 1
    for idx, l in enumerate(ldx):
        if idx < len(ldx)-1: bottom = ldx[idx+1]
        else:                bottom = len(lines)
        length          = int(len(lines[l].split('---')[1].split()))   ## Number of freq blocks in this line
        if debug > 0: print(f"Processing block {idx}/{len(ldx)} with {length} frequencies starting at line {l} and ending at {bottom}")

        # Gets the number of atoms in the first block of frequencies
        if witheigen and idx == 0:
            inieigen, found1 = search_string("Coord Atom Element:", lines, typ='first', lowlim=l)
            if not found1:
                print("Error: Could not find eigenvector block in VNM parsing")
                return None
            len_block_eigen = (bottom-2 - inieigen) 
            natoms = int(len_block_eigen/3)
            if debug > 0: print(f"Number of atoms in eigenvector block: {natoms}")

        for kdx in range(length):
            index_from_file = int(lines[l-2].split()[kdx])
            assert index == index_from_file                

            sym  = str(lines[l-1].split()[kdx])
            freq = float(lines[l].split('---')[1].split()[kdx])
            rm   = float(lines[l+1].split('---')[1].split()[kdx])
            fo   = float(lines[l+2].split('---')[1].split()[kdx])
            ir   = float(lines[l+3].split('---')[1].split()[kdx])
            new_vnm = VNM(index, freq, rm, fo, ir, sym)

            if witheigen:
                if debug > 0 and kdx == 0: print(f"Parsing eigenvector block from {inieigen} to {bottom}")
                atom_idx = []
                atnum = []
                x = []
                y = []
                z = []
    
                run_line = 0
                for jdx in range(natoms):
                    atom_idx.append(int(lines[l+5+run_line].split()[1]))
                    atnum.append(int(lines[l+5+run_line].split()[2]))
                    for cdx in range(3):
                        if   cdx == 0: x.append(float(lines[l+5+run_line].split()[3+kdx])); run_line += 1
                        elif cdx == 1: y.append(float(lines[l+5+run_line].split()[3+kdx])); run_line += 1
                        elif cdx == 2: z.append(float(lines[l+5+run_line].split()[3+kdx])); run_line += 1
                new_vnm.set_mode(atom_idx, atnum, x, y, z)  
            vnms.append(new_vnm)
            index += 1
    return vnms
   
def parse_vnms_from_step(lines: list, witheigen: bool=False, debug: int=0):
    from scope.classes_qc import VNM 
    vnms = []

    ldx, found1 = search_string("Frequencies --", lines, typ='all')
    kdx, found2 = search_string("X      Y      Z", lines, typ='all')
    assert len(ldx) == len(kdx)
    if not found1 or not found2: return None

    index = 1
    for idx, l in enumerate(ldx):
        
        if witheigen and idx == 0:
            width = ldx[1] - kdx[0] - 3  # Should be the number of atoms, useful when retrieving the eigenvectors
        
        length = len(lines[l].split())
        
        f1 = float(lines[l].split()[2])            
        index_from_file = int(lines[l-2].split()[0])
        assert index == index_from_file            
        sym1 = str(lines[l-1].split()[0])
        
        rm_l, found = search_string("Red. masses", lines, typ='first', lowlim=l, uplim=kdx[idx])
        if found: r1 = float(lines[rm_l].split()[3])
        else: r1 = float(0)
            
        fc_l, found = search_string("Frc consts", lines, typ='first', lowlim=l, uplim=kdx[idx])
        if found: fo1 = float(lines[fc_l].split()[3])
        else: fo1 = float(0)
            
        ir_l, found = search_string("IR Inten", lines, typ='first', lowlim=l, uplim=kdx[idx])
        if found: ir1 = float(lines[ir_l].split()[3])
        else: ir1 = float(0)

        new_vnm = VNM(index, f1, r1, fo1, ir1, sym1)
        
        if witheigen:
            atom_idx = []
            atnum = []
            x = []
            y = []
            z = []
            for l2 in range(kdx[idx]+1,kdx[idx]+width+1):
                line = lines[l2]
                atom_idx.append(int(line.split()[0]))
                atnum.append(int(line.split()[1]))
                x.append(float(line.split()[2]))
                y.append(float(line.split()[3]))
                z.append(float(line.split()[4]))
            new_vnm.set_mode(atom_idx, atnum, x, y, z)       

        vnms.append(new_vnm)
        index += 1

        if length > 3:
            f2 = float(lines[l].split()[3])              
            index_from_file = int(lines[l-2].split()[1])
            assert index == index_from_file                
            sym2 = str(lines[l-1].split()[1])
            
            rm_l, found = search_string("Red. masses", lines, typ='first', lowlim=l, uplim=kdx[idx])
            if found: r2 = float(lines[rm_l].split()[4])
            else: r2 = float(0)

            fc_l, found = search_string("Frc consts", lines, typ='first', lowlim=l, uplim=kdx[idx])
            if found: fo2 = float(lines[fc_l].split()[4])
            else: fo2 = float(0)

            ir_l, found = search_string("IR Inten", lines, typ='first', lowlim=l, uplim=kdx[idx])
            if found: ir2 = float(lines[ir_l].split()[4])
            else: ir2 = float(0)
            
            new_vnm = VNM(index, f2, r2, fo2, ir2, sym2)
            
            if witheigen:
                atom_idx = []
                atnum = []
                x = []
                y = []
                z = []
                for l2 in range(kdx[idx]+1,kdx[idx]+width+1):   
                    line = lines[l2]
                    atom_idx.append(int(line.split()[0]))
                    atnum.append(int(line.split()[1]))
                    x.append(float(line.split()[5]))
                    y.append(float(line.split()[6]))
                    z.append(float(line.split()[7]))
                new_vnm.set_mode(atom_idx, atnum, x, y, z)  
            
            vnms.append(new_vnm)
            index += 1
                        
        if length > 4:
            f3 = float(lines[l].split()[4])                
            index_from_file = int(lines[l-2].split()[2])
            assert index == index_from_file                
            sym3 = str(lines[l-1].split()[2])
            
            rm_l, found = search_string("Red. masses", lines, typ='first', lowlim=l, uplim=kdx[idx])
            if found: r3 = float(lines[rm_l].split()[5])
            else: r3 = float(0)

            fc_l, found = search_string("Frc consts", lines, typ='first', lowlim=l, uplim=kdx[idx])
            if found: fo3 = float(lines[fc_l].split()[5])
            else: fo3 = float(0)

            ir_l, found = search_string("IR Inten", lines, typ='first', lowlim=l, uplim=kdx[idx])
            if found: ir3 = float(lines[ir_l].split()[5])
            else: ir3 = float(0)
            
            new_vnm = VNM(index, f3, r3, fo3, ir3, sym3)
            
            if witheigen:
                atom_idx = []
                atnum = []
                x = []
                y = []
                z = []
                for l2 in range(kdx[idx]+1,kdx[idx]+width+1): 
                    line = lines[l2]
                    atom_idx.append(int(line.split()[0]))
                    atnum.append(int(line.split()[1]))
                    x.append(float(line.split()[8]))
                    y.append(float(line.split()[9]))
                    z.append(float(line.split()[10]))
                new_vnm.set_mode(atom_idx, atnum, x, y, z)  
            
            vnms.append(new_vnm)
            index += 1
    return vnms

################
### TD / TDA ###
################
def parse_exc_states(lines, nstates: int=10, debug: int=0):
    from scope.classes_qc import ExcitedState
    if not parse_status_finished(lines): 
        if debug > 0: print("PARSE_TD: Job not finished, excited states may be incomplete or inaccurate.")
        return None

    if debug > 0: print(f"PARSE_TD: Parsing excited states, looking for {nstates} states in the output lines.")
    es_list = []
    for st_num in range(1,nstates+1):

        # Format changes depending on the state number
        if st_num < 10:    
            line_nums, found = search_string(f"Excited State   {st_num}:", lines, typ='first')
        elif st_num >= 10 and st_num < 100: 
            line_nums, found = search_string(f"Excited State  {st_num}:", lines, typ='first')
        elif st_num >= 100 and st_num < 1000: 
            line_nums, found = search_string(f"Excited State {st_num}:", lines, typ='first')
        if not found: 
            print(f"PARSE_TD: Excited State {st_num} not found in lines") 
            return None

        # Parses relevant data
        _, _, idx, _, energy, _, wavelength, _, fosc, s2 = lines[line_nums].split() 
        idx         = int(idx.replace(':',''))
        s2          = float(s2.replace('<S**2>=',''))
        fosc        = float(fosc.replace('f=',''))
        wavelength  = float(wavelength)
        energy      = float(energy)
        if debug > 0: print(f"PARSE_TD: Parsed data for state {st_num}: energy={energy} eV, wavelength={wavelength} nm, fosc={fosc}, s2={s2}")

        if not st_num == idx:  
            raise ValueError(f"PARSE_TD: State number mismatch in TD parsing: expected {st_num}, found {idx} in line: {lines[line_nums]}")
        new_es = ExcitedState(st_num, energy, wavelength, fosc, s2)
        es_list.append(new_es)

    return es_list

############
### TIME ###
############
def G16_time_to_sec(time_list: list):
    floats = []
    strings = []
    for idx, entry in enumerate(time_list):
        try: floats.append(float(entry))
        except: strings.append(str(entry))
    total = 0
    for z in zip(floats, strings):
        if z[1] == 'days': multiplier = 86400
        elif z[1] == 'hours': multiplier = 3600
        elif z[1] == 'minutes': multiplier = 60
        elif z[1] == 'seconds' or z[1] == 'seconds.': multiplier = 1
        else: 
            print("Could not understand label", z[1]) 
            multiplier = 0
        total += z[0]*multiplier
    return total
