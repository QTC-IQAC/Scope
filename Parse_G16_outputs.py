import os
import numpy as np
import sys

from Test_V3.Parse_General import search_string, read_lines_file
from Test_V3.Classes_QC import orbital_set, orbital, VNM 
from Test_V3 import Constants

def search_G16_last_geom_lines(filepath: str):
    warning = False
    try:
        lines = read_lines_file(filepath)
        for idx, l in enumerate(lines):
            lines[idx] = l.strip('\n')

        good_key = False
        key_str = ["Normal termination", "SCF Done", "Population analysis using the SCF Density"]
        key_line = []
        key_found = []
        for sx, string in enumerate(key_str):
            search, found = search_string(string, lines, typ="last")
            key_line.append(search)
            key_found.append(found)
            if not found: warning = True
    except Exception as exc:
        print("Warning searching last geometry in file:", filepath)
        warning = True
    if not warning: return lines[key_line[1]:key_line[0]], warning
    else: return [], True 

def read_G16_orbitals(filepath: str, geomlines):
    orb_eigen_strings = ["Alpha  occ. eigenvalues --", "Alpha virt. eigenvalues --",  "Beta  occ. eigenvalues --", " Beta virt. eigenvalues --"]
    alpha_occupied = orbital_set("Alpha Occupied")
    alpha_virtual = orbital_set("Alpha Virtual")
    beta_occupied = orbital_set("Beta Occupied")
    beta_virtual = orbital_set("Beta Virtual")
    
    all_orbital_sets = []
    all_orbital_sets.append(alpha_occupied)
    all_orbital_sets.append(alpha_virtual)
    all_orbital_sets.append(beta_occupied)
    all_orbital_sets.append(beta_virtual)    

    for idx, orbset in enumerate(all_orbital_sets):
        orbset.get_eigenvalues(orb_eigen_strings[idx], geomlines)
        for line in orbset.eigen_search_lines:
            en = line.split("--")[1].split()
            for e in en:
                if "Occupied" in orbset.name: occ = int(2)
                elif "Virtual" in orbset.name: occ = int(0)
                else: occ = int(10)
                orbset.eigenvalues.append(e)
                orbset.add_orbital(orbital(len(orbset.eigenvalues), occ, e))
    return all_orbital_sets 
    
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

def G16_get_last_geom(lines, debug: int=0):
    coord = []
    ldx, found1 = search_string("Coordinates (Angstroms)", lines, typ='last')
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
                x = l.split()[3]
                y = l.split()[4]
                z = l.split()[5]
                tmp.append(float(x))
                tmp.append(float(y))
                tmp.append(float(z))
                coord.append(tmp)
    return coord

def G16_get_last_energy(lines, debug: int=0):
    ldx, found = search_string("SCF Done", lines, typ='last')        
    if found: ener = float(lines[ldx].split()[4])
    else: ener = float(0.0)
    return ener

def G16_get_freqs(lines):
    freqs = []
    ldx, found = search_string("Frequencies --", lines, typ='all')
    if found:
        for l in ldx:
            length = len(lines[l].split())
            f1 = float(lines[l].split()[2])
            freqs.append(f1)
            if length > 3:
                f2 = float(lines[l].split()[3])
                freqs.append(f2)
            if length > 4:
                f3 = float(lines[l].split()[4])
                freqs.append(f3)
    return freqs

def G16_get_VNM(lines: list, witheigen: bool=False, debug: int=0):
    vnms = []

    ldx, found1 = search_string("Frequencies --", lines, typ='all')
    kdx, found2 = search_string("X      Y      Z", lines, typ='all')
    assert len(ldx) == len(kdx)
    if found1 and found2:
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
                new_vnm.eigenvec(atom_idx, atnum, x, y, z)       

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
                    new_vnm.eigenvec(atom_idx, atnum, x, y, z)  
                
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
                        x.append(float(line.split()[5]))
                        y.append(float(line.split()[6]))
                        z.append(float(line.split()[7]))
                    new_vnm.eigenvec(atom_idx, atnum, x, y, z)  
                
                vnms.append(new_vnm)
                index += 1
    return vnms
