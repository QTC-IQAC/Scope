import os
import numpy as np
import sys

scopepath = '/home/g4vela/SCOPE/Database_SCO/Scripts'
sys.path.append(scopepath)

import Parse_General
from Parse_General import search_string, read_lines_file

import Scope_Classes
from Scope_Classes import orbital_set, orbital, VNM 

import Constants

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

def G16_get_last_geom(lines):
    coord = []
    ldx, found1 = search_string("Coordinates (Angstroms)", lines, typ='last') 
    if found1:
        init_geom_line = ldx + 2
    ldx, found2 = search_string("Distance matrix", lines, typ='last') 
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

def G16_get_VNM(lines: list, witheigen: str=False):
    vnms = []
    freqs = []
    red_masses = []
    force_cnt = [] 
    IR_int = [] 
    atomidxs = [] 
    atnums = [] 
    xs = []
    ys = []
    zs = [] 

    ldx, found = search_string("Frequencies --", lines, typ='all')
    if found:
        index = 1
        for l in ldx:
            length = len(lines[l].split())
            f1 = float(lines[l].split()[2])
            freqs.append(f1*Constants.cm2har)
            r1 = float(lines[l+1].split()[3])
            red_masses.append(r1)
            fo1 = float(lines[l+2].split()[3])
            force_cnt.append(fo1)
            ir1 = float(lines[l+3].split()[3])
            IR_int.append(ir1)

            new_vnm = VNM(index, f1, r1, fo1, ir1)
            vnms.append(new_vnm)
            index += 1

            if length > 3:
                f2 = float(lines[l].split()[3])
                freqs.append(f2*Constants.cm2har)
                r2 = float(lines[l+1].split()[4])
                red_masses.append(r2)
                fo2 = float(lines[l+2].split()[4])
                force_cnt.append(fo2)
                ir2 = float(lines[l+3].split()[4])
                IR_int.append(ir2)

                new_vnm = VNM(index, f2, r2, fo2, ir2)
                vnms.append(new_vnm)
                index += 1

            if length > 4:
                f3 = float(lines[l].split()[4])
                freqs.append(f3*Constants.cm2har)
                r3 = float(lines[l+1].split()[5])
                red_masses.append(r3)
                fo3 = float(lines[l+2].split()[5])
                force_cnt.append(fo3)
                ir3 = float(lines[l+3].split()[5])
                IR_int.append(ir3)

                new_vnm = VNM(index, f3, r3, fo3, ir3)
                vnms.append(new_vnm)
                index += 1

    return vnms
