#!/usr/bin/env python

import numpy as np
import os

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def getangle(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def geom_sco_from_xyz(labels, pos, debug=0):   
    #Computes Structural Variables for all Fe atoms of a given structure
    #Returns two lists, one for FeN and another for NFeN, with all values.
    # Epsylon is eq. 1 of Halcrow's Chem. Soc. Rev., 2011,40, 4119-4142
    # Angle is the average of cis angles, without the deviation from 90º
    
    alldist = []
    allangle = []
    allepsylon = []
    
    #Computes Fe-N distances and vectors
    NFe = labels.count("Fe")
    for idx, p1 in enumerate(pos):
        if labels[idx] == 'Fe':
            
            if debug >= 2: print(f"Atom {idx} is an Iron: {p1}")
            dist_FeX = []
            vec_FeX = []
            for jdx, p2 in enumerate(pos):
                if labels[jdx] == 'N':
                    if debug >= 2: print(f"Atom {jdx} is a Nitrogen: {p2}")
                    dist_FeX.append(np.linalg.norm(np.subtract(p1,p2)))
                    vec_FeX.append(np.subtract(p1,p2))
                    
            if debug >= 2: print("dist_FeX", dist_FeX)
            
            # Checks that it detected 6 N atoms
            if len(vec_FeX) < 6: 
                print("less than 6 N atoms in the Fe coord sphere. Found", len(vec_FeX))
                print("printing coordinates to debug:")
                for idx, at in enumerate(labels):
                    print("%s   %.6f   %.6f   %.6f" % (labels[idx], pos[idx][0], pos[idx][1], pos[idx][2]))
                 
            # Computes average FeN distance for the closest 6 N atoms, and prepares their vectors (vec_FeN)
            dist_FeN = []
            vec_FeN = []
            for idx in range(0,6):
                targ = np.argpartition(dist_FeX, idx)[idx]
                dist_FeN.append(dist_FeX[targ])
                vec_FeN.append(vec_FeX[targ])        
                
            if debug >= 1: print(f"dist_FeN=", dist_FeN)
            AvFeN = np.mean(dist_FeN)
            
            # Takes vec_FeN vectors and computes angles
            NFeN_angles = []
            for idx, v1 in enumerate(vec_FeN):
                for jdx, v2 in enumerate(vec_FeN):
                    if jdx > idx:
                        if debug >= 2: print("Computing angles between N atoms",idx,"and",jdx)
                        angle = getangle(v1, v2)*180.0/3.141592 
                        if angle < 120.0 and angle > 60.0: 
                            NFeN_angles.append(angle)
                        else: 
                            if debug >= 2: print("Angle between N atoms",idx,"and",jdx,"discarded")

            assert len(NFeN_angles) == 12
            epsylon = 0
            for a in NFeN_angles:
                epsylon += np.abs(90-a)

            AvNFeN_angle = np.mean(NFeN_angles)
            if debug >= 2: print(f"{NFeN_angles=}")
            
            alldist.append(AvFeN)
            allangle.append(AvNFeN_angle)
            allepsylon.append(epsylon)
    
    return alldist, allangle, allepsylon

def readxyz(file):
    labels = []
    pos = []
    xyz = open(file, "r")
    n_atoms = xyz.readline()
    title = xyz.readline()
    for line in xyz:
        line_data = line.split()
        if len(line_data) == 4:
            label, x, y, z = line.split()
            pos.append([float(x), float(y), float(z)])
            labels.append(label)
        else:
            print("I can't read the xyz. It has =/ than 4 columns")
    xyz.close()

    return labels, pos

def guess_spin_state(ox_state: str, dist: float, debug: int=0):
    if debug > 0: print(f"Received Ox_state={ox_state} and dist={dist}")
    if ox_state == 2 or ox_state == '2':
        if float(dist) >= 1.8 and float(dist) <= 2.04: guess_spin = 'LS'
        elif float(dist) > 2.04 and float(dist) <= 2.10: guess_spin = 'Unk'
        elif float(dist) > 2.10 and float(dist) <= 2.23: guess_spin = 'HS'
        elif float(dist) > 2.23 : guess_spin = 'Unk'
        else: guess_spin = 'Unk'
    else: guess_spin = 'Unk'
    if debug > 0: print(f"Guess_spin={guess_spin}") 
    return guess_spin
