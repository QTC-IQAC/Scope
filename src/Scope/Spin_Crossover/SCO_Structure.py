import numpy as np
from Scope.Geometry import *

######
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
                targ = sorted(range(len(dist_FeX)), key=lambda k: dist_FeX[k])[idx]
                #targ = np.argpartition(dist_FeX, idx)[idx]
                dist_FeN.append(dist_FeX[targ])
                vec_FeN.append(vec_FeX[targ])        
                if debug >= 1: print(f"Target: {targ=} {dist_FeX[targ]=} {vec_FeX[targ]}") 
                
            if debug >= 1: print(f"dist_FeN=", dist_FeN)
            AvFeN = np.round(np.mean(dist_FeN),3)
            
            # Takes vec_FeN vectors and computes angles
            all_NFeN_angles = []
            for idx, v1 in enumerate(vec_FeN):
                for jdx, v2 in enumerate(vec_FeN):
                    if jdx > idx:
                        angle = get_angle(v1, v2)*180.0/3.141592 
                        all_NFeN_angles.append(angle)
            if debug >= 2: print(f"GEOM_SCO: {all_NFeN_angles=}")

            # Stores the 12th angles closest to 90 degrees
            all_NFeN_angles = np.array(all_NFeN_angles) - 90
            if debug >= 2: print(f"GEOM_SCO: {all_NFeN_angles=} - 90")

            NFeN_angles = [] 
            epsylon = 0
            for idx in range(0,12):
                targ = sorted(range(len(all_NFeN_angles)), key=lambda k: all_NFeN_angles[k])[idx]
                NFeN_angles.append(all_NFeN_angles[targ] + 90)
                epsylon += np.abs(all_NFeN_angles[targ])
                if debug >= 1: print(f"Target: {targ=} {all_NFeN_angles[targ]=} {epsylon=}") 
            epsylon = np.round(epsylon,3)
            if debug >= 2: print(f"GEOM_SCO: {NFeN_angles=}")

            AvNFeN_angle = np.round(np.mean(NFeN_angles),3)
            if debug >= 2: print(f"{NFeN_angles=}")
            
            alldist.append(AvFeN)
            allangle.append(AvNFeN_angle)
            allepsylon.append(epsylon)
    
    return alldist, allangle, allepsylon

######
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
