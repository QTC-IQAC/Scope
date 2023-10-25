import pickle
import numpy as np
import os
import shutil
from ast import literal_eval

#from Scope.Other import str_to_list

def save_binary(variable, pathfile, backup: bool=False):
    pathfile = pathfile.replace("lustre","home")
    if not backup:
        try: 
            file = open(pathfile,'wb')
            pickle.dump(variable,file)
            file.close()
        except Exception as exc:
            print("Error Saving Binary for pathfile:", pathfile)
            print(exc)
    if backup:
        print("Backup not implemented, file not saved")
#        if os.path.isfile(pathfile):
#            saved = False
#            idx = 1
#            while not saved and idx < 11:
#                new_pathfile = pathfile+"_backup"+str(idx)
#                if not os.path.isfile(new_pathfile): 
#                    shutil.copy(pathfile, new_pathfile) 
#                    file = open(pathfile,'wb')
#                    pickle.dump(variable,file)
#                    file.close()
#                    saved = True 

#######################
def load_binary(pathfile):
    with open(pathfile, "rb") as pickle_file:
        binary = pickle.load(pickle_file)
    return binary

#######################
def save_list_as_text(inplist: list, pathfile: str=os.getcwd()+"outfile.txt"):
    with open(pathfile, "w") as fil:
        for l in inplist:
            print(l, file=fil)

#######################
def writexyz(fdir, fname, labels, coord, charge: int=0, spin: int=1):
    if fdir[-1] != "/":
        fdir = fdir + "/"
    natoms = len(labels)
    fullname = fdir + fname
    with open(fullname, "w") as fil:
        print(natoms, file=fil)
        print(charge, spin, file=fil)
        for idx, l in enumerate(labels):
            print("%s  %.6f  %.6f  %.6f" % (l, coord[idx][0], coord[idx][1], coord[idx][2]),file=fil)

def print_xyz(labels, coord):
    for idx, l in enumerate(labels):
        print("%s  %.6f  %.6f  %.6f" % (l, coord[idx][0], coord[idx][1], coord[idx][2]))


#######################
def read_xyz(xyz_file):
    assert(xyz_file[-4:] == ".xyz")
    
    labels = []
    coord = []
    
    try:    xyz = open(xyz_file, "r")
    except: print("Could not read xyz file: {0}".format(xyz_file))

    n_atoms = xyz.readline()
    title = xyz.readline()
    for line in xyz:
        line_data = line.split()
        if len(line_data) == 4:
            label, x, y, z = line.split()
            coord.append([float(x), float(y), float(z)])
            labels.append(label)
        else:
            print("I can't read the xyz. It has =/ than 4 columns")

    xyz.close()
    return labels, coord

def read_user_input(message: str, rtext: bool=False, rtext_options: list=[], rtype: bool=False, rtype_options: list=[], limit_attempts: bool=True, attempts: int=3, debug: int=0):
    att = 0
    correct = False
    if limit_attempts:
        while att < attempts and not correct:
            opt = input(message)
            if debug > 0: print(f"Read opt={opt}, of type={type(opt)}, rtype={rtype}, rtext={rtext}")
            if rtext and not rtype:
                if opt in rtext_options:                                correct = True
                else:                                                   correct = False; att += 1
            elif rtype and not rtext:
                try:        opt = literal_eval(opt); isstr = False
                except:     isstr = True
                if debug > 0: print(f"isstr={isstr}, opt={opt}, type={type(opt)}")

                if type(opt) in rtype_options:                          correct = True
                else:                                                   correct = False; att += 1

            elif rtype and rtext:
                if type(opt) in rtype_options and opt in rtext_options: correct = True
                else:                                                   correct = False; att += 1
            else: correct = True
            if debug > 0: print(f"correct={correct}, #attempt={att}")
            
            if not correct: print(f"Please, try again. Options are: {rtext_options}")
        if correct: return opt
        else:       return None
