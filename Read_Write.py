import pickle
import numpy as np
import os
import shutil

def save_binary(variable, pathfile, backup: bool=False):
    if not backup:
        file = open(pathfile,'wb')
        pickle.dump(variable,file)
        file.close()
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
