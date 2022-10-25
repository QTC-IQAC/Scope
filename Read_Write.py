import pickle
import numpy as np

def save_binary(variable, pathfile):
    file = open(pathfile,'wb')
    pickle.dump(variable,file)
    file.close()

def load_binary(pathfile):
    with open(pathfile, "rb") as pickle_file:
        binary = pickle.load(pickle_file)
    return binary

def writexyz(fdir, fname, labels, pos, charge: int=0, spin: int=1):
    if fdir[-1] != "/":
        fdir = fdir + "/"
    natoms = len(labels)
    fullname = fdir + fname
    with open(fullname, "w") as fil:
        print(natoms, file=fil)
        print(charge, spin, file=fil)
        for idx, l in enumerate(labels):
            print("%s  %.6f  %.6f  %.6f" % (l, pos[idx][0], pos[idx][1], pos[idx][2]),file=fil)
