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
