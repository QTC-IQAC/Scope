import numpy as np

def write_vnm_dyn(vnm: object, maxfactor: int=10, index: int=0, outfolder: str='./'):
    filename: str="dyn_vnm_"+str(index)+".xyz"
    with open(outfolder+filename, "w") as output_file:    
        for f in range(-maxfactor,maxfactor,1):
            natoms = len(vnm.xs)
            print(natoms, file=output_file)
            print("", file=output_file)
            for idx in range(natoms):
                vector = np.array([vnm.xs[idx], vnm.ys[idx], vnm.zs[idx]])
                coord= geom[idx]+vector*f*0.1
                label = elemdatabase.elementsym[vnm.atnums[idx]]
                print(f"{label}  {coord[0]:.5f}  {coord[1]:.5f}  {coord[2]:.5f}", file=output_file)

def vnm_displacement(VNMs: list, initial_coord: list, which: list=[], which_side: str='positive', maxfactor: int=10, debug: int=0):
    # Applies a displacement from the initial geometry using either 'all negative normal modes whose index is in 'which'
    if len(which) == 0: which = [vnm.index for vnm in VNMs]
    new_coord = initial_coord.copy()    
    if debug >= 1: print("initial coord:", new_coord[0])
    for vnm in VNMs:
        if vnm.freq < 0.0 and vnm.index in which:
            if debug >= 1: print("displacing VNM:", vnm.freq_cm)
            for idx in range(len(vnm.xs)):
                vector = np.array([vnm.xs[idx], vnm.ys[idx], vnm.zs[idx]])
                if which_side.lower() == 'positive': displacement = vector*maxfactor*0.1
                elif which_side.lower() == 'negative': displacement = -vector*maxfactor*0.1
                new_coord[idx] = new_coord[idx]+displacement
                if idx == 0 and debug >= 1: print("displaced coord:", new_coord[idx])
    return new_coord
