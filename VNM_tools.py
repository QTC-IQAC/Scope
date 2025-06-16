import numpy as np


def write_vnm_dyn(initial_coord: list, vnm: object, amplitude: int=10, index: int=0, outfolder: str='./', labels: None=list):
    filename: str="dyn_vnm_"+str(index)+".xyz"
    initial_coord = np.array(initial_coord)
    with open(outfolder+filename, "w") as output_file:    
        for f in range(-amplitude,amplitude,1):
            natoms = len(vnm.xs)
            print(natoms, file=output_file)
            print("", file=output_file)
            for idx in range(natoms):
                vector = np.array([vnm.xs[idx], vnm.ys[idx], vnm.zs[idx]])
                coord  = initial_coord[idx] + vector*f*0.1
                if labels is None: label = elemdatabase.elementsym[vnm.atnums[idx]]
                else:              label = vnm.labels[idx]
                print(f"{label}  {coord[0]:.5f}  {coord[1]:.5f}  {coord[2]:.5f}", file=output_file)

def vnm_displacement(VNMs: list, initial_coord: list, which: list=[], which_side: str='positive', amplitude: int=6, debug: int=0):
    ### This function applies a displacement from the initial geometry 
    ### using either 'all' negative normal modes whose index is in 'which'
    ### Apparently, the normal modes in gaussian come in sqrt(amu)·bohr units
    ### So I had to do some changes

    ## Frequencies must have eigenvectors stored
    if not hasattr(VNMs[0],'xs'): return None
    if not VNMs[0].haseigenvec: print("VNMs do not have eigenvectors. Stopping"); return None 

    ## Store Initial Coord
    new_coord = initial_coord.copy()    
    if debug >= 1: print("initial coord:", new_coord[0])

    ## If empty which, then takes all
    if len(which) == 0: which = [vnm.index for vnm in VNMs]

    ## Applies Displacement
    for vnm in VNMs:
        if vnm.index in which:
            if debug >= 1: print("displacing VNM with frequency:", vnm.freq_cm)
            for idx in range(len(vnm.xs)):
                vector = np.array([vnm.xs[idx], vnm.ys[idx], vnm.zs[idx]])*np.sqrt(vnm.masses[idx])*Constants.bohr2angs
                #vector = np.array([vnm.xs[idx], vnm.ys[idx], vnm.zs[idx]])  as it was before

                ## The actual amount of displacement depends on the amplitude defined by the user
                ## ... and the freq_factor, which depends on the frequency
                if   vnm.freq < -20:  freq_factor = 0.1
                else:                 freq_factor = 0.2
                if debug >= 1: print(f"using factor={amplitude*freq_factor}")

                ## Apply displacement to coordinates
                if   which_side.lower() == 'positive': displacement = vector*amplitude*freq_factor
                elif which_side.lower() == 'negative': displacement = -vector*amplitude*freq_factor
                new_coord[idx] = new_coord[idx]+displacement
                if idx == 0 and debug >= 1: print("displaced coord:", new_coord[idx])
    return new_coord
