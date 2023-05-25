def apply_coord_displacement(coord, atom, axis, displacement: float=0.01, units: str='angstrom'):
    mod_coord = coord.copy()
    mod_coord[atom][axis] += displacement
    return mod_coord

def findiff_displacements(coord, displacement: float=0.01, units: str='angstrom'):
    geoms = []
    names = []
    for idx, atom_coord in enumerate(coord):
        for jdx, axis in enumerate(atom_coord):
            ## Positive_Displacement
            geoms.append(apply_coord_displacement(coord, idx, jdx, displacement, units))
            names.append("_"+idx+"_"+jdx+"_pos")
            ## Negative_Displacement
            geoms.append(apply_coord_displacement(coord, idx, jdx, -displacement, units))
            names.append("_"+idx+"_"+jdx+"_neg")
    return geoms, names

#def read_forces(job: object, debug: int=0)
#    
#    ## Read the number of atoms from the gmol object
#    gmol   = job._recipe.subject
#    natoms = gmol.natoms
#
#    fmatrix = np.zeros((natoms,natoms))
#    for idx, comp in enumerate(job.computations.sort(key=lambda x: (x.index))):
#
#        # this could go to a different function that checks that all computations are good
#        if comp.isregistered and comp.isgood:
#        
#        # reads successive (in terms of index) computations
#        # and hessian is computed as the difference between forces


