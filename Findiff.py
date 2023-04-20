def apply_coord_displacement(coord, atom, axis, displacement: float=0.01, units: str='angstrom'):
    mod_coord = coord.copy()
    mod_coord[atom][axis] += displacement
    return mod_coord

#def findiff_displacements(orig_coord, displacement: float=0.01, units: str='angstrom'):
#    geoms = []
#    names = []
#    for idx, atom_coord in enumerate(orig_coord):
#        for jdx, axis in enumerate(atom_coord):
#            ## Positive_Displacement
#            mod_coord = orig_coord.copy()
#            mod_coord[idx][jdx] = orig_coord[idx][jdx]+displacement
#            geoms.append(mod_coord)
#            names.append("_"+idx+"_"+jdx+"_pos")
#            ## Negative_Displacement
#            mod_coord = orig_coord.copy()
#            mod_coord[idx][jdx] = orig_coord[idx][jdx]-displacement
#            geoms.append(mod_coord)
#            names.append("_"+idx+"_"+jdx+"_neg")
#    return geoms, names

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
