from scope import constants
from copy import deepcopy
import numpy as np

def apply_coord_displacement(coord, atom: int, axis: int, displacement: float=0.01, units: str='angstrom'):
    mod_coord = deepcopy(coord)
    if units == 'angstrom':    mod_coord[atom][axis] += displacement
    elif units == 'bohr':      mod_coord[atom][axis] += displacement*Constants.bohr2angs
    else: print(f"FINDIFF: APPLY COORD DISPLACEMENT: could not understand the specified {units=}"); return None
    return mod_coord

def get_central_difference(f1, f2, displacement: float=0.01, units: str='angstrom', debug: int=0):
    ## displacement is how much the coordinate is displaced from the center point, so the difference between the positive and negative displacements is twice
    ## Units of the displacement, not units of the forces
    ## f1 must be the forces of the positive displacement
    ## f2 must be the forces of the negative displacement
    assert np.shape(f1) == np.shape(f2)    
    if units == 'angstrom': displacement = displacement*Constants.angs2bohr
    elif units == 'bohr':   pass
    else: print("FINDIFF: GET CENTRAL: could not understand the specified units"); return None
    cdiff = - ( f1 - f2 ) / (2.0 * displacement)
    return cdiff

def findiff_displacements(coord, displacement: float=0.01, units: str='bohr'):
    ## Default should be bohr instead
    if   units == "bohr": pass
    elif units == "angstrom": displacement = displacement * Constants.angs2bohr
    geoms = []
    names = []
    for idx, atom_coord in enumerate(coord):
        for jdx, axis in enumerate(atom_coord):
            ## Positive_Displacement
            geoms.append(apply_coord_displacement(coord, idx, jdx, displacement, units))
            names.append("_"+str(idx)+"_"+str(jdx)+"_pos")
            ## Negative_Displacement
            geoms.append(apply_coord_displacement(coord, idx, jdx, -displacement, units))
            names.append("_"+str(idx)+"_"+str(jdx)+"_neg")
    return geoms, names

def mass_weight_hessian(hessian, masses, debug: int=0):
    assert np.shape(hessian)[0] == 3*len(masses)
    natoms = len(masses)
    dim = np.shape(hessian)[0]
    mw_hessian = hessian.copy()
    ij = 0
    for idx in range(1,natoms+1):
        for jdx in range(1,3+1):
            ij += 1
            kl = 0
            for kdx in range(1,natoms+1):
                for ldx in range(1,3+1):
                    kl += 1
                    factor = 1.0/np.sqrt(masses[idx-1]*masses[kdx-1])
                    mw_hessian[ij-1,kl-1] = hessian[ij-1,kl-1]*factor
    return mw_hessian

def project_out(hessian, proj_rot: bool=False, proj_tra: bool=True, debug: int=0):
    from scope.other import gram_schmidt

    if not proj_rot and not proj_tra: return hessian
    
    dim = np.shape(hessian)[0]    
    um = np.identity(dim)
    
    if proj_rot: h = np.zeros((dim,6))
    else:        h = np.zeros((dim,3))
    
    ## translations ##
    for k in range(3):
        a = np.zeros((dim))
        for i in range(k, dim, 3):
            a[i] = 1./np.sqrt(dim)
        h[:,k] = a[:]
        
    ## Rotations: from 559 to 676 of vibanalysis.f 
    if proj_rot: pass # when implemented, this would create the 3 other columns of h. That is: h[:,3:6]
    
    # orthogonalize with gram_schmidt
    new_h = gram_schmidt(h)
    
    if not proj_rot: 
        ## we add the remaining 3 columns to h
        b = np.zeros((dim,3))
        final_h = np.concatenate((new_h, b), axis=1)
    else:
        final_h = new_h.copy()
                
    # projection matrix
    for j in range(dim):
        for i in range(dim):
            for k in range(6):
                um[i,j] = um[i,j] - final_h[i,k]*final_h[j,k]
    
    # correct hessian
    f = np.zeros((dim,dim))
    for j in range(dim):
        for i in range(dim):
            for k in range(dim):
                f[i,j] = f[i,j] + um[k,i]*hessian[k,j]
    
    proj_hessian = np.zeros((dim,dim))
    for j in range(dim):
        for i in range(dim):
            for k in range(dim):
                proj_hessian[i,j] = proj_hessian[i,j] + f[i,k]*um[k,j]
                
    return proj_hessian

def get_VNM_from_findiff(job: object, proj_rot: bool=False, proj_tra: bool=True, debug: int=0):
    from scope.operations.vecs_and_mats import symmetrize

    factor_cminv = 5140.487   ## must figure out where does it come from

    import numpy as np
    from scope.classes_qc import VNM

    ## Read the number of atoms from the source
    source = job._workflow.source
    natoms  = source.natoms
    masses = [elemdatabase.elementweight[l] for l in labels]     # For the mass-weighted hessian
    atnums  = [elemdatabase.elementnr[l] for l in labels]        # For the creation of VNM
    atomidxs = [i+1 for i in range(natoms)]                      # For the creation of VNM

    hessian = np.zeros((natoms*3,natoms*3))

    ## We collect all computations in a simpler variable
    comps = job.computations.sort(key=lambda x: (x.index))

    ## Checks that all computations are good and are registered, otherwise quits
    for idx, c in enumerate(comps):
        if c.isregistered and c.isgood: pass
        else:
            print(f"get_VNM_findiff: WARNING. Computation {idx} is reg={c.isregistered} good={c.isgood}")
            return None

    # Collects data
    for idx in range(3*natoms):
        # reads successive computations
        f1 = np.array(comps[2*idx].forces)
        f2 = np.array(comps[2*idx+1].forces)

        # and the dentral difference is computed as the difference between forces
        cdiff = get_central_difference(f1, f2)  # Nx3 components

        # insert cdiff into hessian
        for jdx, atom in enumerate (cdiff):
            for kdx, axis in enumerate(atom):
                row = 3*jdx + kdx
                column = idx
                hessian[column, row] = axis

    ## Projects rotations and/or translations, symmetrizes, and mass-weights
    if proj_rot or proj_tra:
        proj_hessian = project_out(hessian, proj_rot=proj_rot, proj_tra=proj_tra, debug=debug)
    else:
        proj_hessian = hessian.copy()

    sym_hessian = symmetrize(hessian)
    mw_hessian = mass_weight_hessian(sym_hessian, masses)

    # Diagonalize the mass-weighted hessian. Eigenval should have mass units in it
    eigenval, eigenvec = np.linalg.eig(mw_hessian)

    nvnm   = np.shape(eigenval)[0]
    assert nvnm % 3 == 0 
    assert natoms == int(nvnm / 3)

    ## Transform the eigenvalues to cm-1
    e_cm = []
    for e in enumerate(eigenval):
        if e >= 0: sign =  1
        else:      sign = -1
        e_cm.append(float(np.sqrt(abs(e))*factor_cminv)*sign)

    ## Transpose the matrix of eigenvectors, since np.linalg.eig returns it transposed
    eigenvec = eigenvec.T

    VNMs = []
    for idx, e in enumerate(e_cm):
        new_VNM = VNM(idx, e)
        evec = np.reshape(eigenvec[i],[nvnm,3])
        xs = evec[:,0]
        ys = evec[:,1]
        zs = evec[:,2]
        new_VNM.eigenvec(atomidxs,atnums,xs,ys,zs)
        VNMs.append(new_VNM)

    # So far, VNMs have been created in reverse order (from largest to smallest frequency). Here we reverse
    VNMs.reverse()

    return VNMs

