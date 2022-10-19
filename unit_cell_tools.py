import numpy as np

def get_unit_cell_volume(a, b, c, alpha, beta, gamma):
    return float(a*b*c*np.sin(np.deg2rad(beta)))

def cellvec_2_cellparam(cellvec):
    a = np.linalg.norm(cellvec[0])
    b = np.linalg.norm(cellvec[1])
    c = np.linalg.norm(cellvec[2])
    
    tmp1 = cellvec[0][0]*cellvec[1][0]+cellvec[0][1]*cellvec[1][1]+cellvec[0][2]*cellvec[1][2]
    tmp2 = cellvec[0][0]*cellvec[2][0]+cellvec[0][1]*cellvec[2][1]+cellvec[0][2]*cellvec[2][2]
    tmp3 = cellvec[1][0]*cellvec[2][0]+cellvec[1][1]*cellvec[2][1]+cellvec[1][2]*cellvec[2][2]    
    tmp4 = a * b
    tmp5 = a * c
    tmp6 = b * c
    tmp7 = tmp1/tmp4
    tmp8 = tmp2/tmp5
    tmp9 = tmp3/tmp6
    
    alpha = np.arccos(tmp7)/(2*np.pi/360)
    beta = np.arccos(tmp8)/(2*np.pi/360)
    gamma = np.arccos(tmp9)/(2*np.pi/360)
    
    return list([a, b, c, alpha, beta, gamma])
