import sys
import numpy as np

def get_overlap(type1, type2, exp1, exp2, r1, r2):
    if type1 == 0 and type2 ==0:
        #over = over_ss(type1, type2,exp1, exp2, r1, r2)
        over = over_ss(exp1, exp2, r1, r2)
    elif type1 > 0  and type2 > 0:
        over = over_pp(type1, type2,exp1, exp2, r1, r2)
    
    elif type1 > 0  and type2 ==0:
        over = over_ps(type1, type2, exp1, exp2, r1, r2)
    elif type1 == 0 and type2 > 0:
        over = over_ps(type2, type1, exp2, exp1, r2, r1)
    else:
        print("error orbital combination does not exist!")
        sys.exit()
    return over

def gab(exp1, exp2, r1, r2):
    return np.exp(- (exp1 * exp2)/(exp1 +exp2) * np.linalg.norm(r1-r2)**2)

def over_ss(exp1, exp2, r1, r2, norm =  True):
    over = (np.pi/(exp1 + exp2))**1.5 * gab(exp1, exp2, r1, r2)
    if norm:
        over *= (2.0 * exp1 / np.pi)**0.75 
        over *= (2.0 * exp2 / np.pi)**0.75 
    return over

def over_pp(type1, type2, exp1, exp2, r1, r2 , norm =  True):
    i = type1 - 1
    j = type2 - 1
    if i == j:
        #over = gab(exp1, exp2, r1, r2) * np.pi**1.5 / (exp1 + exp2)**2.5 * ( 0.5 - exp1 * exp2 / (exp1 + exp2) (r1[i] - r2[i]) * (r1[j] - r2[j])
        over =   0.5 - exp1 * exp2 / (exp1 + exp2) * (r1[i] - r2[i]) * (r1[j] - r2[j])
    else:    
        over = - exp1 * exp2 / (exp1 + exp2) * (r1[i] - r2[i]) * (r1[j] - r2[j])
    over *= gab(exp1, exp2, r1, r2) * np.pi**1.5 / (exp1 + exp2)**2.5 
    if norm == True:
        over *= (128 *exp1**5 / np.pi**3)**0.25
        over *= (128 *exp2**5 / np.pi**3)**0.25
    return over

def over_ps(type1, type2, exp1, exp2, r1, r2, norm =  True):
    i  =  type1 - 1
    over = -gab(exp1, exp2, r1, r2) * exp2 * np.pi**1.5 / (exp1 + exp2)**2.5  * ( r1[i] - r2[i])
    if norm == True:
        over *= (128 *exp1**5 / np.pi**3)**0.25
        over *= (2.0 * exp2 / np.pi)**0.75 
    return over


def get_overlap_mat(noo, orb_coord, orb_type, coef_mat, exp_mat):
    over_mat = np.zeros((noo, noo))
    for i in range(noo):
        for j in range(noo):
            for k in range(3):
                for l in range(3):
                    #print(get_overlap(orb_type[i], orb_type[j], exp_mat[i,k], exp_mat[j,l], orb_coord[i,:], orb_coord[j,:]))
#
                    over_mat[i,j] += coef_mat[i,k] *  coef_mat[j,l] * get_overlap(orb_type[i], orb_type[j], exp_mat[i,k], exp_mat[j,l], orb_coord[i,:], orb_coord[j,:]) 
    return over_mat
