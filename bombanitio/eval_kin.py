import sys
import numpy as np

def get_kin(type1, type2, exp1, exp2, r1, r2):
    if type1 == 0 and type2 ==0:
        #kin = kin_ss(type1, type2,exp1, exp2, r1, r2)
        kin = kin_ss(exp1, exp2, r1, r2)
    elif type1 > 0  and type2 > 0:
        kin = kin_pp(type1, type2,exp1, exp2, r1, r2)
    
    elif type1 > 0  and type2 ==0:
        kin = kin_ps(type1, type2, exp1, exp2, r1, r2)
    elif type1 == 0 and type2 > 0:
        kin = kin_ps(type2, type1, exp2, exp1, r2, r1)
    else:
        print("error orbital combination does not exist!")
        sys.exit()
    return kin

def gab(exp1, exp2, r1, r2):
    return np.exp(- (exp1 * exp2)/(exp1 +exp2) * np.linalg.norm(r1-r2)**2)

def kin_ss(exp1, exp2, r1, r2):
    kin =  gab(exp1, exp2, r1, r2) * exp1 * exp2 * np.pi**1.5 / (exp1 + exp2)**2.5 * (3 - 2 * exp1 * exp2 / ( exp1 + exp2) * np.linalg.norm( r1 - r2)**2)
    kin *= (2.0 * exp1 / np.pi)**0.75 
    kin *= (2.0 * exp2 / np.pi)**0.75 
    return kin

def kin_pp(type1, type2, exp1, exp2, r1, r2):
    i = type1 - 1
    j = type2 - 1
    if i == j:
        #kin = gab(exp1, exp2, r1, r2) * np.pi**1.5 / (exp1 + exp2)**2.5 * ( 0.5 - exp1 * exp2 / (exp1 + exp2) (r1[i] - r2[i]) * (r1[j] - r2[j])
        #kin =   0.5 - exp1 * exp2 / (exp1 + exp2) * (r1[i] - r2[i]) * (r1[j] - r2[j])
        kin =  2.5 - exp1 * exp2 / ( exp1 + exp2) * np.linalg.norm(r1-r2)**2 
        kin += exp1 * exp2 / (exp1 + exp2) * ( 2 * exp1 * exp2 / (exp1 + exp2) * np.linalg.norm( r1 - r2)**2  - 7 )* (r1[i] - r2[i]) * (r1[j] - r2[j])
    else:    
        #kin = - exp1 * exp2 / (exp1 + exp2) * (r1[i] - r2[i]) * (r1[j] - r2[j])
        kin = exp1 * exp2 / (exp1 + exp2) * ( 2 * exp1 * exp2 / (exp1 + exp2) * np.linalg.norm( r1 - r2)**2  - 7 )* (r1[i] - r2[i]) * (r1[j] - r2[j])
    #kin *= gab(exp1, exp2, r1, r2) * np.pi**1.5 / (exp1 + exp2)**2.5 
    kin *= gab(exp1, exp2, r1, r2) *   exp1 * exp2  * np.pi**1.5 / (exp1 + exp2)**3.5 
    kin *= (128 *exp1**5 / np.pi**3)**0.25
    kin *= (128 *exp2**5 / np.pi**3)**0.25
    return kin

def kin_ps(type1, type2, exp1, exp2, r1, r2):
    i  =  type1 - 1
    kin = gab(exp1, exp2, r1, r2) * exp1 * exp2**2  * np.pi**1.5 / (exp1 + exp2)**3.5 * (2 * exp1 * exp2 / ( exp1 + exp2) * np.linalg.norm(r1-r2)**2 -5)  * ( r1[i] - r2[i])
    kin *= (128 *exp1**5 / np.pi**3)**0.25
    kin *= (2.0 * exp2 / np.pi)**0.75 
    return kin


def get_kin_mat(noo, orb_coord, orb_type, coef_mat, exp_mat):
    kin_mat = np.zeros((noo, noo))
    for i in range(noo):
        for j in range(noo):
            for k in range(3):
                for l in range(3):
                    #print(get_kinlap(orb_type[i], orb_type[j], exp_mat[i,k], exp_mat[j,l], orb_coord[i,:], orb_coord[j,:]))
#
                    kin_mat[i,j] += coef_mat[i,k] *  coef_mat[j,l] * get_kin(orb_type[i], orb_type[j], exp_mat[i,k], exp_mat[j,l], orb_coord[i,:], orb_coord[j,:]) 
    return kin_mat
