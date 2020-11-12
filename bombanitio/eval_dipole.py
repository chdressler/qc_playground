import sys
import numpy as np
from bombanitio import eval_over
#import eval_four.delta
#from bombanitio import eval_four.delta 
from bombanitio.eval_four import delta 

def get_dip(type1, type2, exp1, exp2, r1, r2, k):
    if type1 == 0 and type2 ==0:
        #dip = dip_ss(type1, type2,exp1, exp2, r1, r2)
        dip = dip_ss(exp1, exp2, r1, r2,k)
    elif type1 > 0  and type2 > 0:
        dip = dip_pp(type1, type2,exp1, exp2, r1, r2,k)
    
    elif type1 > 0  and type2 ==0:
        dip = dip_ps(type1, type2, exp1, exp2, r1, r2,k)
    elif type1 == 0 and type2 > 0:
        dip = dip_ps(type2, type1, exp2, exp1, r2, r1,k)
    else:
        print("error orbital combination does not exist!")
        sys.exit()
    return dip

def gab(exp1, exp2, r1, r2):
    return np.exp(- (exp1 * exp2)/(exp1 +exp2) * np.linalg.norm(r1-r2)**2)



#def over_ss(exp1, exp2, r1, r2):
#def over_pp(type1, type2, exp1, exp2, r1, r2):


def dip_ss(exp1, exp2, r1, r2, k):
    #dip = eval_over.over_ps(k+1 , 0, exp1, exp2, r1, r2, norm = False)  + r1[k] *eval_over.over_ss(exp1, exp2, r1, r2, norm = False)
    dip = eval_over.over_ps(k+1 , 0, exp1, exp2, r1, r2, norm = False)  + r1[k] *eval_over.over_ss(exp1, exp2, r1, r2, norm = False)
    #dip = eval_over.over_ps(k+1 , 0, exp1, exp2, r1, r2)  * r1[k] *eval_over.over_ss(exp1, exp2, r1, r2)
    dip *= (2.0 * exp1 / np.pi)**0.75 
    dip *= (2.0 * exp2 / np.pi)**0.75 
    #print(dip, exp1, exp2, r1, r2, k)
    return dip

def dip_pp(type1, type2, exp1, exp2, r1, r2, k):
    i = type1 - 1
    j = type2 - 1
    #dip = r1[k] * eval_over.over_pp(type1, k+1 , exp1, exp2, r1, r2, norm = False) +  1/(2*exp1) * eval_over.over_ps(0 , type2, exp1, exp2, r1, r2, norm = False) * delta(i,k)
    dip = r1[k] * eval_over.over_pp(type1, type2 , exp1, exp2, r1, r2, norm = False) +  1/(2*exp1) * eval_over.over_ps(0 , type2, exp1, exp2, r1, r2, norm = False) * delta(i,k)
    dip += - exp2 * np.pi**1.5/(exp1 + exp2)**3.5*gab(exp1, exp2, r1, r2)*(0.5*(r1[i] - r2[i])*delta(j,k) + 0.5 * (r1[j]-r2[j])* delta(i,k) + 0.5 * (r1[k]-r2[k])*delta(i,j) - exp1*exp2/(exp1+exp2)*(r1[i]- r2[i])*(r1[j]- r2[j])*(r1[k]- r2[k]))
    dip *= (128 *exp1**5 / np.pi**3)**0.25
    dip *= (128 *exp2**5 / np.pi**3)**0.25
    #print(dip, exp1, exp2, r1, r2, k)
    return dip

def dip_ps(type1, type2, exp1, exp2, r1, r2,k):
    i  =  type1 - 1
    dip = eval_over.over_pp(type1, k+1 , exp1, exp2, r1, r2, norm = False) + r2[k] * eval_over.over_ps(type1, type2, exp1, exp2, r1, r2, norm = False)
    dip *= (128 *exp1**5 / np.pi**3)**0.25
    dip *= (2.0 * exp2 / np.pi)**0.75
    #print(dip, exp1, exp2, r1, r2, k)
    return dip


def get_dip_mat(noo, orb_coord, orb_type, coef_mat, exp_mat, kk):
    dip_mat = np.zeros((noo, noo))
    for i in range(noo):
        for j in range(noo):
    #        print("drin",i,j, noo,)
            #if (i ==  5 or i ==6 or i == 0 or i ==1) and (j ==  5 or j ==6 or j == 0 or j ==1) and kk = 0:
            for k in range(3):
                for l in range(3):
                    #print(get_diplap(orb_type[i], orb_type[j], exp_mat[i,k], exp_mat[j,l], orb_coord[i,:], orb_coord[j,:]))
#                   
                    dip_mat[i,j] += coef_mat[i,k] *  coef_mat[j,l] * get_dip(orb_type[i], orb_type[j], exp_mat[i,k], exp_mat[j,l], orb_coord[i,:], orb_coord[j,:], kk) 
                    #print(dip_mat
                    #print
    return dip_mat


def  get_nuc_pole(coord, zoa):
    nuc_pole = np.zeros((3))
    noa = zoa.shape[0]
    for i in range(noa):
        nuc_pole += coord[i] * zoa[i]
    return nuc_pole


def get_el_dipole(p, dip_mat):
        dipole = dip_mat * p
        dipole = dipole.sum(axis=(1,2))
        return dipole

