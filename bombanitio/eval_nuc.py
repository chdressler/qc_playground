import sys
import numpy as np
import scipy.special
from bombanitio.special import specialx

def get_nuc(type1, type2, exp1, exp2, r1, r2, r3):
    if type1 == 0 and type2 ==0:
        #nuc = nuc_ss(type1, type2,exp1, exp2, r1, r2)
        nuc = nuc_ss(exp1, exp2, r1, r2, r3)
    elif type1 > 0  and type2 > 0:
        nuc = nuc_pp(type1, type2,exp1, exp2, r1, r2, r3)
    
    elif type1 > 0  and type2 ==0:
        nuc = nuc_ps(type1, type2, exp1, exp2, r1, r2, r3)
    elif type1 == 0 and type2 > 0:
        nuc = nuc_ps(type2, type1, exp2, exp1, r2, r1, r3)
    else:
        print("error orbital combination does not exist!")
        sys.exit()
    return nuc


#specialx.specialfunc(3,1)
#def nu(exp1, exp2):
#    return exp1 +  exp2
#def t
#
#nu = exp1 +  exp2
#rr = exp1 * r1 + exp2 * r2 / ( exp1 + exp2) - r3
#t = nu**0.5 * np.linalg.norm(rr)

#np.longdouble


def s1(t):
    return specialx.specialfunc(1,t)[0]
    #t = np.float128(t)
    #return scipy.special.erf(t)/t

def s2(t): 
    return specialx.specialfunc(2,t)[1]
    #t = np.float128(t)
    #return (2 * np.pi**-0.5 * t *  np.exp(-t**2)- scipy.special.erf(t) )/ t**3

def s3(t):
    return specialx.specialfunc(3,t)[2]
    #t = np.float128(t)
    #return (3* scipy.special.erf(t) - 2*np.pi**1.5*(3*t +  2* t**3)*np.exp(-t**2))/(t**5)

def gab(exp1, exp2, r1, r2):
    return np.exp(- (exp1 * exp2)/(exp1 +exp2) * np.linalg.norm(r1-r2)**2)

def nuc_ss(exp1, exp2, r1, r2, r3):
    nu = exp1 +  exp2
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - r3
    t = nu**0.5 * np.linalg.norm(rr)
    nuc = - gab(exp1, exp2, r1, r2)* np.pi**1.5 / nu *  s1(t)
    nuc *= (2.0 * exp1 / np.pi)**0.75 
    nuc *= (2.0 * exp2 / np.pi)**0.75 
    #print(exp1, exp2, r1,r2, r3, nuc, nu,rr, t)
    return nuc

def nuc_pp(type1, type2, exp1, exp2, r1, r2, r3):
    nu = exp1 +  exp2
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - r3
    t = nu**0.5 * np.linalg.norm(rr)
    i = type1 - 1
    j = type2 - 1
    if i == j:
        #nuc = gab(exp1, exp2, r1, r2) * np.pi**1.5 / (exp1 + exp2)**2.5 * ( 0.5 - exp1 * exp2 / (exp1 + exp2) (r1[i] - r2[i]) * (r1[j] - r2[j])
        #nuc =   0.5 - exp1 * exp2 / (exp1 + exp2) * (r1[i] - r2[i]) * (r1[j] - r2[j])
        nuc = nu * rr[i]*rr[j] * s3(t) + (1 + 2* exp1 *(r1[j]-r2[j])*rr[i] - 2*exp2 * (r1[i] - r2[i])*rr[j]) * s2(t) + (2 - (4*exp1*exp2*(r1[i] - r2[i]) *(r1[j] -r2[j]))/nu)*s1(t)  
    else:    
        nuc = nu * rr[i]*rr[j] * s3(t) + (2* exp1 *(r1[j]-r2[j])*rr[i] - 2*exp2 * (r1[i] - r2[i])*rr[j]) * s2(t) + (-1* (4*exp1*exp2*(r1[i] - r2[i]) *(r1[j] -r2[j]))/nu)*s1(t)  
    nuc *= -gab(exp1, exp2, r1, r2) * np.pi**1.5 / (4*nu**2) 
    nuc *= (128 *exp1**5 / np.pi**3)**0.25
    nuc *= (128 *exp2**5 / np.pi**3)**0.25
    return nuc

def nuc_ps(type1, type2, exp1, exp2, r1, r2, r3):
    nu = exp1 +  exp2
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - r3
    t = nu**0.5 * np.linalg.norm(rr)
    i  =  type1 - 1
    nuc = -gab(exp1, exp2, r1, r2)  * np.pi**1.5 / (2*nu) * (rr[i]*s2(t) - (2* exp2 * (r1[i] - r2[i])) / nu * s1(t))
    nuc *= (128 *exp1**5 / np.pi**3)**0.25
    nuc *= (2.0 * exp2 / np.pi)**0.75 
    return nuc

#noo, orb_coord, orb_type, coef_mat, exp_mat, zoa
def get_nuc_mat(noo, orb_coord, orb_type, coef_mat, exp_mat, coord, zoa):
    nuc_mat = np.zeros((noo, noo))
    for i in range(noo):
        for j in range(noo):
            #if (i == 6 and j == 0) or (i == 0 and j == 6):
            # print("jetzt!!!!!!!!!!!!!!", i,j,)
             for k in range(3):
                for l in range(3):
                    for m in range(zoa.shape[0]):
                    #print(get_nuc(orb_type[i], orb_type[j], exp_mat[i,k], exp_mat[j,l], orb_coord[i,:], orb_coord[j,:]))
                        tmp1 = zoa[m] * coef_mat[i,k] *  coef_mat[j,l] * get_nuc(orb_type[i], orb_type[j], exp_mat[i,k], exp_mat[j,l], orb_coord[i,:], orb_coord[j,:], coord[m,:]) 
                        #if m == 0:
                        # print(i,j, k, l, m, tmp1)
                        ## print(orb_coord[i,:], orb_coord[j,:], coord[m,:])
                        nuc_mat[i,j] += tmp1
           # for m in range(zoa.shape[0]):
           #         nuc_mat[i,j]
    return nuc_mat
