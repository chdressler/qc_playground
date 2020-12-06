import sys
import numpy as np
import scipy.special
from bombanitio.special import specialx
import sys

def get_four(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    if type1 == 0 and type2 ==0 and type3 == 0 and type4 == 0:
        four = four_ssss(exp1, exp2, exp3, exp4,  r1, r2, r3, r4)

    elif type1 > 0  and type2 > 0 and type3 > 0 and type4 > 0:
        four = four_pppp(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)
    
    elif type1 > 0  and type2 ==0 and type3 == 0 and type4 == 0:
        four = four_psss(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)
    elif type1 == 0 and type2 > 0 and type3 == 0 and type4 == 0:
        four = four_psss(type2, type1, type3, type4,  exp2, exp1, exp3, exp4,  r2, r1, r3, r4)
    elif type1 == 0 and type2 == 0 and type3 > 0 and type4 == 0:
        #four = four_psss(type3, type2, type1, type4,  exp3, exp2, exp1, exp4,  r3, r2, r1, r4)
        four = four_psss(type3, type4, type1, type2,  exp3, exp4, exp1, exp2,  r3, r4, r1, r2)
    elif type1 == 0  and type2 ==0 and type3 == 0 and type4 > 0:
        #four = four_psss(type4, type2, type3, type1,  exp4, exp2, exp3, exp1,  r4, r2, r3, r1)
        four = four_psss(type4, type3, type1, type2,  exp4, exp3, exp1, exp2,  r4, r3, r1, r2)
    
   
    elif type1 > 0  and type2 > 0 and type3 == 0 and type4 == 0:
        four = four_ppss(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)
    elif type1 == 0  and type2 == 0 and type3 > 0 and type4 > 0:
        four = four_ppss(type3, type4, type1, type2,  exp3, exp4, exp1, exp2,  r3, r4, r1, r2)


    elif type1 > 0  and type2 == 0 and type3 > 0 and type4 == 0:
        four = four_psps(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)
    elif type1 == 0  and type2 >  0 and type3 > 0 and type4 == 0:
        four = four_psps(type2, type1, type3, type4,  exp2, exp1, exp3, exp4,  r2, r1, r3, r4)
    elif type1 > 0  and type2 == 0 and type3 ==  0 and type4 > 0:
        four = four_psps(type1, type2, type4, type3,  exp1, exp2, exp4, exp3,  r1, r2, r4, r3)
    elif type1 == 0  and type2 > 0 and type3 ==  0 and type4 > 0:
        four = four_psps(type2, type1, type4, type3,  exp2, exp1, exp4, exp3,  r2, r1, r4, r3)
  
   
    elif type1 > 0  and type2 > 0 and type3 > 0 and type4 == 0:
        four = four_ppps(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)
    elif type1 > 0  and type2 > 0 and type3 == 0 and type4 > 0:
        four = four_ppps(type1, type2, type4, type3,  exp1, exp2, exp4, exp3,  r1, r2, r4, r3)
    elif type1 > 0  and type2 == 0 and type3 > 0 and type4 > 0:
        #four = four_ppps(type1, type4, type3, type2,  exp1, exp4, exp3, exp2,  r1, r4, r3, r2)
        four = four_ppps(type3, type4, type1, type2,  exp3, exp4, exp1, exp2,  r3, r4, r1, r2)
    elif type1 == 0  and type2 > 0 and type3 > 0 and type4 > 0:
        #print("hallo korr type")
        #four = four_ppps(type4, type2, type3, type1,  exp4, exp2, exp3, exp1,  r4, r2, r3, r1)
        four = four_ppps(type4, type3, type2, type1,  exp4, exp3, exp2, exp1,  r4, r3, r2, r1)

        
    else:
        print("error orbital combination does not exist!")
        sys.exit()
    return four


        


#specialx.specialfunc(3,1)
#def nu(exp1, exp2):
#    return exp1 +  exp2
#def t
#
#nu = exp1 +  exp2
#rr = exp1 * r1 + exp2 * r2 / ( exp1 + exp2) - r3
#t = nu**0.5 * np.linalg.norm(rr)

#np.longdouble


def delta(i,j):
    #if i.is_integer() and j.is_integer():
        if i == j: 
            k = 1
        else:
            k = 0
    #else:
    #    print(" i and j has to be integers")
    #    sys.exit()
        return k

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

def s4(t):
    return specialx.specialfunc(4,t)[3]


def s5(t):
    return specialx.specialfunc(5,t)[4]


def gab(exp1, exp2, r1, r2):
    return np.exp(- (exp1 * exp2)/(exp1 +exp2) * np.linalg.norm(r1-r2)**2)

def four_ssss(exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4) 
    t = q * np.linalg.norm(rr)
    four = gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) * q * np.pi**3 / (nu *  theta)**1.5  *  s1(t)
    four *= (2.0 * exp1 / np.pi)**0.75 
    four *= (2.0 * exp2 / np.pi)**0.75 
    four *= (2.0 * exp3 / np.pi)**0.75 
    four *= (2.0 * exp4 / np.pi)**0.75 
    #print(exp1, exp2, r1,r2, r3, four, nu,rr, t)
    return four


def four_psss(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    four = gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) / 2  * q * np.pi**3 / (nu**2.5 *  theta**1.5 ) * (q**2 * rr[i] * s2(t)  - 2* exp2 * (r1[i] - r2[i]) * s1(t))
    four *= (128 *exp1**5 / np.pi**3)**0.25
    four *= (2.0 * exp2 / np.pi)**0.75 
    four *= (2.0 * exp3 / np.pi)**0.75
    four *= (2.0 * exp4 / np.pi)**0.75
    return four



def four_ppss(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    #four = gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) / 2  * q * np.pi**3 / (nu**2.5 *  theta**1.5 ) * (q**2 * rr[i] * s2(t)  - 2* exp2 * (r1[i] - r2[i]) * s1(t))
    four = gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) / 4  * q * np.pi**3 / (nu**3.5 *  theta**1.5 ) * (q**4 * rr[i] * rr[j] * s3(t) + q**2 * ( delta(i,j) + 2 * exp1 * (r1[j] - r2[j])* rr[i] - 2 * exp2 * (r1[i] - r2[i]) * rr[j]) * s2(t) +  (2*nu*delta(i,j) -4*exp1*exp2*(r1[i] - r2[i]) * (r1[j] - r2[j]))*s1(t)    )
    four *= (128 *exp1**5 / np.pi**3)**0.25
    four *= (128 *exp2**5 / np.pi**3)**0.25
    four *= (2.0 * exp3 / np.pi)**0.75
    four *= (2.0 * exp4 / np.pi)**0.75
    return four



def four_psps(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4, norm = True):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    #four = gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) / 2  * q * np.pi**3 / (nu**2.5 *  theta**1.5 ) * (q**2 * rr[i] * s2(t)  - 2* exp2 * (r1[i] - r2[i]) * s1(t))
    four = gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) / 4  * q * np.pi**3 / (nu**2.5 *  theta**2.5 ) * (-q**4 * rr[i] * rr[k] * s3(t) + q**2 * (2*exp2 *(r1[i]-r2[i])* rr[k]- 2*exp4*(r3[k]-r4[k])*rr[i]-delta(i,k)) *s2(t) + 4 *exp2 * exp4 * (r1[i] -r2[i])*(r3[k]-r4[k])*s1(t))           
    if norm:
        four *= (128 *exp1**5 / np.pi**3)**0.25
        four *= (2.0 * exp2 / np.pi)**0.75
        four *= (128 *exp3**5 / np.pi**3)**0.25
        four *= (2.0 * exp4 / np.pi)**0.75
    return four

def u1(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    tmp = -4 * nu * exp4 * (r3[k]-r4[k]) * delta(i,j)
    return tmp



def u2(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    tmp  = 2 *(2*exp2 *exp4 * (r1[i] - r2[i])*(r3[k] -r4[k])*rr[j] + exp2 * (r1[i] -r2[i]) *  delta(j,k) - nu * rr[k]*delta(i,j) -  exp4 * (r3[k]-r4[k]) * delta(i,j)) 
    #tmp = -  exp4 * (r3[k]-r4[k]) * delta(i,j)
    return tmp


def u3(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    l = type4 - 1
    #tmp  =  2 * theta * rr[i] * rr[j] * delta(k,l) +  2 * exp4 * (r3[k] - r4[k]) * (rr[i]* delta(j,l) + rr[j] * delta(i,l)) - 2 * exp2 * (r1[i] - r2[i])*(rr[k]*delta(j,l) +rr[j]*delta(k,l)) + delta(i,j)*delta(k,l) + delta(i,l) * delta(j,k) + delta(i,k) * delta(j,l) 
    tmp  =  2 * exp2 * (r1[i] - r2[i]) * rr[j] * rr[k] - 2 * exp4 * (r3[k]- r4[k])*rr[i]*rr[j] - rr[k]*delta(i,j) - rr[i] * delta(j,k) - rr[j]*delta(i,k) 
    return tmp

def u4(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    l = type4 - 1
    #tmp  = rr[i] * rr[k] * delta(j,l) + rr[i] * rr[j] * delta(k,l) + rr[j] * rr[k] * delta(i,l)
    tmp  = - rr[i] * rr[j] * rr[k]
    return tmp

def four_ppps(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4, norm = True):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    l = type4 - 1
    #four = exp1 / nu * (r1[j] - r2[j]) * four_psps(type1, type4, type3, type4,  exp1, exp4, exp3, exp4,  r1, r4, r3, r4) + gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) / 8  * q * np.pi**3 / (nu**3.5 *  theta**2.5 ) * (q**6 * u4(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) * s4(t) + q**4 * u3(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) * s3(t) + q**2 * u2(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)  * s2(t) + u1(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) * s1(t))
    #four = exp1 / nu * (r1[j] - r2[j]) * four_psps(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) 
    four = exp1 / nu * (r1[j] - r2[j]) * four_psps(type1, 0, type3, 0,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4, norm = False) 
    #print(four)
    four +=  gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) / 8  * q * np.pi**3 / (nu**3.5 *  theta**2.5 ) * (q**6 * u4(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) * s4(t) + q**4 * u3(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) * s3(t) + q**2 * u2(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)  * s2(t) + u1(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) * s1(t))
    #print(u1(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4))
    #print(u2(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4))
    #print(u3(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4))
    #print(u4(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4))
    #four = gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) / 2  * q * np.pi**3 / (nu**2.5 *  theta**1.5 ) * (q**2 * rr[i] * s2(t)  - 2* exp2 * (r1[i] - r2[i]) * s1(t))
    #four = gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) / 4  * q * np.pi**3 / (nu**2.5 *  theta**2.5 ) * (-q**4 * rr[i] * rr[k] * s3(t) + q**2 * (2*exp2 *(r1[i]-r2[i])* rr[k]- 2*exp4*(r3[k]-r4[k])*rr[i]-delta(i,k)) *s2(t) + 4 *exp2 * exp4 * (r1[i] -r2[i])*(r3[k]-r4[k])*s1(t))
    #four = 
    if norm:
        four *= (128 *exp1**5 / np.pi**3)**0.25
        four *= (128 *exp2**5 / np.pi**3)**0.25
        four *= (128 *exp3**5 / np.pi**3)**0.25
        four *= (2.0 * exp4 / np.pi)**0.75
    #print(four)
    return four


def v1(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    l = type4 - 1
    tmp = 4* nu * theta * delta(i,j) * delta(k,l)
    return tmp



def v2(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    l = type4 - 1
    tmp = 2*((nu + theta) *delta(i,j)*delta(k,l) - 2 * exp2 * theta *(r1[i] - r2[i])*rr[j]*delta(k,l) - 2*exp2*exp4*(r1[i]-r2[i])*(r3[k]-r4[k])*delta(j,l))
    return tmp

def v3(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    l = type4 - 1
    tmp  =  2 * theta * rr[i] * rr[j] * delta(k,l) +  2 * exp4 * (r3[k] - r4[k]) * (rr[i]* delta(j,l) + rr[j] * delta(i,l)) - 2 * exp2 * (r1[i] - r2[i])*(rr[k]*delta(j,l) +rr[j]*delta(k,l)) + delta(i,j)*delta(k,l) + delta(i,l) * delta(j,k) + delta(i,k) * delta(j,l) 
    return tmp





def v4(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    l = type4 - 1
    tmp  = rr[i] * rr[k] * delta(j,l) + rr[i] * rr[j] * delta(k,l) + rr[j] * rr[k] * delta(i,l)
    #tmp  = - rr[i] * rr[j] * rr[k]
    return tmp



def four_pppp(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4):
    nu = exp1 +  exp2
    theta = exp3 + exp4
    q = ( nu * theta / (nu + theta))**0.5
    rr = (exp1 * r1 + exp2 * r2) / ( exp1 + exp2) - (exp3 * r3 + exp4 * r4) / ( exp3 + exp4)
    t = q * np.linalg.norm(rr)
    i  =  type1 - 1
    j = type2 - 1
    k = type3 - 1
    l = type4 - 1
    #four = exp1 / nu * (r1[j] - r2[j]) * four_ppps(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)
    #changed 20.11.20
    #four = exp1 / nu * (r1[j] - r2[j]) * four_ppps(type1, 0, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)
    #four = exp1 / nu * (r1[j] - r2[j]) * four_ppps(type3, type4, type1, 0,  exp3, exp4, exp1, exp2,  r3, r4, r1, r2)
    four = exp1 / nu * (r1[j] - r2[j]) * four_ppps(type3, type4, type1, 0,  exp3, exp4, exp1, exp2,  r3, r4, r1, r2, norm = False)
    four += gab(exp1, exp2, r1, r2) * gab(exp3, exp4, r3, r4) / 16  * q * np.pi**3 / (nu**3.5 *  theta**3.5 ) * (-q**8*u4(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) * rr[l]*s5(t) + q**6*(v4(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) + 2*exp3*(r3[l] - r4[l])*u4(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)- u3(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)*rr[l])*s4(t) + q**4 * (v3(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) +  2 * exp3 * (r3[l]-r4[l]) * u3(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) - u2(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)*rr[l]) * s3(t) + q**2 * (v2(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) +  2* exp3 * (r3[l] -r4[l])*u2(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) - u1(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)*rr[l])*s2(t) + (v1(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4) + 2* exp3*(r3[l] - r4[l])* u1(type1, type2, type3, type4,  exp1, exp2, exp3, exp4,  r1, r2, r3, r4)) * s1(t))
    four *= (128 *exp1**5 / np.pi**3)**0.25
    four *= (128 *exp2**5 / np.pi**3)**0.25
    four *= (128 *exp3**5 / np.pi**3)**0.25
    four *= (128 *exp4**5 / np.pi**3)**0.25
    #print(four)
    return four







#noo, orb_coord, orb_type, coef_mat, exp_mat, zoa
def get_four_mat(noo, orb_coord, orb_type, coef_mat, exp_mat):
    four_mat = np.zeros((noo, noo, noo, noo))
    for i in range(noo):
        print(i/noo) 
        for j in range(noo):
            for k  in range(noo):
                for l in range(noo):
                    #5553
                    #if (i ==0  and (j ==0 or j == 1) and k == 0 and l == 0) or (i == 5 and j == 5 and k ==5 and l == 3): 
                    #if (i ==0  and (j ==0 or j == 1) and k == 0 and l == 0) or (i == 5 and j == 5 and k ==5 and l == 3) or (i == 6 and j == 6 and k ==6 and (l == 3 or l ==2)): 
                    #6,3,2,2 6,3,3,2
                    #if (i ==0  and (j ==0 or j == 1) and k == 0 and l == 0) or (i == 5 and j == 5 and k ==5 and l == 3) or (i == 6 and j == 6 and k ==6 and (l == 3 or l ==2)) or ((i == 6 or i == 4) and ( j== 6) and ( k ==6) and (l == 4 or l ==3)) or ( i==6 and j == 4 and k ==6 and l==4) or ( i==6 and j == 6 and k ==4 and l==4) or (i == 6 and j == 3 and (k == 2 or k == 3) and l == 2):
                 
                     #if (i ==0  and (j ==0 or j == 1) and k == 0 and l == 0) or (i == 5 and j == 5 and k ==5 and l == 3) or (i == 6 and j == 6 and k ==6 and (l == 3 or l ==2)) or ((i == 6 or i == 4) and ( j== 6) and ( k ==6) and (l == 4 or l ==3)) or ( i==6 and j == 4 and k ==6 and l==4) or ( i==6 and j == 6 and k ==4 and l==4) or (i == 6 and j == 3 and (k == 2 or k == 3) and l == 2) or (i == 3 and j ==3 and k ==3 and l ==3):
                     for m in range(3):
                        for n in range(3):
                            for o  in range(3):
                                for p in range(3):
                                    four_mat[i,j,k,l] +=  coef_mat[i,m] *  coef_mat[j,n] * coef_mat[k,o] *  coef_mat[l,p] * get_four(orb_type[i], orb_type[j], orb_type[k], orb_type[l],  exp_mat[i,m], exp_mat[j,n], exp_mat[k,o], exp_mat[l,p], orb_coord[i,:], orb_coord[j,:], orb_coord[k,:], orb_coord[l,:]   ) 
                                    #print("hello", four_mat[i,j,k,l])
    return four_mat
