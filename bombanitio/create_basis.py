# coding: utf-8
#get_ipython().system('cat basisset.py')
#import bse
import basis_set_exchange as bse
#get_ipython().run_line_magic('pinfo', 'bse.get_basis')
#bse.get_basis("sto-3g")
#get_ipython().run_line_magic('pinfo', 'bse.get_basis')
#bse.get_basis("sto-3g", elements="H")
#bse.get_basis("sto-3g", elements="H")["exponents"]
#bse.get_basis("sto-3g", elements="H")["elements"]
#bse.get_basis("sto-3g", elements="H")["elements"]["1"]
#bse.get_basis("sto-3g", elements="H")["elements"]["1"]["electron_shells"]
#bse.get_basis("sto-3g", elements="H")["elements"]["1"]["electron_shells"]["exponents"]
#bse.get_basis("sto-3g", elements="H")["elements"]["1"]["electron_shells"][0]["exponents"]
#exps = bse.get_basis("sto-3g", elements="H")["elements"]["1"]["electron_shells"][0]["exponents"]
import numpy as np
#np.fromstring(exps)
#np.asfarray(exps)



def orbitals(j, z):
    """[input] j: orbital type, z: atomic number, [outpout] coef: contraction coefficients, exp: contraction exponents"""
    #1: 1s, 2:2s, 3:2px, 4:2py, 5:2pz, 6:3s, 7:3px, 8:3py, 9:3pz
    if j == 1:
    #if j <= 2:
        exp = bse.get_basis("sto-3g")["elements"][str(z)]["electron_shells"][0]["exponents"]
        coef = bse.get_basis("sto-3g")["elements"][str(z)]["electron_shells"][0]["coefficients"][0]
    elif j > 1 and j <= 2:
        exp = bse.get_basis("sto-3g")["elements"][str(z)]["electron_shells"][1]["exponents"]
        coef = bse.get_basis("sto-3g")["elements"][str(z)]["electron_shells"][1]["coefficients"][0]
    elif j > 2 and j <= 5:
        exp = bse.get_basis("sto-3g")["elements"][str(z)]["electron_shells"][1]["exponents"]
        coef = bse.get_basis("sto-3g")["elements"][str(z)]["electron_shells"][1]["coefficients"][1]
    elif j > 5 and j <= 6:
        exp = bse.get_basis("sto-3g")["elements"][str(z)]["electron_shells"][2]["exponents"]
        coef = bse.get_basis("sto-3g")["elements"][str(z)]["electron_shells"][2]["coefficients"][0]
    elif j > 6 and j <= 9:
        exp = bse.get_basis("sto-3g")["elements"][str(z)]["electron_shells"][2]["exponents"]
        coef = bse.get_basis("sto-3g")["elements"][str(z)]["electron_shells"][2]["coefficients"][1]
    else:
        print("error: orbital number does not exist")
    return np.array(coef), np.array(exp) 


if __name__ == '__main__':
   #print("hello123")
   #elements=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
   #elements=[10,11,12,13,14,15,16,17,18]
   elements=[11,12,13,14,15,16,17,18]
   #elements=[1,2,3]
   orbs = [1,2,3,4,5,6,7,8,9]
   #orbs = [1]
   for i in elements:
       for j in orbs:
           print(i,j)
           print(orbitals(j, i))

#def orbitals(j, z):

     
#     return a,d

def get_noo(zoa):
    noo = 0
    for i in list(zoa):
        if i <= 2:
            noo += 1
        elif i > 2 and i <= 10: 
            noo += 5
        elif i > 10 and i <= 18:
            noo += 9
        else: 
            print("error: nucleus not supported")
    return noo


def basisset(noa, zoa, coord):
    """[input] noa: number of atoms, zoa: 1d array of atomic number, coord: 3darray of  coords [outpout] noo: number of orbitals, orb_coord: 1d array assingment of nuclei to orbitals , coeff_mat: 2d mat, matrix with contraction coefficients in each row, exp_mat: matrix containing in eacht row the contraction exponents"""
    #lqm= {}
    #qml[1] = 1
    #qml[2] = 1
    #qml[3] = 5
    #qml[4] = 5
    #qml[5] = 5
    #qml[6] = 5
    #qml[8] = 5
    #qml[9] = 5
    #qml[11] = 9
    #qml[12] = 9
    #qml[13] = 9
    #qml[14] = 9
    #qml[15] = 9
    #qml[16] = 9
    #qml[17] = 9
    #qml[18] = 9
    #noo = 0
    #for i in list(zoa):
    #    if i <= 2:
    #        noo += 1
    #    elif i > 2 and i <= 10: 
    #        noo += 5
    #    elif i > 10 and i <= 18:
    #        noo += 9
    #    else: 
    #        print("error: nucleus not supported")
    noo = get_noo(zoa)
    coef_mat = []
    exp_mat = []
    orb_coord = []
    orb_type = [] # 0:s,  1:px, 2:py, 3pz  
    for i, ele  in enumerate(list(zoa)):
        if ele <= 2:
            qml = 1
        elif ele > 2 and ele <= 10: 
            qml = 5 
        elif ele > 10 and ele <= 18:
            qml = 9
        else: 
            print("error: nucleus not supported")
        for k in range(qml):
            coef, exp = orbitals(k+1, ele)
            coef_mat.append(coef)
            exp_mat.append(exp)
            orb_coord.append(coord[i,:])
            if k+1  == 1 or k+1 == 2 or k+1  ==  6:
                orb_temp = 0
            elif k+1  == 3 or k+1 == 7:
                orb_temp = 1
            elif k+1  == 4 or k+1 == 8:
                orb_temp = 2
            elif k+1  == 5 or k+1 == 9:
                orb_temp = 3
            else:
                print("error: either s nor p orbital")
            orb_type.append(orb_temp)
    coef_mat = np.array(coef_mat)
    exp_mat = np.array(exp_mat)
    orb_coord = np.array(orb_coord)
    orb_type = np.array(orb_type)
    return noo, orb_coord, orb_type, coef_mat.astype('float64'), exp_mat.astype('float64')





        

#        for k in range(5):

#
#        if ele  <= 2:
#            #for k in range(
#            coef, exp = orbitals(1, ele)
#            coef_mat.append(coef)
#            exp_mat.append(exp)
#            orb_coord.append(coord[i,:])
 #           orb_type.append(0)
#        elif ele > 2 and ele <= 10:  
#            for k in range(5):
#                if k == 0:
#                    coef, exp = orbitals(1, ele)
#                    coef_mat.append(coef)
#                    exp_mat.append(exp)
###                    orb_coord.append(coord[i,:])
  #                  orb_type.append(0)
  #              elif k == 1:
  #                  coef, exp = orbitals(2, ele)
  #                  coef_mat.append(coef)
  #                  exp_mat.append(exp)
  #                  orb_coord.append(coord[i,:])
  #                  orb_type.append(0)
   ##             elif k > 1 and k < 5:
   #                 coef, exp = orbitals(k+1, ele)
    ##                coef_mat.append(coef)
    ##                exp_mat.append(exp)
    #                orb_coord.append(coord[i,:])
     #           
##
 #   return orb_coord, orb_type, coef_mat, exp_mat
    











































