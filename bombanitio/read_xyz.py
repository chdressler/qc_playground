import numpy as np
#import ipdb
import os
import time


def plot_initial_orbs(noo, orb_coord, orb_type, coef_mat, exp_mat):
    bohr_in_a = 0.529177210903
    np.set_printoptions(suppress=True)
    orb_coord = orb_coord * bohr_in_a
    #numpy.set_printoptions(precision=2)
    print("orbital number    center of orbital in Angstrom")
    for i in range(noo):
        print(i, "                ", orb_coord[i])
    print()
    print("orbital number    center of orbital in au")
    #np.set_printoptions(precision=3)
    for i in range(noo):
        print(i, "                ", orb_coord[i]/bohr_in_a)
    print()
    print("orbital number    position of orbital")
    for i in range(noo):
        print(i, "                    ", orb_type[i])
    print()
    #print("orbital number    position of orbital in au")
    #for i in range(noo):
    #    print(i, "                    ", orb_type[i]/bohr_in_a)
    #print()
    print("orbital number    contraction exponents of orbital")
    for i in range(noo):
        #np.set_printoptions(suppress=True)
        #print(exp_mat)
        print(i, "          ", exp_mat[i].astype('float64'))
    print()
    print("orbital number    contraction coefficients of orbital")
    for i in range(noo):
        print(i, "          ", coef_mat[i].astype('float64'))
    print()
        #print("orbital number    orbital centra")
        #ORBITAL NUMBER     ORBITAL CONTRACTION EXPONENTS


def xyz_to_bohr(coord):
    bohr_in_a = 0.529177210903
    return coord/bohr_in_a


def read_trajec(line, coords, atoms, noa):
    i = len(coords) -1
    tmp1 = line.rstrip('\n').split(" ")
    while '' in tmp1: tmp1.remove('')
    if not tmp1:
        pass
    elif tmp1[0] == str(noa):
        i += 1
        coords.append([])
    elif tmp1[0] == 'i':
        pass
    else:
        atoms.append(tmp1[0])
        coords[i].append([float(tmp1[1]),float(tmp1[2]),float(tmp1[3])])
    return coords, atoms


def traj_to_ar(path1, noa):
 atoms = []
 coords = []
# with open(path1) as f:
#    first_line = f.readline()
# noa = int(first_line.strip())
# print("number of detected atoms is: ", noa)
 with open(path1) as f:
     for line in f:
        coords, atoms = read_trajec(line, coords, atoms, noa)

 atom_list = atoms[:noa]
 atom_ar = np.array([atom_list])
 coord_ar = np.array(coords)
 atom_ar = atom_ar[0,:]
 return coord_ar, atom_ar

def wrap_to_zero(coord, pbc_mat):
    for i in range(1,coord.shape[0], 1):
        disp = coord[i,:,:] - coord[i-1,:,:]
        #for j in range(coord.shape[2]):
        #   #disp = coord[i,j,:] - coord[i-1,j,:]
        #   if disp[j
        #tmp1 = np.where(disp[:,0] > pbc_mat[0,0]/2) 
        #coord[i,:,0][tmp1] -= pbc_mat[0,0]
        #tmp2 = np.where(disp[:,0] < -pbc_mat[0,0]/2) 
        #coord[i,:,0][tmp1] += pbc_mat[0,0]
        for k in range(3):
            tmp1 = np.where(disp[:,k] > pbc_mat[k,k]/2)
            coord[i,:,k][tmp1] -= pbc_mat[k,k]
            tmp2 = np.where(disp[:,k] < -pbc_mat[k,k]/2)
            coord[i,:,k][tmp1] += pbc_mat[k,k]
    return coord



def remove_com(coord, atoms, pbc_mat, zero = True):
   mass_dict = dict()
   mass_dict["Si"] = 28
   mass_dict["Li"] = 0
   mass_dict["O"] = 16
   mass_dict["H"] = 1
   mass_dict["C"] = 12
   mass_dict["N"] = 1
   mass_dict["Na"] = 23
   mass_dict["F"] = 1
   mass_dict["S"] = 32
   for i in mass_dict.keys():
       if mass_dict[i] == 0:
          print("WARNING mass of "+ i +" is equal to zero.")
   mges = 0
   for j in range(coord.shape[1]):
       mges += mass_dict[atoms[j]]
   #for i  in range(1,1,coord.shape[0]):
   #for i  in range(1,coord.shape[0],1):
   for i  in range(0,coord.shape[0]):
       #ii#
       #com = 0
       #ipdb.set_trace()
       com = [0,0,0]
       for j in range(coord.shape[1]):
            com += coord[i,j,:] * mass_dict[atoms[j]]
       com /= mges
      #if i == 1:
       if i == 0:
           if zero == True:
               com_init =  np.array([0.0,0.0,0.0])
           else:
               com_init = com
       #else:
       #    coord[i,:,:]  = coord[i,:,:]  - com + com_init
       coord[i,:,:]  = coord[i,:,:]  - com + com_init
   return coord






def easy_read(path1, pbc, com = True, wrap = True):
#def easy_read(path1, pbc, com, wrap) :
    """reads trajectory and check if fast npz exist, if not npz file is created

    input: 
         string path1, 
         number of atoms noa,
         2-dimensional array with pbcs pbc[xaxis,yaxis,zaxis],
         bool com remove center of mass?,
         bool wrap wrap molecules into box?
    output:
         3-dimensional array of coordinates coord[time step, atom type, xyz coordinates],
         1-dimensional array with atom types
    """
    #auxiliary_file = path1 + ".npz"
    print("new")
    start_time = time.time()
    auxiliary_file = path1 + "_"+ str(com) + "_" +str(wrap) + "_"+ ".npz"
    with open(path1) as f:
        first_line = f.readline()
    noa = int(first_line.strip())
    print("number of detected atoms is: ", noa)
    if os.path.exists(auxiliary_file) == False:
        print("auxiliary file " + auxiliary_file  +" does not exist")
        coord, atom = traj_to_ar(path1, noa)
        nof = coord.shape[0]
        print("number of frames: " +str(nof))
        print("--- %s seconds ---" % (time.time() - start_time))
        if wrap:
            print("wrap")
            coord = wrap_to_zero(coord, pbc)
            print("--- %s seconds ---" % (time.time() - start_time))
        if com:
            print("remove com")
            coord = remove_com(coord, atom, pbc)
            print("--- %s seconds ---" % (time.time() - start_time))
        np.savez(auxiliary_file, coord, atom)
    else:
        print("auxiliary file " + auxiliary_file  +"  exists")
        ar2 = np.load(auxiliary_file)
        coord = ar2['arr_0']
        atom = ar2['arr_1']
    print("--- %s seconds ---" % (time.time() - start_time))
    return coord, atom


def atom_to_zoa(atom):
    dict1={}
    dict1["H"] = 1
    dict1["He"] = 2
    dict1["Li"] = 3
    dict1["Be"] = 4
    dict1["B"] = 5
    dict1["C"] = 6
    dict1["N"] = 7
    dict1["O"] = 8
    dict1["F"] = 9
    dict1["Ne"] = 10
    dict1["Na"] = 11
    dict1["Mg"] = 12
    dict1["Al"] = 13
    dict1["Si"] = 14
    dict1["P"] = 15
    dict1["S"] = 16
    dict1["Cl"] = 17
    dict1["Ar"] = 18
    zoa = []
    for i in list(atom):
        zoa.append(dict1[i])
    return np.array(zoa)
