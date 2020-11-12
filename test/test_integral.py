from bombanitio import eval_over
from bombanitio import eval_kin
from bombanitio import eval_nuc
from bombanitio import eval_four
from bombanitio import read_xyz
import pytest
from bombanitio import create_basis
import numpy as np

def test_over():
    coord, atom = read_xyz.easy_read("test_coord.xyz", None, False, False)
    coord = coord[0,:,:]
    coord = read_xyz.xyz_to_bohr(coord)
    zoa = read_xyz.atom_to_zoa(atom)
    noa = atom.shape[0]
    noo, orb_coord, orb_type, coef_mat, exp_mat = create_basis.basisset(noa, zoa, coord)
    #read_xyz.plot_initial_orbs(noo, orb_coord, orb_type, coef_mat, exp_mat)
    over_mat = eval_over.get_overlap_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    over_ref =  np.array([[1.00000,   0.23670,   0.00000,   0.00000,   0.00000,   0.05397,   0.05397],
    [0.23670,   1.00000,   0.00000,   0.00000,   0.00000,   0.47475,   0.47475],
    [0.00000,   0.00000,   1.00000,   0.00000,   0.00000,   0.31112,  -0.31112],
    [ 0.00000,   0.00000,   0.00000,  1.00000,   0.00000,   0.24081,   0.24081],
    [ 0.00000,   0.00000,   0.00000,   0.00000,   1.00000,   0.00000,   0.00000],
    [ 0.05397,   0.47475,   0.31112,   0.24081,   0.00000,   1.00000,   0.25167],
    [ 0.05397,   0.47475,  -0.31112,   0.24081,   0.00000,   0.25167,   1.00000]])
    assert over_mat == pytest.approx(over_ref.astype('float64'), rel = 0.001)

def test_kin():
    coord, atom = read_xyz.easy_read("test_coord.xyz", None, False, False)
    coord = coord[0,:,:]
    coord = read_xyz.xyz_to_bohr(coord)
    zoa = read_xyz.atom_to_zoa(atom)
    noa = atom.shape[0]
    noo, orb_coord, orb_type, coef_mat, exp_mat = create_basis.basisset(noa, zoa, coord)
    kin_mat = eval_kin.get_kin_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    kin_ref =  np.array([[29.00324,  -0.16801,   0.00000,   0.00000,   0.00000,  -0.00252,  -0.00252],
    [-0.16801,   0.80813,   0.00000,   0.00000,   0.00000,   0.12857,   0.12857],
    [ 0.00000,  0.00000,   2.52873,   0.00000,  0.00000,   0.22476,  -0.22476],
    [ 0.00000,   0.00000,   0.00000,   2.52873,   0.00000,   0.17396,   0.17396],
    [ 0.00000,   0.00000,   0.00000,   0.00000,   2.52873,   0.00000,   0.00000],
    [-0.00252,   0.12857,   0.22476,   0.17396,   0.00000,   0.76003,   0.00848],
    [ -0.00252,   0.12857,  -0.22476,   0.17396,   0.00000,   0.00848,   0.76003]])
    #assert kin_mat == pytest.approx(kin_ref.astype('float64'), rel = 0.01)
    #np.set_printoptions(precision=4, suppress=True)
    #print(kin_mat-kin_ref)
    #print()
    #print(kin_mat)
    #print()
    #print(kin_ref)
    #print()
    #print((kin_mat-kin_ref)/kin_mat)
    #print()
    assert kin_mat == pytest.approx(kin_ref.astype('float64') , rel=1e-3, abs=1e-4)






def test_nuc():
    coord, atom = read_xyz.easy_read("test_coord.xyz", None, False, False)
    coord = coord[0,:,:]
    coord = read_xyz.xyz_to_bohr(coord)
    zoa = read_xyz.atom_to_zoa(atom)
    noa = atom.shape[0]
    noo, orb_coord, orb_type, coef_mat, exp_mat = create_basis.basisset(noa, zoa, coord)
#NUCLEAR-ATTRACTION MATRIX
    nuc_mat = eval_nuc.get_nuc_mat(noo, orb_coord, orb_type, coef_mat, exp_mat, coord, zoa)
    #print("nuclear attraction  matrix: ")
    nuc_ref = np.array([[-61.72413,  -7.44477,   0.00000,  -0.01900,   0.00000,  -1.74672,  -1.74672],
 [-7.44477, -10.14287,   0.00000,  -0.22350,   0.00000,  -3.86944,  -3.86944],
 [ 0.00000,   0.00000, -10.14262,   0.00000,   0.00000,  -2.25477,   2.25477],
 [-0.01900,  -0.22350,   0.00000, -10.07996,   0.00000,  -1.81819,  -1.81819],
 [ 0.00000,   0.00000,   0.00000,   0.00000,  -9.98632,   0.00000,   0.00000],
 [-1.74672,  -3.86944,  -2.25477,  -1.81819,   0.00000,  -5.83642,  -1.61650],
 [-1.74672,  -3.86944,   2.25477,  -1.81819,   0.00000,  -1.61650,  -5.83642]])
    
    np.set_printoptions(precision=4, suppress=True)
    print(nuc_ref - nuc_mat)
    print(nuc_mat)
    assert nuc_mat == pytest.approx(nuc_ref.astype('float64'), rel = 0.001)

def test_four_ss():
    coord, atom = read_xyz.easy_read("test_coord.xyz", None, False, False)
    coord = coord[0,:,:]
    coord = read_xyz.xyz_to_bohr(coord)
    zoa = read_xyz.atom_to_zoa(atom)
    noa = atom.shape[0]
    noo, orb_coord, orb_type, coef_mat, exp_mat = create_basis.basisset(noa, zoa, coord)
    four_mat = eval_four.get_four_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    assert four_mat[0,0,0,0] == pytest.approx(4.78508, rel = 0.001) #ssss
    assert four_mat[0,1,0,0] == pytest.approx(0.74138, rel = 0.001) #ssss
    assert four_mat[5,5,5,3] == pytest.approx(0.19847, rel = 0.001) #psss
    assert four_mat[6,6,6,3] == pytest.approx(0.19847, rel = 0.001) #psss
    assert four_mat[6,6,6,2] == pytest.approx(-0.25641, rel = 0.001)#psss
#0.19847
#7 7 5 5     0.48275
    assert four_mat[6,6,4,4] == pytest.approx(0.48275, rel = 0.001)#ppss
    assert four_mat[6,6,4,3] == pytest.approx(0.0, rel = 0.001)#ppss
#7 5 7 5     0.02343
    assert four_mat[6,4,6,4] == pytest.approx(0.02343, rel = 0.001)#psps
#test_nuc()
#7 4 3 3     0.15098
##7 4 4 1     0.01162
##7 4 4 2     0.07176
#7 4 4 3    -0.01229
    assert four_mat[6,3,2,2] == pytest.approx(0.15098, rel = 0.001)#ppps
    assert four_mat[6,3,3,2] == pytest.approx(-0.01229, rel = 0.001)#ppps
#4 4 4 4     0.88016
    assert four_mat[3,3,3,3] == pytest.approx(0.88016, rel = 0.001)#ppps
#def test_four_full():
    four_pre = np.loadtxt("../data/FOUR_CENTRE_TEST_DATA.dat")
    for i in range(four_pre.shape[0]):
#    for i in range(3):
#        print("hallllllllllllllllllllloooooooooooooooooooo")
        print(i)
        #assert four_mat[int(four_pre[i,0])-1,int(four_pre[i,1])-1,int(four_pre[i,2])-1,int(four_pre[i,3])-1] == pytest.approx(four_pre[i,4], rel = 0.05)
        assert four_mat[int(four_pre[i,0])-1,int(four_pre[i,1])-1,int(four_pre[i,2])-1,int(four_pre[i,3])-1] == pytest.approx(four_pre[i,4], rel=1e-3, abs=1e-5)






#When subroutines SMATRIX, KMATRIX and VMATRIX are properly coded, 
#the test program TWO_CENTRE_TEST should give the following output print:
#

#TEST OF SUBROUTINES SMATRIX, KMATRIX AND VMATRIX
#MOLECULE: WATER


#OVERLAP MATRIX
#  1.00000   0.23670   0.00000   0.00000   0.00000   0.05397   0.05397
#  0.23670   1.00000   0.00000   0.00000   0.00000   0.47475   0.47475
#  0.00000   0.00000   1.00000   0.00000   0.00000   0.31112  -0.31112
#  0.00000   0.00000   0.00000   1.00000   0.00000   0.24081   0.24081
#  0.00000   0.00000   0.00000   0.00000   1.00000   0.00000   0.00000
#  0.05397   0.47475   0.31112   0.24081   0.00000   1.00000   0.25167
#  0.05397   0.47475  -0.31112   0.24081   0.00000   0.25167   1.00000
#
#
#KINETIC-ENERGY MATRIX
# 29.00324  -0.16801   0.00000   0.00000   0.00000  -0.00252  -0.00252
# -0.16801   0.80813   0.00000   0.00000   0.00000   0.12857   0.12857
#  0.00000   0.00000   2.52873   0.00000   0.00000   0.22476  -0.22476
#  0.00000   0.00000   0.00000   2.52873   0.00000   0.17396   0.17396
#  0.00000   0.00000   0.00000   0.00000   2.52873   0.00000   0.00000
# -0.00252   0.12857   0.22476   0.17396   0.00000   0.76003   0.00848
# -0.00252   0.12857  -0.22476   0.17396   0.00000   0.00848   0.76003
#
#
#NUCLEAR-ATTRACTION MATRIX
#-61.72413  -7.44477   0.00000  -0.01900   0.00000  -1.74672  -1.74672
# -7.44477 -10.14287   0.00000  -0.22350   0.00000  -3.86944  -3.86944
#  0.00000   0.00000 -10.14262   0.00000   0.00000  -2.25477   2.25477
# -0.01900  -0.22350   0.00000 -10.07996   0.00000  -1.81819  -1.81819
#  0.00000   0.00000   0.00000   0.00000  -9.98632   0.00000   0.00000
# -1.74672  -3.86944  -2.25477  -1.81819   0.00000  -5.83642  -1.61650
# -1.74672  -3.86944   2.25477  -1.81819   0.00000  -1.61650  -5.83642

