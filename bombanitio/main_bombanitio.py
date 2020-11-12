import numpy as np
from  bombanitio import read_xyz
from bombanitio import create_basis
#import flit
from bombanitio import eval_over
from bombanitio import eval_kin
from bombanitio import eval_nuc
from bombanitio import eval_four
from bombanitio import eval_dipole
from bombanitio import scf




def main():
    print("brrrr")
    print("Performing a  self-consistent HF calculation of the molecule given by the file  coord.xyz in this directory ... ")
    print("Reading coordinates ... ")
    coord, atom = read_xyz.easy_read("coord.xyz", None, False, False)
    coord = coord[0,:,:]
    coord = read_xyz.xyz_to_bohr(coord)
    zoa = read_xyz.atom_to_zoa(atom)
    noa = atom.shape[0] 
    noo, orb_coord, orb_type, coef_mat, exp_mat = create_basis.basisset(noa, zoa, coord)
    read_xyz.plot_initial_orbs(noo, orb_coord, orb_type, coef_mat, exp_mat)

    over_mat = eval_over.get_overlap_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    print("overlap matrix: ")
    np.set_printoptions(precision=4, suppress=True)
    print(over_mat.astype('float64'))
    

    kin_mat = eval_kin.get_kin_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    print("kinetic energy  matrix: ")
    print(kin_mat.astype('float64'))

    nuc_mat = eval_nuc.get_nuc_mat(noo, orb_coord, orb_type, coef_mat, exp_mat, coord, zoa)
    print("nuclear attraction  matrix: ")
    print(nuc_mat.astype('float64'))

    dip_matx = eval_dipole.get_dip_mat(noo, orb_coord, orb_type, coef_mat, exp_mat, 0)
    dip_maty = eval_dipole.get_dip_mat(noo, orb_coord, orb_type, coef_mat, exp_mat, 1)
    dip_matz = eval_dipole.get_dip_mat(noo, orb_coord, orb_type, coef_mat, exp_mat, 2)
    dip_mat = np.array([dip_matx, dip_maty, dip_matz])
    print("dipole moment matrix: ")
    print(dip_mat.astype('float64'))
    
    print("calculate four center  integrals.... ")
    four_mat = eval_four.get_four_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    #print(four_mat[0,0,0,0])
    #print(four_mat[5,5,5,3])
    
    #create core hamiltonian 
    h =  kin_mat + nuc_mat

    x = scf.ortho_basis(over_mat)
    noe = zoa.sum()
    print("number of electron:", noe) 
    en, p, coeff, energy_orb = scf.scf(four_mat, noe, noo, x, h)
    print()
    print("#######################################################")
    print("final orbital energies:", energy_orb)
    print("#######################################################")
    #print()
    print()
    print("#######################################################")
    print("final electronic energy:", en)
    print("#######################################################")
    print()

    etot = scf.calc_tot_en(noa, en, coord, zoa)
    print("#######################################################")
    print("final total energy:", etot)
    print("#######################################################")
    #print(etot)
    dipole = eval_dipole.get_el_dipole(p, dip_mat)
    print()
    print()
    print()
    print("#######################################################")
    print("electronic contribution dipole moment:", dipole)
    nuc_pole = eval_dipole.get_nuc_pole(coord, zoa)
    print("nucear contribution dipole moment:", nuc_pole)
    print("total dipole moment:", nuc_pole - dipole, np.linalg.norm(nuc_pole - dipole))    
    print("#######################################################")
    
    np.savetxt("coefficients", coeff[:int(noe/2),:]) 



