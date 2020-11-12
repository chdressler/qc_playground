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
import sys

def main():
    original_stdout = sys.stdout
    print("Calculation of the first three moment expanded states")
    print()
    print("initalize system and calculating one and two electron matrices...")
    with open('detailed_output_moment_expansion.txt', 'w') as f:
        sys.stdout = f
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
        diff_dens_list = []
        #for i in range([0,1,2,3]):
        sys.stdout = original_stdout
        print("performing single point calculation with and without external potentials")
        sys.stdout = f
        for i in [0,1,2,3]:
            if i == 0:
                ext = None
            else:
                ext = 0.0001* dip_mat[i-1, :,:]
        
            print()
            print()
            print()
            print()
            print()
            print("#######################################################", file = f)
            print("start of new single point calculation for the " + str(i)+"-th moment expanded state",file =f )
            print("#######################################################")
            print()
             
            h =  kin_mat + nuc_mat
            en, p, coeff, energy_orb = scf.scf(four_mat, noe, noo, x, h, ext)
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
            print("total dipole moment:", nuc_pole - dipole)
            print("#######################################################")
            
            str1 = "coefficients" + str(i)
            #np.savetxt("coefficients", coeff[:int(noe/2),:])
            np.savetxt(str1, coeff[:int(noe/2),:])
            str4 = "density_mat" + str(i)
            #np.savetxt("coefficients", coeff[:int(noe/2),:])
            np.savetxt(str4, p)
            
            print()
            print("#######################################################")
            #print("coefficients matrix of  " + str(i)+"th single point")
            #print(coeff[:int(noe/2),:])
            print("density matrix obtained by  " + str(i)+"-th single point calculation")
            print(p)
            print("#######################################################")
            print()


            if i == 0:
                equi_coeff = np.copy(coeff[:int(noe/2),:])
                equi_p = np.copy(p)
            else:     
                str2 = "coefficients_mes" + str(i)
                str3 = "denstiy_mes" + str(i)
                diff = coeff[:int(noe/2),:] - equi_coeff
                diff_p = p - equi_p
                np.savetxt(str2, diff)
                np.savetxt(str3, diff_p)
                diff_dens_list.append(diff_p)
                print()
                #print("#######################################################")
                #print("coefficients matrix of equilibrium single point")
                #print(equi_coeff)
                #print("#######################################################")
                #print("coefficients matrix of" + str(i)+"th moment expanded state")
                #print(diff)
                #print("#######################################################")
                print("#######################################################")
                print("density matrix from obtained from single point calculation without applied external potentials")
                print(equi_p)
                print("#######################################################")
                #print("density matrix of" + str(i)+"-th moment expanded state")
                print("density matrix of electron density difference obtained by  single point calculations with the "+ str(i)+"-th  and without an external potential")
                print(diff_p)
                print("#######################################################")
                
            print()
            print("#######################################################")
            print("end of single point calculation for the " + str(i)+"-th moment expanded state")
            print("#######################################################")
            print()
            print("#########################################################################################################################")
            print()
    sys.stdout = original_stdout
    print("post processing of response densities obtained by the external potential")
    diff_dens_list = np.array(diff_dens_list)
    diff_dens_list2 = np.copy(diff_dens_list)
    #print(diff_dens_list.shape, dip_mat.shape)
    #print(
    dip_mat2 = np.copy(dip_mat)
    tilde_mom = diff_dens_list[np.newaxis, :,:,:] * dip_mat[:, np.newaxis ,:,:]   
    tilde_mom = tilde_mom.sum(axis = (2,3))
    print()
    print("#######################################################")
    print("moment matrix of response densities obtained by the electron density difference of  single point calculation with and without an external potential ")
    print(tilde_mom*10000)
    print("#######################################################")
    tilde_mom = diff_dens_list2[np.newaxis, :,:,:] * dip_mat2[:, np.newaxis ,:,:]   
    tilde_mom = tilde_mom.sum(axis = (2,3))
    tilde_mom = (tilde_mom + tilde_mom.T)/2
    chol_left = np.linalg.cholesky(tilde_mom * -1)
    #chol_right = chol_left
    chol_right = chol_left.T
    inv_chol_right = np.linalg.inv(chol_right)
    print("#######################################################")
    print("MES moment  matrix from choleky decompostion:")
    print("#######################################################")
    print(chol_right*10000)
    #print()
    #print(inv_chol_right)
    mom_states = np.zeros(diff_dens_list2.shape)
    #for k in range(3)
    for i in range(3):
        for j in range(3):
            #mom_states[i] = inv_chol_right[i,j] * diff_dens_list[j]
            mom_states[i] += inv_chol_right[j,i] * diff_dens_list2[j]
    print("MES moment matrix obtaned by integrating over MES")
    #mom_mat = mom_states[np.newaxis,:,:,:] * dip_mat2[:, np.newaxis ,:,:]
    mom_mat =  dip_mat2[np.newaxis,: ,:,:] * mom_states[:,np.newaxis,:,:] 
    mom_mat = mom_mat.sum(axis = (2,3))
    print(mom_mat*10000)

    #tilde_mom = diff_dens_list2[np.newaxis, :,:,:] * dip_mat2[:, np.newaxis ,:,:]   
    #tilde_mom = tilde_mom.sum(axis = (2,3))
    #print("HIST")
    #print(tilde_mom*10000)
    print()
    print("saving density matrix of moment expanded states ... ")
    for j in range(mom_states.shape[0]):
        fn_mes = "mes_" + str(j+1)
        np.savetxt(fn_mes, mom_states[j])
        print(fn_mes)
