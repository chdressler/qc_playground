from bombanitio import eval_over
from bombanitio import eval_kin
from bombanitio import eval_nuc
from bombanitio import eval_four
from bombanitio import eval_dipole
from bombanitio import read_xyz
from bombanitio import create_basis
from bombanitio import scf
import pytest
import numpy as np


def test_orca_water():
    coord, atom = read_xyz.easy_read("test_coord.xyz", None, False, False)
    coord = coord[0, :, :]
    coord = read_xyz.xyz_to_bohr(coord)
    zoa = read_xyz.atom_to_zoa(atom)
    noa = atom.shape[0]
    noo, orb_coord, orb_type, coef_mat, exp_mat = create_basis.basisset(noa, zoa, coord)
    read_xyz.plot_initial_orbs(noo, orb_coord, orb_type, coef_mat, exp_mat)

    over_mat = eval_over.get_overlap_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    np.set_printoptions(precision=4, suppress=True)
    kin_mat = eval_kin.get_kin_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    nuc_mat = eval_nuc.get_nuc_mat(
        noo, orb_coord, orb_type, coef_mat, exp_mat, coord, zoa
    )

    dip_matx = eval_dipole.get_dip_mat(noo, orb_coord, orb_type, coef_mat, exp_mat, 0)
    dip_maty = eval_dipole.get_dip_mat(noo, orb_coord, orb_type, coef_mat, exp_mat, 1)
    dip_matz = eval_dipole.get_dip_mat(noo, orb_coord, orb_type, coef_mat, exp_mat, 2)
    dip_mat = np.array([dip_matx, dip_maty, dip_matz])

    four_mat = eval_four.get_four_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    h = kin_mat + nuc_mat
    ##diagonalisation of S
    # eig, u = np.linalg.eig(over_mat)
    ##creation of matrix x via canonical orthogonalisation
    # x = u / (eig) ** 0.5
    ## xt is transposed of x
    # xt = x.T
    x = scf.ortho_basis(over_mat)
    noe = zoa.sum()
    # print("number of electron:", noe)
    # en, p, c = scf.scf(four_mat, noe, noo, x, xt, h)
    en, p, c, coeff = scf.scf(four_mat, noe, noo, x, h)
    etot = scf.calc_tot_en(noa, en, coord, zoa)
    # print(etot)
    # dipole = dip_mat * p
    # dipole = dipole.sum(axis=1)
    # dipole = dipole.sum(axis=1)
    dipole = eval_dipole.get_el_dipole(p, dip_mat)
    nuc_pole = eval_dipole.get_nuc_pole(coord, zoa)
    tot_dip = nuc_pole - dipole
    orca_en = -74.962919685003
    orca_dip = np.array([0.000000, 0.678999, -0.000000])
    assert etot == pytest.approx(orca_en, rel=1e-5, abs=1e-4)
    assert tot_dip == pytest.approx(orca_dip, rel=1e-5, abs=1e-3)
