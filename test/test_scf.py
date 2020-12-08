import pkg_resources

from bombanitio import (
    eval_over,
    eval_kin,
    eval_nuc,
    eval_four,
    eval_dipole,
    read_xyz,
    create_basis,
    scf,
)


import pytest
import numpy as np


def test_form_g_mat(coord, atom, zoa):
    p = np.array(
        [
            [2.10625, -0.44603, 0.00000, 0.10859, 0.00000, -0.02843, -0.02843],
            [-0.44603, 1.96730, 0.00000, -0.61771, 0.00000, -0.03406, -0.03406],
            [0.00000, 0.00000, 0.73554, 0.00000, 0.00000, 0.53974, -0.53974],
            [0.10859, -0.61771, 0.00000, 1.24023, 0.00000, 0.47272, 0.47272],
            [0.00000, 0.00000, 0.00000, 0.00000, 2.00000, 0.00000, 0.00000],
            [-0.02843, -0.03406, 0.53974, 0.47272, 0.00000, 0.60093, -0.19120],
            [-0.02843, -0.03406, -0.53974, 0.47272, 0.00000, -0.19120, 0.60093],
        ]
    )
    coord = read_xyz.xyz_to_bohr(coord)
    noa = atom.shape[0]
    # basisset(noa, zoa, coord)
    noo, orb_coord, orb_type, coef_mat, exp_mat = create_basis.basisset(noa, zoa, coord)
    # print(noo)
    # print(orb_coord)
    # print(orb_type)
    # print(coef_mat)
    # print(exp_mat)
    read_xyz.plot_initial_orbs(noo, orb_coord, orb_type, coef_mat, exp_mat)
    #    print(coef_mat.dtype)

    over_mat = eval_over.get_overlap_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    np.set_printoptions(precision=4, suppress=True)
    kin_mat = eval_kin.get_kin_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    nuc_mat = eval_nuc.get_nuc_mat(
        noo, orb_coord, orb_type, coef_mat, exp_mat, coord, zoa
    )

    print("pre_noo:", noo)
    four_mat = eval_four.get_four_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)

    h = kin_mat + nuc_mat
    g = scf.formg(noo, four_mat, p)
    f = h + g
    fock_ref = np.loadtxt(
        pkg_resources.resource_filename("bombanitio", "data/fock.data")
    )
    # print()
    # print(f - fock_ref)
    # print()
    # print((f - fock_ref)/fock_ref)
    # print()
    # print(f )
    # print()
    # print(fock_ref)

    # ../data/orca_dens_eq ../data/orca_fock_eq
    assert f == pytest.approx(fock_ref, rel=1e-3, abs=1e-5)
    orca_dens = np.loadtxt(
        pkg_resources.resource_filename("bombanitio", "data/orca_dens_eq")
    )
    # print(orca_dens - p)
    # print(

    orca_fock = np.loadtxt(
        pkg_resources.resource_filename("bombanitio", "data/orca_fock_eq")
    )
    orca_fock2 = np.zeros(orca_fock.shape)
    orca_dens2 = np.zeros(orca_dens.shape)

    dict1 = {}
    dict1[0] = 0
    dict1[1] = 1
    dict1[2] = 4
    dict1[3] = 2
    dict1[4] = 3
    dict1[5] = 5
    dict1[6] = 6

    for i in range(orca_fock2.shape[0]):
        for j in range(orca_fock2.shape[1]):
            # k = i
            # l = j
            # if i  =
            orca_fock2[dict1[i], dict1[j]] = orca_fock[i, j]
            orca_dens2[dict1[i], dict1[j]] = orca_dens[i, j]

    print(f - orca_fock)
    print()
    print(f)
    print()
    print(orca_fock)
    assert f == pytest.approx(orca_fock2, rel=1e-3, abs=1e-5)
    assert p == pytest.approx(orca_dens2, rel=1e-3, abs=1e-5)
