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


def test_form_g_mat():
    np.set_printoptions(precision=3, suppress=True)
    # dens_xxx_01  fock_xxx_01 coord_xxx_01.xyz
    orca_dens = np.loadtxt(
        pkg_resources.resource_filename("bombanitio", "data/dens_xxx_03")
    )
    orca_fock = np.loadtxt(
        pkg_resources.resource_filename("bombanitio", "data/fock_xxx_03")
    )
    orca_h = np.loadtxt(pkg_resources.resource_filename("bombanitio", "data/h_xxx_03"))
    orca_fock2 = np.zeros(orca_fock.shape)
    orca_dens2 = np.zeros(orca_dens.shape)
    orca_h2 = np.zeros(orca_h.shape)

    dict1 = {}
    dict1[0] = 0
    dict1[1] = 1
    dict1[2] = 4
    dict1[3] = 2
    dict1[4] = 3
    dict1[5] = 5
    dict1[6] = 6
    dict1[7] = 9
    dict1[8] = 7
    dict1[9] = 8

    for i in range(orca_fock2.shape[0]):
        for j in range(orca_fock2.shape[1]):
            # k = i
            # l = j
            # if i  =
            orca_fock2[dict1[i], dict1[j]] = orca_fock[i, j]
            orca_dens2[dict1[i], dict1[j]] = orca_dens[i, j]
            orca_h2[dict1[i], dict1[j]] = orca_h[i, j]
    p = orca_dens2
    # p =np.array([[2.10625, -0.44603, 0.00000, 0.10859, 0.00000, -0.02843, -0.02843],
    # [-0.44603, 1.96730, 0.00000,-0.61771, 0.00000, -0.03406, -0.03406],
    # [0.00000,  0.00000, 0.73554, 0.00000, 0.00000,  0.53974, -0.53974],
    # [0.10859, -0.61771, 0.00000, 1.24023, 0.00000,  0.47272,  0.47272],
    # [0.00000,  0.00000, 0.00000, 0.00000, 2.00000,  0.00000,  0.00000],
    # [-0.02843,-0.03406, 0.53974, 0.47272, 0.00000,  0.60093, -0.19120],
    # [-0.02843,-0.03406,-0.53974, 0.47272, 0.00000, -0.19120,  0.60093]])
    coord, atom = read_xyz.easy_read(
        pkg_resources.resource_filename("bombanitio", "data/coord_xxx_03.xyz"),
        None,
        False,
        False,
    )
    coord = coord[0, :, :]
    coord = read_xyz.xyz_to_bohr(coord)
    zoa = read_xyz.atom_to_zoa(atom)
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
    np.set_printoptions(precision=5, suppress=True)
    kin_mat = eval_kin.get_kin_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)
    nuc_mat = eval_nuc.get_nuc_mat(
        noo, orb_coord, orb_type, coef_mat, exp_mat, coord, zoa
    )

    print("pre_noo:", noo)
    h = kin_mat + nuc_mat
    # print( h - orca_h2)
    # print()
    # print( (h - orca_h2)/h)
    # print()
    # print(orca_h2)
    # print()
    # print(h)
    assert h == pytest.approx(orca_h2, rel=1e-3, abs=1e-4)

    four_mat = eval_four.get_four_mat(noo, orb_coord, orb_type, coef_mat, exp_mat)

    g = scf.formg(noo, four_mat, p)
    f = h + g
    # fock_ref = np.loadtxt("../data/fock.data")
    # print()
    # print(f - fock_ref)
    # print()
    # print((f - fock_ref)/fock_ref)
    # print()
    # print(f )
    # print()
    # print(fock_ref)

    # ../data/orca_dens_eq ../data/orca_fock_eq
    # assert f == pytest.approx(fock_ref, rel=1e-3, abs=1e-5)

    print(f - orca_fock2)
    print()
    # print((f - orca_fock2)/f)
    print()
    print(f)
    print()
    print(orca_fock2)

    assert f == pytest.approx(orca_fock2, rel=1e-3, abs=1e-4)
    assert p == pytest.approx(orca_dens2, rel=1e-3, abs=1e-4)
