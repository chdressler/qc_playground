# import basisset
from bombanitio import create_basis
import numpy as np
import pytest
from bombanitio import read_xyz

# import basis_set_exchange as bse
# bse.get_basis("sto-3g", )["elements"]["12"]["electron_shells"][1]["coefficients"][1]
# ['0.1559162750E+00', '0.6076837186E+00', '0.3919573931E+00']
# bse.get_basis("sto-3g", )["elements"]["12"]["electron_shells"][1]["exponents"]
# ['0.1512182352E+02', '0.3513986579E+01', '0.1142857498E+01']


def test_orbitals():
    # assert create_basis.orbitals(9,12) == np.array([0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00]) ,np.array([0.2070156070E+03, 0.3770815124E+02, 0.1020529731E+02])
    coef, exp = create_basis.orbitals(4, 12)
    # assert coef == np.array([0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00])
    assert coef.astype("float64") == pytest.approx(
        np.array([0.1559162750e00, 0.6076837186e00, 0.3919573931e00])
    )
    assert exp.astype("float64") == pytest.approx(
        np.array([0.1512182352e02, 0.3513986579e01, 0.1142857498e01])
    )
    # assert coef == pytest.approx(np.array([0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00]))
    # assert create_basis.orbitals(4,12) == np.array([0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00]) ,np.array([0.1512182352E+02, 0.3513986579E+01, 0.1142857498E+01])
    # assert create_basis.orbitals(4,12) == np.array([0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00]) ,np.array([0.1512182352E+02, 0.3513986579E+01, 0.1142857498E+01])


#    1         0.2070156070E+03       0.1543289673E+00
# 2         0.3770815124E+02       0.5353281423E+00
# 3         0.1020529731E+02       0.4446345422E+00


# ORBITAL NUMBER     ORBITAL CONTRACTION EXPONENTS
# 1                  6.44364  23.80886 130.70929
# 2                  0.38039   1.16959   5.03315
# 3                  0.38039   1.16959   5.03315
# 4                  0.38039   1.16959   5.03315
# 5                  0.38039   1.16959   5.03315
# 6                  0.16886   0.62391   3.42525
# 7                  0.16886   0.62391   3.42525


def test_create_basis():
    coord, atom = read_xyz.easy_read("test_coord.xyz", None, False, False)
    coord = coord[0, :, :]
    zoa = read_xyz.atom_to_zoa(atom)
    noa = atom.shape[0]
    # basisset(noa, zoa, coord)
    noo, orb_coord, orb_type, coef_mat, exp_mat = create_basis.basisset(noa, zoa, coord)

    ref_exp = np.array(
        [
            [6.44364, 23.80886, 130.70929],
            [0.38039, 1.16959, 5.03315],
            [0.38039, 1.16959, 5.03315],
            [0.38039, 1.16959, 5.03315],
            [0.38039, 1.16959, 5.03315],
            [0.16886, 0.62391, 3.42525],
            [0.16886, 0.62391, 3.42525],
        ]
    )
    # print(ref_exp[:,::-1])
    # assert ref_exp[:,::-1].astype('float64') == pytest.approx(coef_mat)
    assert ref_exp[:, ::-1] == pytest.approx(exp_mat.astype("float64"), rel=0.0001)
    # assert ref_exp[:,::-1] == pytest.approx(coef_mat.astype('float64'), rel = 0.01)


# WATER
# 10
# 3
# 0.00000    0.00000     0.0000 8
# 0.75690    0.58585     0.0000 1
# -0.75690    0.58585     0.0000 1


def test_easy_read():
    ref_coord = np.array(
        [
            [0.00000, 0.00000, 0.0000],
            [0.75690, 0.58585, 0.0000],
            [-0.75690, 0.58585, 0.0000],
        ]
    )
    ref_atom = np.array(["O", "H", "H"])
    ref_z = np.array([8, 1, 1])
    coord, atom = read_xyz.easy_read("test_coord.xyz", None, False, False)
    coord = coord[0, :, :]
    zoa = read_xyz.atom_to_zoa(atom)
    assert coord == pytest.approx(ref_coord)
    # print(ref_atom == atom)
    # assert
    assert (atom == ref_atom).all()
    # assert len(actual) == len(expected)
    # assert all([a == b for a, b in zip(actual, expected)])
    # assert atom.all()  == ref_atom.all()
    assert (ref_z == zoa).all()
    # assert coord == pytest.approx(ref_coord)
    bohr_in_a = 0.529177210903
    coord = read_xyz.xyz_to_bohr(coord)
    ref_coord /= bohr_in_a
    assert coord == pytest.approx(ref_coord)
