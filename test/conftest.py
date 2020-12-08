import pkg_resources

import pytest

from bombanitio import read_xyz


@pytest.fixture
def coord_atom():
    filename = pkg_resources.resource_filename("bombanitio", "data/test_coord.xyz")
    return read_xyz.easy_read(filename, None, False, False)


@pytest.fixture
def coord(coord_atom):
    coord, _ = coord_atom
    coord = coord[0, :, :]
    return coord


@pytest.fixture
def atom(coord_atom):
    _, atom = coord_atom
    return atom


@pytest.fixture
def zoa(atom):
    return read_xyz.atom_to_zoa(atom)