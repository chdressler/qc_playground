import setuptools
from numpy.distutils.core import setup, Extension


lib = Extension(
    name="bombanitio.special.special.f", sources=["bombanitio/special/special.f"]
)

setup(
    name="bombanitio",
    packages=["bombanitio"],
    ext_modules=[lib],
    install_requires=["basis_set_exchange", "scipy"],
    entry_points={
        "console_scripts": [
            "write_to_cube = bombanitio.dens_to_cube:dens_to_cube",
            "check_numeric = bombanitio.dens_to_cube:calc_moments_numerically",
            "bombanitio = bombanitio.main_bombanitio:main",
            "moment_expansion = bombanitio.mom_exp:main",
        ],
    },
)
