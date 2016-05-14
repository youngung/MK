from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("func_fld_cy.pyx")
)
