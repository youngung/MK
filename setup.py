## Dependents
from distutils.core import setup
from numpy.distutils.core import setup as setup_numpy
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from numpy.distutils.core import Extension

from distutils.core import setup

# setup(name='MK',
#       version='0.0',
#       description='FLD using Marciniak-Kuczynski and macro-mechanical constitutive model',
#       author='Youngung Jeong',
#       author_email='youngung.jeong@gmail.com',
#       packages=['YF','MK','HD'])



## Fortran subroutines with f2py
ext_modules = []
ext_modules += [
    Extension(
        name="mk_for",
        sources=['src/gauss.f'],
        extra_compile_args=['-O3']
    )]

setup_numpy(ext_modules=ext_modules)
