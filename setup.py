## Dependents
from distutils.core import setup
from numpy.distutils.core import setup as setup_numpy
from distutils.extension import Extension
# from Cython.Build import cythonize
# from Cython.Distutils import build_ext
from numpy.distutils.core import Extension

from distutils.core import setup
import os, shutil, site

setup(name='MK',
      version='0.0',
      description='FLD using Marciniak-Kuczynski and macro-mechanical constitutive model',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',
      url='https://github.com/youngung/mk.git',

      packages=['mk',
                'mk.yieldFunction',
                'mk.library',
                'mk.materials',
                'mk.tests'],
      package_dir={'mk':'mk',
                   'mk.yieldFunction':'mk/yieldFunction',
                   'mk.library':'mk/library',
                   'mk.materials':'mk/materials',
                   'mk.tests':'mk/tests'
      })


## Compile the stand-alone binary of yld2000
path_site = site.getsitepackages()[0]
os.system('bash setup_for.sh')
shutil.copy('fortran/yld2000_sa',path_site)
# os.remove('fortran/yld2000_sa')


## Fortran subroutines with f2py
ext_modules = []
ext_modules += [
        Extension(
            name="yf_for",
            sources=['fortran/for.f'],
            # extra_compile_args=['-O3']
        )]
ext_modules += [
        Extension(
            name="yf_yld2000",
            sources=['fortran/yld2000.f'],
            depends=['fortran/yld2000_lib.f'],
            # extra_compile_args=['-O3']
        )]
setup_numpy(ext_modules=ext_modules)
