## Dependents
# from numpy.distutils.core import setup as setup_numpy
# from distutils.extension import Extension
# from Cython.Build import cythonize
# from Cython.Distutils import build_ext
# from numpy.distutils.core import Extension

from distutils.core import setup

setup(name='MK',
      version='0.0',
      description='Macromechanical model-based Marciniak-Kucynski FLD',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',
      packages=['YF','MK','HD'])
