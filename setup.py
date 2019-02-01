from distutils.core import setup, Extension
import numpy as np

sourcefiles = ['PolSolver.cpp']

c_ext = Extension("PolSolver", sources=sourcefiles,
                  extra_compile_args=["-Wno-deprecated","-O3"])
             #     libraries=['cfitsio','gsl','fftw3'],
             #     extra_link_args=["-Xlinker", "-export-dynamic"])

setup(ext_modules=[c_ext], include_dirs=[np.get_include()])






