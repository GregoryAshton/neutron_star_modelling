#!/usr/bin/env python 

from distutils.core import setup
from Cython.Distutils import Extension
from Cython.Distutils import build_ext
import cython_gsl
import os

os.system("mv nsmod/*.so .")

name_list = ["nsmod_cython", 
             "nsmod_one_component_model",
             "nsmod_two_component_model",
             "nsmod_two_component_model2",
             "nsmod_one_component_model_with_Euler",
             "switching_torque_with_Euler"
             ]
extension_list = [Extension(name,
                ["src/"+name+".pyx"],
                libraries=cython_gsl.get_libraries(),
                library_dirs=[cython_gsl.get_library_dir()],
                include_dirs=[cython_gsl.get_cython_include_dir()])
                for name in name_list]

setup(
    include_dirs=[cython_gsl.get_include()],
    cmdclass={'build_ext': build_ext},
    ext_modules=extension_list
    )

os.system("mv *.so nsmod/")
