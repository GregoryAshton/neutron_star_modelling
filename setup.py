from distutils.core import setup
from Cython.Distutils import Extension
from Cython.Distutils import build_ext
import cython_gsl


setup(
    include_dirs=[cython_gsl.get_include()],
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("nsmod_cython",
                ["src/nsmod_cython.pyx"],
                libraries=cython_gsl.get_libraries(),
                library_dirs=[cython_gsl.get_library_dir()],
                include_dirs=[cython_gsl.get_cython_include_dir()])]
    )

setup(
    include_dirs=[cython_gsl.get_include()],
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("nsmod_two_component_model",
                 ["src/nsmod_two_component_model.pyx"],
                 libraries=cython_gsl.get_libraries(),
                 library_dirs=[cython_gsl.get_library_dir()],
                 include_dirs=[cython_gsl.get_cython_include_dir()])]
    )

setup(
    include_dirs=[cython_gsl.get_include()],
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("nsmod_two_component_model2",
                 ["src/nsmod_two_component_model2.pyx"],
                 libraries=cython_gsl.get_libraries(),
                 library_dirs=[cython_gsl.get_library_dir()],
                 include_dirs=[cython_gsl.get_cython_include_dir()])]
    )



# python setup.py build_ext --inplace

