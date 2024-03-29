# Greg Ashton: Neutron star modelling

This is a repository to collect and organise some of the code and writeup from
my PhD reseach on *timing noise in neutron stars*.

This project aims to provide a flexible set of tools to model neutron stars
and appropriate methods to plot the results. This is done through a set of
modules accesed via `nsmod.py`. The modules should be accessible from the
command line (see `nsmod.py -h` provided the directory has been added to
`PATH`), and from a python interpretor such as `IPython`. 

Dependencies 
------------

* Python 2.7 
* Python modules : numpy, scipy, matplotlib, h5py
* Cython 
* CythonGSL
* hdf5 
* C compiler such as `gcc`

All modules should be installable via `pip` or a package manager.





