============
Overview
============

nsmod is an attempt to bring a messy collection of code used in my PhD studies
into a usable format, the primary aim being to remind me what does what. I am
studying timing noise in neutron stars, `nsmod` contains code to:

Model a neutron star
====================

We provide `nsmod.Model` which collects together models including

- A simple elastic ellipsoid acted on by an electromagnetic torque

- A two component model which contains a spherical core that couples to the
  outer ellisoid (acted on by the EM torque)

There will output files saved in the `hdf5py` format, this removes the need
to define our own file structure.

Plotting
--------

The results of these models are often plotted in a similar manner, as such a
collection of functions is available in the `nsmod.Plot` module

Nonlinear dynamics
====================

In the progress of the work an effort was made to understand the results from
the models mentioned above in the context of chaos and attractors (see Duncan
& Lorimer 2010). We provide a module `nsmod.NLD_Functions` including the
associated functions




