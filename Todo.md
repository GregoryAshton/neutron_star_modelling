# Todo list for nsmod

## `src/ pyx` files

1. Merge `nsmod_two_component2` and `nsmod_two_component` so that initial data is
specified by the user and not the function.

2. Update `lib/Model.py` to use nsmod_one_component model rather than `nsmod_cython`
this should be removed

3. Add anom torque equations to `nsmod_one_component_model`

## `lib/Plot.py`


## `lib/File_Functions.py`

1. Currently the computed time scales such as $\tau_{a}$ only make sense in a
Biaxial case.

## `lib/NLD_Functions.py`

1. URGENT: Dotted variables works only the biaxial case, this needs updating
 before we can progress
