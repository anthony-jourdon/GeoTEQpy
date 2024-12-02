# GeoTEQpy

`geoteqpy` is a Python package designed to link long-term geodynamic models and dynamic rupture models.

Currently, it can perform 3 main tasks:

- Structured mesh extrusion,
- Computation of the medial axis of a 3D arbitrary shape,
- Export of [pTatin3d](https://github.com/laetitialp/ptatin-gene) data to [ASAGI](https://github.com/TUM-I5/ASAGI) using [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) format.

## Installation
Once downloaded using:

> git clone git@github.com:anthony-jourdon/GeoTEQpy.git

or 

> git clone https://github.com/anthony-jourdon/GeoTEQpy.git

`geoteqpy` can be installed using 2 different methods.
### Method 1: Using `pip` (recommended)
```sh
pip install /path/to/geoteqpy
```
This will automatically download and install dependencies.
### Method 2 (not recommended)
First, all the dependencies need to be installed:

- [numpy](https://numpy.org/)
- [pyvista](https://docs.pyvista.org/)
- [netCDF4](https://unidata.github.io/netcdf4-python/)
- [ptatin3d-pyviztools](https://bitbucket.org/ptatin/ptatin3d-pyviztools) branch `anthony_jourdon/post-proc-script`

Once in the `geoteqpy` directory, move to

> geoteqpy/c

and type
```sh
make all
```
This should build the `faulttools.so` shared library.
Once done, you need to add the location of `geoteqpy` to your `PYTHONPATH` environment variable. This can be done with

```sh
export PYTHONPATH=$PYTHONPATH:$PWD
```

if the command is run inside the repository containing the module `geoteqpy`.

### Test installation
Once installed, open a python shell and type
```python
import geoteqpy
```
if there are no errors everything is good!
You should also see the following message indicating that the shared libraries has been succesfully loaded:

```
[load_clib] Successfully loaded shared object: faulttools.so
[load_clib] Successfully loaded shared object: viztools.so
```

## Documentation
An HTML documentation is available.
The documentation can be built and accessed locally (no internet connection required) from the repository using
```sh
sphinx-build -M html docs/source/ docs/build/
```
However, this requires the Python package [sphinx](https://www.sphinx-doc.org/en/master/index.html).
An online documentation is available at this [link](https://geoteqpy.readthedocs.io/en/latest/).