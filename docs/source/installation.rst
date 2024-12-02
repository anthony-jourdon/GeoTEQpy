Installation
============

Downloading the package
-----------------------

To download the package, you can go to the `GitHub repository <https://github.com/anthony-jourdon/GeoTEQpy>`_ 
of the project and download the sources as a zip file or clone the repository using git:

.. code-block:: bash

    git clone https://github.com/anthony-jourdon/GeoTEQpy.git

or 

.. code-block:: bash

    git clone git@github.com:anthony-jourdon/GeoTEQpy.git

Installation
------------

1. Installation using pip (recommended)
.......................................

To install the package, it is recommended to use ``pip`` and the ``setup.py`` file provided in the sources.
To do so, you can run the following command in the root directory of the project:

.. code-block:: bash

    pip install .

or from anywhere else:

.. code-block:: bash

    pip install /path/to/GeoTEQpy/

This method will automatically download and install the required dependencies. 

2. Installation from sources (not recommended)
..............................................

If for any reason you do not want to install the package using ``pip``, 
you can install the package manually from sources.
First, you need the required dependencies:

- `numpy <https://numpy.org/>`_
- `pyvista <https://docs.pyvista.org/>`_
- `netCDF4 <https://unidata.github.io/netcdf4-python/>`_
- `ptatin3d-pyviztools <https://bitbucket.org/ptatin/ptatin3d-pyviztools>`_ branch ``anthony_jourdon/post-proc-script``

Once the dependencies are installed, move to 

.. code-block:: bash

    cd /path/to/GeoTEQpy/geoteqpy/c/

and run:

.. code-block:: bash

    make all

to compile the ``C`` code. Then, you need to add the root directory of the project to your ``PYTHONPATH``:

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:/path/to/GeoTEQpy/

or add the following line to your ``.bashrc`` or ``.bash_profile``:

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:/path/to/GeoTEQpy/

Testing the installation
------------------------

Once the package is installed, you can test the installation in a python environment by running:

.. code-block:: python

    import geoteqpy

You should see the following message:

.. code-block:: bash

    [load_clib] Successfully loaded shared object: faulttools.so
    [load_clib] Successfully loaded shared object: viztools.so

meaning that the ``C`` libraries have been correctly compiled and loaded.