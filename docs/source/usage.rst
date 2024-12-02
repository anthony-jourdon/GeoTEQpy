Usage
=====

GeoTEQpy can perform 3 main tasks:

1. Structured mesh extrusion in any direction
2. Computation of the medial axis of an arbitrary 3D shape
3. Export data from `pTatin3d`_ to `SeisSol`_ using the format required by `ASAGI`_

Mesh extrusion
--------------
Mesh extrusion can be performed using the executable ``scripts/mesh_extrude.py``.
This executable reads a YAML file containing the parameters needed for the extrusion.
Use the following command to run the executable:

.. code-block:: bash

    python scripts/mesh_extrude.py -f path/to/parameters.yaml

Use the following command to see the available options:

.. code-block:: bash

    python scripts/mesh_extrude.py -h

Medial axis computation
-----------------------
The medial axis of a 3D shape can be computed using the executable ``scripts/get_medial_axis.py``.
This executable reads a YAML file containing the parameters needed for the computation.
Use the following command to run the executable:

.. code-block:: bash

    python scripts/get_medial_axis.py -f path/to/parameters.yaml

Use the following command to see the available options:

.. code-block:: bash

    python scripts/get_medial_axis.py -h

Export data from `pTatin3d`_ to `SeisSol`_
------------------------------------------
Exporting data from `pTatin3d`_ to `SeisSol`_ can be performed using the executable ``scripts/ptatin2asagi.py``.
This executable reads a YAML file containing the parameters needed for the export.
Use the following command to run the executable:

.. code-block:: bash

    python scripts/ptatin2asagi.py -f path/to/parameters.yaml

Use the following command to see the available options:

.. code-block:: bash

    python scripts/ptatin2asagi.py -h

Additional executable
---------------------
In addition to the previous functionalities, GeoTEQpy provides an executable to rotate faults surface from a coordinate
system where :math:`y` is the vertical direction (e.g. `pTatin3d`_) to a 
coordinate system where :math:`z` is the vertical direction (e.g. `SeisSol`_).
This executable reads a YAML file containing the parameters needed for the export.
Use the following command to run the executable:

.. code-block:: bash

    python scripts/fault_surface_rotation.py -f path/to/parameters.yaml

Use the following command to see the available options:

.. code-block:: bash

    python scripts/fault_surface_rotation.py -h