import os

version = [1, 0, 0]

from .register_cfunc import *
clib, cfunc = load_clib()

from .ptatin_load       import *
from .asagi             import *
from .rotations         import *
from .ptatin_to_seissol import *
from .field_utils       import *
from .mesh_extrusion    import *
from .medial_axis       import *