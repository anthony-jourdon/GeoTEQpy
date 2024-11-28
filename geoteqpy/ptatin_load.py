import pyvista as pvs
import numpy as np
import pyviztools as pvt

class ModelData:
  """
  Class to read and manipulate data from a pTatin3d simulation.
  It also provides methods to interpolate data from the pTatin3d mesh to a uniform mesh.

  :param str vts:
    Name of the vts file to read. Default: ``None``.

  :Attributes:

  .. py:attribute:: ptatin_mesh
    :type: pvt.cfemesh.CFEMeshQ1

    pTatin3d mesh object. Default: ``None``.

  .. py:attribute:: uniform_mesh
    :type: pvt.cfemesh.CFEMeshQ1

    Uniform mesh object. Default: ``None``.

  .. py:attribute:: O
    :type: np.ndarray

    Minimum coordinates of the domain in each direction. Default: ``None``.

  .. py:attribute:: L
    :type: np.ndarray

    Maximum coordinates of the domain in each direction. Default: ``None``.
  
  :Methods:
  """
  def __init__(self,vts:str|None=None) -> None:
    self.ptatin_mesh:pvt.cfemesh.CFEMeshQ1|None  = None
    self.uniform_mesh:pvt.cfemesh.CFEMeshQ1|None = None
    self.O:np.ndarray|None=None
    self.L:np.ndarray|None=None

    if vts is not None:
      print(f"Creating ptatin_mesh from file {vts}")
      self.ptatin_mesh = self.mesh_create_from_vts(vts)

    return

  def check_ptatin_mesh(self) -> None:
    """
    Verify that a mesh object has been created for the pTatin3d simulation results.

    :raises ValueError: If the mesh object does not exist.
    """
    if self.ptatin_mesh is None:
      serr = f'ptatin_mesh has not been created yet.\n'
      serr += f'Use {self.mesh_create_from_vts.__name__}() or {self.mesh_create_from_options.__name__}() first.'
      raise ValueError(serr)
    return

  def set_default_fields(self) -> list[str]:
    """
    Set default fields to interpolate from the pTatin3d mesh to the uniform mesh.
    By default, all cell and point fields are considered.

    :return: List of keys corresponding to fields to interpolate.
    :rtype: list[str]
    """
    fields = list(self.ptatin_mesh.cell_data.keys())
    fields.extend(list(self.ptatin_mesh.point_data.keys()))
    return fields

  def mesh_create_from_vts(self, fname:str) -> pvt.cfemesh.CFEMeshQ1:
    """
    Create a :py:class:`pyviztools.CFEMeshQ1` object from a vts file.

    :param str fname: Name of the vts file to read.

    :return: Mesh object.
    :rtype: :py:class:`pyviztools.CFEMeshQ1`
    """
    pv_mesh:pvs.StructuredGrid = pvs.read(fname)
    # get bounding box
    self.O  = np.array([ pv_mesh.bounds[0], pv_mesh.bounds[2], pv_mesh.bounds[4] ], dtype=np.float32)
    self.L  = np.array([ pv_mesh.bounds[1], pv_mesh.bounds[3], pv_mesh.bounds[5] ], dtype=np.float32)
    size    = np.asarray(pv_mesh.dimensions)
    
    # create finite element mesh object
    mesh = pvt.cfemesh.CFEMeshQ1(3,size[0]-1,size[1]-1,size[2]-1)
    # attach coords
    mesh.coor = pv_mesh.points
    # create element to vertex connectivity
    mesh.create_e2v()
    # copy fields
    mesh.point_data = {}
    mesh.cell_data  = {}
    for field in pv_mesh.point_data:
      if mesh.is_vertex_field(pv_mesh.point_data[field]):
        mesh.point_data[field] = pv_mesh.point_data[field]
    for field in pv_mesh.cell_data:
      if mesh.is_cell_field(pv_mesh.cell_data[field]):
        mesh.cell_data[field] = pv_mesh.cell_data[field]
    return mesh
  
  def mesh_create_from_options(self, O:np.ndarray, L:np.ndarray, m:np.ndarray) -> pvt.cfemesh.CFEMeshQ1:
    """
    Create a :py:class:`pyviztools.CFEMeshQ1` object from a bounding box and mesh size.

    :param np.ndarray O: Minimum coordinates of the domain in each direction.
    :param np.ndarray L: Maximum coordinates of the domain in each direction.
    :param np.ndarray m: Number of elements in each direction.

    :return: Mesh object.
    :rtype: :py:class:`pyviztools.CFEMeshQ1`

    :Example:

    .. code:: python

      import numpy as np
      
      O = np.array([0,0,0], dtype=np.float32)
      L = np.array([1,1,1], dtype=np.float32)
      m = np.array([8,8,8], dtype=np.int32)

      M = ModelData()
      mesh = M.mesh_create_from_options(O,L,m)
    """
    # create mesh object
    mesh = pvt.cfemesh.CFEMeshQ1(3,m[0],m[1],m[2])
    # create coords
    pvt.utils.femesh_set_uniform_coordinates(mesh,np.float32(O),np.float32(L))
    # create element to vertex connectivity
    mesh.create_e2v()
    mesh.point_data = {}
    mesh.cell_data  = {}
    return mesh

  def mesh_create_from_petscvec(self):
    raise NotImplementedError(f'Method {self.mesh_create_from_petscvec.__name__} not implemented yet')

  def get_strainrate(self,**keys) -> np.ndarray:
    """
    Get strain rate from its key in the cell data of the pTatin3d mesh or compute it from the velocity field such that:

    .. math::

      \\dot{\\varepsilon}_{ij} = \\frac{1}{2} \\left( \\frac{\\partial u_i}{\\partial x_j} + \\frac{\\partial u_j}{\\partial x_i} \\right)

    where :math:`\\dot{\\varepsilon}_{ij}` is the strain rate tensor and :math:`u_i` is the velocity vector.
    
    :param dict keys: Dictionary of keys to identify the strain rate field and the velocity field.
      
      - ``strainrate_key`` (str): Key to identify the strain rate field in the cell data.
      - ``velocity_key`` (str): Key to identify the velocity field in the point data.
    
    The velocity field is mandatory only if the strain rate field is not found in the cell data.
    If the strain rate key is not provided or not found in the data it will be computed from the velocity field.

    :raises ValueError: If both the strain rate and the velocity field are not found or provided.

    :return: Strain rate field. Shape ``(n_cells,9)``.
    :rtype: np.ndarray
    """
    self.check_ptatin_mesh()
    strainrate_key:str = keys.get('strainrate_key',None)
    velocity_key:str   = keys.get('velocity_key',None)

    if strainrate_key is None or strainrate_key not in self.ptatin_mesh.cell_data:
      print(f"strainrate_key: {strainrate_key} not found in {self.ptatin_mesh.cell_data.keys()}. Computing it from velocity field identified as: '{velocity_key}'")
      if velocity_key is None or velocity_key not in self.ptatin_mesh.point_data:
        raise ValueError(f'velocity_key: {velocity_key} not found in {self.ptatin_mesh.point_data.keys()}')
      strain_rate = pvt.utils.compute_strainrate_tensor(self.ptatin_mesh,self.ptatin_mesh.point_data[velocity_key])
    else:
      strain_rate = self.ptatin_mesh.cell_data[strainrate_key]
    return strain_rate

  def get_deviatoric_stress(self,**keys):
    """
    Get deviatoric stress from its key in the cell data of the pTatin3d mesh 
    or compute it from the viscosity and strain rate fields such that:

    .. math::

      \\tau_{ij} = 2 \\eta \\dot{\\varepsilon}_{ij}
    
    where :math:`\\tau_{ij}` is the deviatoric stress tensor, :math:`\\eta` is the viscosity and :math:`\\dot{\\varepsilon}_{ij}` is the strain rate tensor.

    :param dict keys: Dictionary of keys: 
      
      - ``deviatoric_stress_key`` (str): Key to identify the deviatoric stress field in the cell data.
      - ``viscosity_key`` (str): Key to identify the viscosity field in the cell data.
      - ``strainrate_key`` (str): Key to identify the strain rate field in the cell data.
      - ``velocity_key`` (str): Key to identify the velocity field in the point data.
    
    If ``deviatoric_stress_key`` is provided and found in data, the deviatoric stress field is returned and no other 
    field key is required. If not found, the deviatoric stress field is computed from the viscosity and strain rate fields.
    If the strain rate field is not found, it will be computed from the velocity field using the method :py:meth:`get_strainrate <geoteqpy.ModelData.get_strainrate>`. 
    Therefore, the only two mandatory fields are the viscosity and the velocity fields.

    :return: Deviatoric stress field. Shape ``(n_cells,9)``.
    :rtype: np.ndarray
    """
    self.check_ptatin_mesh()
    dev_stress_key:str = keys.get('deviatoric_stress_key',None)
    viscosity_key:str  = keys.get('viscosity_key',None)
    strainrate_key:str = keys.get('strainrate_key',None)
    velocity_key:str   = keys.get('velocity_key',None)

    if dev_stress_key is None or dev_stress_key not in self.ptatin_mesh.cell_data:
      strain_rate = self.get_strainrate(strainrate_key=strainrate_key,velocity_key=velocity_key)
      deviatoric_stress = pvt.utils.compute_deviatoric_stress(self.ptatin_mesh.cell_data[viscosity_key],strain_rate)
      deviatoric_stress.shape = (deviatoric_stress.shape[0],9)
    else:
      deviatoric_stress = self.ptatin_mesh.cell_data[dev_stress_key]
    return deviatoric_stress

  def get_total_stress(self,**keys):
    """
    Get the total stress field from its key or compute it from the deviatoric stress and the pressure fields such that:

    .. math::

      \\sigma_{ij} = \\tau_{ij} - p \\delta_{ij}

    where :math:`\\tau_{ij}` is the deviatoric stress tensor, :math:`p` is the pressure and :math:`\\delta_{ij}` is the Kronecker delta.
    
    :param dict keys: Dictionary of keys:

      - ``total_stress_key`` (str): Key to identify the total stress field in the cell data.
      - ``deviatoric_stress_key`` (str): Key to identify the deviatoric stress field in the cell data.
      - ``viscosity_key`` (str): Key to identify the viscosity field in the cell data.
      - ``strainrate_key`` (str): Key to identify the strain rate field in the cell data.
      - ``velocity_key`` (str): Key to identify the velocity field in the point data.
      - ``pressure_key`` (str): Key to identify the pressure field in the cell data.

      
    If the total stress field is not found in the cell data, it will be computed from 
    the deviatoric stress using the method :py:meth:`get_deviatoric_stress <geoteqpy.ModelData.get_deviatoric_stress>` 
    and the pressure fields. See :py:meth:`get_deviatoric_stress <geoteqpy.ModelData.get_deviatoric_stress>` for more details.

    :return: Total stress field. Shape ``(n_cells,9)``.
    :rtype: np.ndarray
    """
    self.check_ptatin_mesh()
    total_stress_key:str = keys.get('total_stress_key',None)
    dev_stress_key:str   = keys.get('deviatoric_stress_key',None)
    viscosity_key:str    = keys.get('viscosity_key',None)
    strainrate_key:str   = keys.get('strainrate_key',None)
    velocity_key:str     = keys.get('velocity_key',None)
    pressure_key:str     = keys.get('pressure_key',None)

    if total_stress_key is None or total_stress_key not in self.ptatin_mesh.cell_data:
      deviatoric_stress = self.get_deviatoric_stress(deviatoric_stress_key=dev_stress_key,viscosity_key=viscosity_key,strainrate_key=strainrate_key,velocity_key=velocity_key)
      total_stress = self.ptatin_mesh.create_cell_field(dof=9)
      
      if pressure_key is None:
        raise ValueError(f'Providing a pressure_key is mandatory to compute the total stress field')
      if pressure_key not in self.ptatin_mesh.cell_data:
        if pressure_key in self.ptatin_mesh.point_data:
          pressure = self.ptatin_mesh.interp2cell(np.float32(self.ptatin_mesh.point_data[pressure_key]))
        else:
          raise ValueError(f'pressure_key: {pressure_key} not found in {self.ptatin_mesh.cell_data.keys()}or {self.ptatin_mesh.point_data.keys()}')
      else:
        pressure = self.ptatin_mesh.cell_data[pressure_key]

      for i in range(3):
        for j in range(3):
          if i == j:
            total_stress[:,3*i+j] = deviatoric_stress[:,i*3+j] - pressure
          else:
            total_stress[:,3*i+j] = deviatoric_stress[:,i*3+j]
    else:
      total_stress = self.ptatin_mesh.cell_data[total_stress_key]
    return total_stress

  def interpolate_to_uniform_mesh(self,fields:list|None=None) -> None:
    """
    Interpolate fields from the pTatin3d mesh created with
    :py:meth:`mesh_create_from_vts <geoteqpy.ModelData.mesh_create_from_vts>` to a uniform mesh created with
    :py:meth:`mesh_create_from_options <geoteqpy.ModelData.mesh_create_from_options>`.

    .. warning:: During the interpolation, all cell fields are projected to point fields.

    .. note:: Interpolation is linear and uses :math:`Q_1` isoparametric element shape functions.

    :param list fields: List of keys corresponding to fields to interpolate. Default: ``None``.
        If no fields are provided, all cell and point fields present on the ptatin mesh object are considered.

    :raises ValueError: If the uniform mesh has not been created yet.
    """
    self.check_ptatin_mesh()
    if self.uniform_mesh is None:
      serr = f'Error: {self.uniform_mesh} has not been created yet.\n'
      serr += f'Use {self.mesh_create_from_options.__name__} first.'
      raise ValueError(serr)
    
    # if no fields are provided, interpolate all fields
    if fields is None:
      fields = self.set_default_fields()

    # locate points in pTatin3d FE mesh and attribute them an element index and local coords
    el_p,xi_p = pvt.utils.femesh_point_location(self.ptatin_mesh, self.uniform_mesh.coor)
    # first, project all cell fields to point fields
    for field in fields:
      if field in self.ptatin_mesh.cell_data:
        if self.ptatin_mesh.is_cell_field(self.ptatin_mesh.cell_data[field]):
          self.ptatin_mesh.point_data[field] = self.ptatin_mesh.cell_field_to_point_field(self.ptatin_mesh.cell_data[field],volume_scale=True)

    # interpolate from pTatin to uniform grid
    for field in fields:
      self.uniform_mesh.point_data[field] = self.ptatin_mesh.interp_vertex_field(np.float32(self.ptatin_mesh.point_data[field]),el_p,xi_p)
    return
  
def test():
  import os

  mdir = os.path.join(os.environ['M3D_DIR'],'StrikeSlip3D/MiddleRes/EastAnatolia')
  vts  = os.path.join(mdir,'step000825-cellfields.vts')

  ptatin = ModelData()
  ptatin.ptatin_mesh  = ptatin.mesh_create_from_vts(vts)
  ptatin.uniform_mesh = ptatin.mesh_create_from_options(ptatin.O,ptatin.L,np.array([256,64,128],dtype=np.int32))

  ptatin.ptatin_mesh.cell_data['dev_stress'] = ptatin.get_deviatoric_stress(
    viscosity_key='viscosity',
    strainrate_key='StrainRate',
    velocity_key='u')

  ptatin.ptatin_mesh.cell_data['total_stress'] = ptatin.get_total_stress(
    deviatoric_stress_key='dev_stress',
    pressure_key='pressure_p')
  """
  pvt.simple_vts_writer(
    os.path.join(mdir,'step000825-total_stress.vts'),
    ptatin.ptatin_mesh,
    cellData=ptatin.ptatin_mesh.cell_data,
    pointData=ptatin.ptatin_mesh.point_data)
  """
  ptatin.interpolate_to_uniform_mesh(fields=['total_stress','density'])
  pvt.simple_vts_writer(
    os.path.join(mdir,'step000825-interpolated.vts'),
    ptatin.uniform_mesh,
    pointData=ptatin.uniform_mesh.point_data)
  return

if __name__ == '__main__':
  test()