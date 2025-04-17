#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  GeoTEQpy
#  filename: medial_axis.py
#
#  This file is part of GeoTEQpy.
#
#  GeoTEQpy is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, either 
#  version 3 of the License, or any later version.
#
#  GeoTEQpy is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with GeoTEQpy. 
#  If not, see <https://www.gnu.org/licenses/>.
#====================================================================================================

import pyvista as pvs
import numpy as np
import geoteqpy as gte
import pyviztools.utils as fe_utils
from pyviztools import CFEMeshQ1
from time import perf_counter

class MedialAxis:
  """
  Class to compute the medial axis of a 3-dimensional shape.
  The medial axis is a set of points located at the centre of the shape.
  Medial axis is computed using the algorithm described in the paper:

    - Ma, J., Bae, S.W., Choi, S. (2019). 3D medial axis point approximation using nearest neighbors and the normal field. 
      Computer-Aided Design, 114, 1-11. 
      https://doi.org/10.1016/j.cad.2019.05.002

  :param pyvista.UnstructuredGrid mesh:
    Mesh of the shape for which the medial axis is computed.
  :param float radius_ma:
    Radius of research for the medial axis.
  :param float radius_cov:
    Radius of research for the covariance matrix.

  :Attributes:

  .. py:attribute:: medial_axis
    :type: pyvista.PolyData

    Pyvista mesh object of the medial axis points.
    
  .. py:attribute:: radius_ma
    :type: float

    Radius of research for the medial axis.

  .. py:attribute:: radius_cov
    :type: float

    Radius of research for the covariance matrix.

  :Example:

  .. code:: python

    import pyvista as pvs
    import geoteqpy as gte

    # Load a mesh
    mesh = pvs.read("path/to/mesh.vtu")
    # Create a MedialAxis object
    ma = gte.MedialAxis(mesh=mesh,radius_ma=0.5,radius_cov=1.0)
    # Compute the medial axis and the orientation vectors
    medial_axis = ma.get_medial_axis_mesh()

  :Methods:

  """
  def __init__(self,mesh:pvs.UnstructuredGrid,radius_ma:float,pca_method:str="sphere",cov_radius:float=1e4,cov_knearest:int=100) -> None:
    self.mesh                     = mesh
    self.radius_ma                = radius_ma
    self.medial_axis:pvs.PolyData = None
    # compute the medial axis
    self.compute_medial_axis()

    self.pca:gte.PCA = None
    if pca_method == "sphere":
      self.pca = gte.SpherePCA(points=self.medial_axis.points,radius=cov_radius)
    elif pca_method == "knearest":
      self.pca = gte.KNearestPCA(points=self.medial_axis.points,knearest=cov_knearest)
    self.get_orientation_vectors()
    return
  
  def get_medial_axis_mesh(self) -> pvs.PolyData:
    """
    Get the medial axis mesh.

    :param bool get_eigv:
      If ``True``, the orientation vectors are computed and added to the mesh.
      Default is ``True``.

    :return:
      Pyvista mesh object of the medial axis points.
    :rtype: pyvista.PolyData
    """
    if self.medial_axis is None:
      self.compute_medial_axis()
    self.get_orientation_vectors()
    return self.medial_axis

  def compute_medial_axis(self) -> pvs.PolyData:
    """
    Compute the medial axis of the shape.
    Attach the result to the class attribute :py:attr:`medial_axis <geoteqpy.MedialAxis.medial_axis>`. 
    """
    t0 = perf_counter()
    # get points coordinates
    points = self.mesh.points
    # get points normals
    normals = self.mesh.point_data['Normals']
    # compute the medial axis
    # Verify the array dimensions
    if points.ndim != 2:
      raise RuntimeError(f'points coordinates must be of the shape (npoints, 3). Given shape: {points.shape}')
    # Verify the spatial dimension
    ndim = points.shape[1] 
    if ndim != 3:
      raise RuntimeError(f'points coordinates must be of the shape (npoints, 3). Given shape: {points.shape}')
    
    # Get number of points
    npoints = np.int64(points.shape[0])
    # Verify that there are the same number of points and normals
    if normals.shape[0] != npoints:
      raise RuntimeError(f'The number of points and the number of normals are different: (npoints, nnormals) = ({npoints,normals.shape[0]})') 

    # flatten arrays for c memory
    points_v = np.float64( np.reshape(points,  (ndim*npoints)) )
    normal_v = np.float64( np.reshape(normals, (ndim*npoints)) )
    medial_axis_v = np.zeros(shape=(ndim*npoints), dtype=np.float64)
    # call the c function that computes the medial axis
    gte.cfunc["utils"]["medial_axis"](
      points_v, # flattened points coordinates
      normal_v, # flattened points normals
      npoints, # number of points
      np.float64(self.radius_ma), # radius of research for the medial axis
      np.float64(-1),
      np.float64(-1),
      medial_axis_v # returned coordinates of the points in the medial axis
    )
    # reshape medial axis into (npoints, 3)
    medial_axis = np.reshape(medial_axis_v, (npoints,ndim))
    # create a pyvista object of the medial axis points
    self.medial_axis = pvs.PolyData(medial_axis)
    t1 = perf_counter()
    print(f"-> compute_medial_axis() execution time: {t1-t0:g} (sec)")
    return

  def compute_covariance_eigv_sphere(self) -> np.ndarray:
    """
    Compute the covariance matrix and its eigen vectors for each point of the medial axis.

    :return:
      Array of the eigen vectors for each point of the medial axis.
      Shape is ``(npoints, 3, 3)`` where npoints is the number of points in the medial axis,

      - ``[:,0,:]`` is the first eigen vector, 
      - ``[:,1,:]`` is the second eigen vector and 
      - ``[:,2,:]`` is the third eigen vector.
    :rtype: np.ndarray
    """
    if self.medial_axis is None:
      self.compute_medial_axis()
    t0 = perf_counter()
    points = np.array(self.medial_axis.points, dtype=np.float32) # check if np.float64 needed
    if points.ndim != 2:
      raise RuntimeError(f'points coordinates must be of the shape (npoints, 3). Given shape: {points.shape}')
    # Get the number of points
    npoints = points.shape[0]

    # define min and max coords of our points
    O = np.array([ np.min(points[:,0]), np.min(points[:,1]), np.min(points[:,2]) ], dtype=np.float32)
    L = np.array([ np.max(points[:,0]), np.max(points[:,1]), np.max(points[:,2]) ], dtype=np.float32)

    # define how much cell we need to get cells of length radius
    # transform into an integer
    # warning python int() function returns the integer part of a number so add 1 to be sure to contain every point
    m = np.int64(( L - O )/self.radius_cov + 1)
    # create a mesh around our points
    mesh = CFEMeshQ1(3, m[0], m[1], m[2])
    # create the connectivity table
    mesh.create_e2v()
    # give it coords
    # WARNING: to ensure that all points are in the mesh maybe we need (mesh, O - radius, L + radius) To check
    fe_utils.femesh_set_uniform_coordinates(mesh, O, L)
    # Locate points in the mesh and attribute them an element index and local coords
    t2 = perf_counter()
    el_p,xi_p = fe_utils.femesh_point_location(mesh, points)
    t3 = perf_counter()
    print(f"femesh_point_location() located npoints: {npoints} in execution time: {t3-t2:g} (sec)")

    # Prepare function arguments
    p_coords_v = np.reshape(points, (npoints*3)) # flattened coordinates array
    msize      = np.array([mesh.m, mesh.n, mesh.p], dtype=np.int32) # mesh size
    e_vectors  = np.zeros(shape=(npoints*9), dtype=np.float64) # initialized flat array for eigen vectors
    
    # Call the c function that computes the covariance matrix and its eigen vectors
    gte.cfunc["utils"]["sphere_cov_eig_vectors"](
      np.int64(npoints), # number of points
      np.int32(mesh.ne), # number of cells
      msize,             # mesh size
      np.float64(p_coords_v), # flattened coordinates array
      np.int32(el_p), # cell indices of the points
      np.float32(self.radius_cov), # radius of research for the coviariance matrix
      e_vectors # eigen vectors returned by the function
    )
    # Reshape to get (npoints, 3, 3)
    e_vectors = np.reshape(e_vectors, (npoints,3,3))
    t1 = perf_counter()
    print(f"-> compute_covariance_eigv() execution time: {t1-t0:g} (sec)")
    return e_vectors
  
  def compute_covariance_eigv_knearest(self) -> np.ndarray:
    t0 = perf_counter()
    points  = np.array(self.medial_axis.points, dtype=np.float64)
    npoints = points.shape[0]
    # Prepare function arguments
    p_coords_v = np.reshape(points, (npoints*3)) # flattened coordinates array
    e_vectors  = np.zeros(shape=(npoints*9), dtype=np.float64) # initialized flat array for eigen vectors
    gte.cfunc["utils"]["cov_eig_vectors_knearest"](
      np.int64(npoints), # number of points
      np.int64(self.cov_knearest), # number of nearest neighbors to compute the covariance matrix
      np.float64(p_coords_v), # flattened coordinates array
      e_vectors # eigen vectors returned by the function
    )
    # Reshape to get (npoints, 3, 3)
    e_vectors = np.reshape(e_vectors, (npoints,3,3))
    t1 = perf_counter()
    print(f"-> compute_covariance_eigv() for {npoints} points, execution time: {t1-t0:g} (sec)")
    return e_vectors

  def compute_covariance_eigv(self) -> np.ndarray:
    return self.compute_covariance_eigv_knearest()

  def get_orientation_vectors(self):
    """
    Register the eigen vectors of the covariance matrix on the 
    :py:attr:`medial_axis <geoteqpy.MedialAxis.medial_axis>` mesh 
    if they are not already present.
    """
    if self.pca is None:
      return
    keys = list(self.medial_axis.point_data.keys())
    keys.extend(self.medial_axis.cell_data.keys())
    # check if the orientation vectors are already on the mesh
    if "eigv_0" in keys and "eigv_1" in keys and "eigv_2" in keys:
      # if they are already present, do nothing
      return
    if self.pca is not None:
      # compute the PCA
      eigvec = self.pca.compute_pca()
      # get the orientation vectors
      self.medial_axis["eigv_0"] = eigvec[:,0,:]
      self.medial_axis["eigv_1"] = eigvec[:,1,:]
      self.medial_axis["eigv_2"] = eigvec[:,2,:]
    return
