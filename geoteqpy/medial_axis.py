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
  :param str pca_method:
    Method to compute the covariance matrix and its eigenvectors.
    Available methods are ``"sphere"`` or ``"knearest"``.
    
    - ``"sphere"``: compute the covariance matrix using a sphere of radius ``cov_radius``.
    - ``"knearest"``: compute the covariance matrix using the ``cov_knearest`` nearest neighbors.
  
  :param float cov_radius:
    Radius of the sphere to select the points for the covariance matrix. Default is 1e4.
  :param int cov_knearest:
    Number of nearest neighbors to select the points for the covariance matrix. Default is 100.

  :Attributes:

  .. py:attribute:: medial_axis
    :type: pyvista.PolyData

    Pyvista mesh object of the medial axis points.
    
  .. py:attribute:: radius_ma
    :type: float

    Radius of research for the medial axis.

  .. py:attribute:: pca
    :type: geoteqpy.PCA

    PCA object to compute the covariance matrix and its eigenvectors.

  :Example:

  .. code:: python

    import pyvista as pvs
    import geoteqpy as gte

    # Load a mesh
    mesh = pvs.read("path/to/mesh.vtu")
    # Create a MedialAxis object, compute the medial axis and the orientation vectors
    ma = gte.MedialAxis(mesh=mesh,radius_ma=1e4,pca_method="knearest",cov_knearest=100)

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
