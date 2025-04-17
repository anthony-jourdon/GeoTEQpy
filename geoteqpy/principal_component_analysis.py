import numpy as np
import geoteqpy as gte
import pyviztools.utils as fe_utils
from pyviztools import CFEMeshQ1
from time import perf_counter

class PCA:
  """
  Class to compute the Principal Component Analysis (PCA) of a set of points in 3D space.
  The PCA is computed using the covariance matrix of the points.

  This class is the Parent class of the PCA classes. It is not meant to be used directly. 
  """
  def __init__(self, points):
    self.points  = np.array(points, dtype=np.float64)
    if self.points.ndim != 2:
      raise RuntimeError(f'points coordinates must be of the shape (npoints, 3). Given shape: {self.points.shape}')
    self.npoints = self.points.shape[0]
    return
  
  def compute_pca(self) -> np.ndarray:
    raise NotImplementedError(f"Class {self.__class__.__name__} does not contain compute_pca() method. Must use one of the derived classes.")
  
class SpherePCA(PCA):
  """
  Class to compute the Principal Component Analysis (PCA) of a set of points in 3D space 
  using a sphere to select the points that contribute to the covariance matrix.
  """
  def __init__(self, points, radius):
    super().__init__(points)
    self.radius = radius
    return
  
  def compute_pca(self) -> np.ndarray:
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
    t0 = perf_counter()
    # define min and max coords of our points
    O = np.array([ np.min(self.points[:,0]), np.min(self.points[:,1]), np.min(self.points[:,2]) ], dtype=np.float32)
    L = np.array([ np.max(self.points[:,0]), np.max(self.points[:,1]), np.max(self.points[:,2]) ], dtype=np.float32)

    # define how much cell we need to get cells of length radius
    # transform into an integer
    # warning python int() function returns the integer part of a number so add 1 to be sure to contain every point
    m = np.int64(( L - O )/self.radius + 1)
    # create a mesh around our points
    mesh = CFEMeshQ1(3, m[0], m[1], m[2])
    # create the connectivity table
    mesh.create_e2v()
    # give it coords
    # WARNING: to ensure that all points are in the mesh maybe we need (mesh, O - radius, L + radius) To check
    fe_utils.femesh_set_uniform_coordinates(mesh, O, L)
    # Locate points in the mesh and attribute them an element index and local coords
    t2 = perf_counter()
    el_p,_ = fe_utils.femesh_point_location(mesh, np.float32(self.points))
    t3 = perf_counter()
    print(f"femesh_point_location() located npoints: {self.npoints} in execution time: {t3-t2:g} (sec)")

    # Prepare function arguments
    p_coords_v = np.reshape(self.points, (self.npoints*3)) # flattened coordinates array
    msize      = np.array([mesh.m, mesh.n, mesh.p], dtype=np.int32) # mesh size
    e_vectors  = np.zeros(shape=(self.npoints*9), dtype=np.float64) # initialized flat array for eigen vectors
    
    # Call the c function that computes the covariance matrix and its eigen vectors
    gte.cfunc["utils"]["sphere_cov_eig_vectors"](
      np.int64(self.npoints), # number of points
      np.int32(mesh.ne), # number of cells
      msize,             # mesh size
      np.float64(p_coords_v), # flattened coordinates array
      np.int32(el_p), # cell indices of the points
      np.float32(self.radius), # radius of research for the coviariance matrix
      e_vectors # eigen vectors returned by the function
    )
    # Reshape to get (npoints, 3, 3)
    e_vectors = np.reshape(e_vectors, (self.npoints,3,3))
    t1 = perf_counter()
    print(f"-> compute_pca() for {self.npoints} points, execution time: {t1-t0:g} (sec)")
    return e_vectors

class KNearestPCA(PCA):
  """
  Class to compute the Principal Component Analysis (PCA) of a set of points in 3D space 
  using the K-nearest neighbours points to select the points that contribute to the covariance matrix.
  """
  def __init__(self, points, knearest):
    super().__init__(points)
    self.knearest = knearest
    return
  
  def compute_pca(self):
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
    t0 = perf_counter()
    # Prepare function arguments
    p_coords_v = np.reshape(self.points, (self.npoints*3)) # flattened coordinates array
    e_vectors  = np.zeros(shape=(self.npoints*9), dtype=np.float64) # initialized flat array for eigen vectors
    gte.cfunc["utils"]["cov_eig_vectors_knearest"](
      np.int64(self.npoints), # number of points
      np.int64(self.knearest), # number of nearest neighbors to compute the covariance matrix
      np.float64(p_coords_v), # flattened coordinates array
      e_vectors # eigen vectors returned by the function
    )
    # Reshape to get (npoints, 3, 3)
    e_vectors = np.reshape(e_vectors, (self.npoints,3,3))
    t1 = perf_counter()
    print(f"-> compute_pca() for {self.npoints} points, execution time: {t1-t0:g} (sec)")
    return e_vectors