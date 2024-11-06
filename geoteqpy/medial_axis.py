import pyvista as pvs
import numpy as np
import geoteqpy as gte
import pyviztools.utils as fe_utils
from pyviztools import CFEMeshQ1
from time import perf_counter

class MedialAxis:
  def __init__(self,mesh:pvs.UnstructuredGrid,radius_ma:float,radius_cov:float) -> None:
    self.mesh        = mesh
    self.radius_ma   = radius_ma
    self.radius_cov  = radius_cov
    self.medial_axis = None
    return
  
  def get_medial_axis_mesh(self,get_eigv:bool=True) -> pvs.PolyData:
    if self.medial_axis is None:
      self.compute_medial_axis()
    if get_eigv:
      self.get_orientation_vectors()
    return self.medial_axis

  def compute_medial_axis(self) -> pvs.PolyData:
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
    points_v = np.float64( np.reshape(points,  newshape=(ndim*npoints)) )
    normal_v = np.float64( np.reshape(normals, newshape=(ndim*npoints)) )
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
    medial_axis = np.reshape(medial_axis_v,newshape=(npoints,ndim))
    # create a pyvista object of the medial axis points
    self.medial_axis = pvs.PolyData(medial_axis)
    return
  
  def compute_covariance_eigv(self):
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
    p_coords_v = np.reshape(points, newshape=(npoints*3)) # flattened coordinates array
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
    e_vectors = np.reshape(e_vectors, newshape=(npoints,3,3))
    t1 = perf_counter()
    print(f"-> compute_covariance_eigv() execution time: {t1-t0:g} (sec)")
    return e_vectors

  def get_orientation_vectors(self):
    keys = list(self.medial_axis.point_data.keys())
    keys.extend(self.medial_axis.cell_data.keys())
    # check if the orientation vectors are already on the mesh
    if ("eigv_0" not in keys or 
        "eigv_1" not in keys or 
        "eigv_2" not in keys):
      eigvec = self.compute_covariance_eigv()
    # get the orientation vectors
    self.medial_axis["eigv_0"] = eigvec[:,0,:]
    self.medial_axis["eigv_1"] = eigvec[:,1,:]
    self.medial_axis["eigv_2"] = eigvec[:,2,:]
    return