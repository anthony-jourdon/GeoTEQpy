#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  GeoTEQpy
#  filename: register_cfunc.py
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

import os
import ctypes
from numpy.ctypeslib import ndpointer

def load_clib():
  
  vt_module_path = os.path.dirname(os.path.abspath(__file__))
  c_library_path = os.path.join(vt_module_path, "c", "lib")
  lib_name = "faulttools.so"
  for file in os.listdir(c_library_path):
    if file.startswith("faulttools"):
      lib_name = file
      break
  c_library_name = os.path.join(c_library_path, lib_name)

  clib = None
  libpath = os.path.dirname(c_library_name)
  libname = os.path.basename(c_library_name)
  
  # Append current working directry to LD_LIBRARY_PATH
  try:
    os.environ["LD_LIBRARY_PATH"] += os.pathsep + libpath
  except:
    os.environ["LD_LIBRARY_PATH"] = libpath

  try:
    clib = ctypes.cdll.LoadLibrary(os.path.join(libpath, libname))
    print('[load_clib]','Successfully loaded shared object:', libname)
  except:
    print('[load_clib]','Failed to load shared object:', libname)

  if clib is None:
    cfuncs = None
    return clib, cfuncs

  cfuncs = dict()
  cfuncs["utils"] = load_c_methods(clib)

  return clib, cfuncs

def load_c_methods(clib):

  function = dict()
  # ===========================================================================================
  c_method = clib.ft_compute_eigV
  c_method.restype = None 
  c_method.argtypes = [
                       ctypes.c_long,                                    # n_points
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # E[]
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # eig_vec[]
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")  # eig_val[]
                       ]
  function["eigv"] = c_method
  # ===========================================================================================
  c_method = clib.ft_compute_eigV_sorted
  c_method.restype = None 
  c_method.argtypes = [
                       ctypes.c_long,                                    # n_points
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # E[]
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # eig_vec[]
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")  # eig_val[]
                       ]
  function["eigv_sorted"] = c_method
  # ===========================================================================================
  c_method = clib.ft_compute_covariance_matrix_3d
  c_method.restype = None 
  c_method.argtypes = [
                       ctypes.c_long,                                   # n_points
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # data[]
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")  # cov_v[]
                       ]
  function["covariance_matrix_3d"] = c_method
  # ===========================================================================================
  c_method = clib.ft_compute_cov_eig_vectors
  c_method.restype = None 
  c_method.argtypes = [
                       ctypes.c_float,                                  # radius
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # point_coor[]
                       ctypes.c_long,                                   # npoints
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")  # vectors[]
                       ]
  function["cov_eig_vectors"] = c_method
  # ===========================================================================================
  c_method = clib.sphere_cov_eig_vectors
  c_method.restype = None 
  c_method.argtypes = [
                       ctypes.c_long,                                    # npoints
                       ctypes.c_int,                                     # ncells
                       ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),    # msize[]
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # p_coor[]
                       ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),    # p_cellidx[]
                       ctypes.c_float,                                   # radius
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")  # e_vectors[]
                       ]
  function["sphere_cov_eig_vectors"] = c_method
  # ===========================================================================================
  c_method = clib.compute_medial_axis_transform
  c_method.restype = None 
  c_method.argtypes = [
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # points[]
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # normals[]
                       ctypes.c_long,                                    # npoints
                       ctypes.c_double,                                  # radius_init
                       ctypes.c_double,                                  # angle_planar
                       ctypes.c_double,                                  # angle_preserve
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")  # medial_axis_points[]
                       ]
  function["medial_axis"] = c_method

  # ===========================================================================================
  c_method = clib.ft_apply_dc_patches
  c_method.restype = None
  c_method.argtypes = [
                       ctypes.c_long,                                    # npoints
                       ctypes.c_int,                                     # npatches
                       ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),  # coor[]
                       ctypes.c_float,                                   # radius0
                       ctypes.c_int,                                     # N0
                       ctypes.c_float,                                   # D
                       ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),  # Dc0[]
                       ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")   # Dc[]
                       ]
  function["apply_dc_patches"] = c_method
  return function
