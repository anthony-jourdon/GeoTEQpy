import os
import numpy as np
import pyvista as pvs

def find_last_false_idx(arr):
  """
  Search backward the first index of the point containing a non nan value
  """
  for i in range(len(arr)-1, -1, -1):
    if arr[i] == False:
      return i
  return None

def fill_nan_j(field:np.ndarray,size:np.ndarray) -> np.ndarray:
  fcopy = field.copy()
  
  # get boolean array: if nan => True else => False
  nan_mask = np.isnan(fcopy)
  # for each point, check wether any of the component evaluate to nan
  mask = np.any(nan_mask, axis=1)
  # get the indices of the nan values in a flatten array
  nan_indices = np.argwhere(mask).ravel()
  # reshape mask to 3d
  mask_3d = np.reshape(mask,newshape=(size[2],size[1],size[0]))
  # create a 2d array to hold indices of the non nan values along y axis 
  idx = np.zeros(shape=(mask_3d[:,0,:].shape), dtype=np.int32)
  # search for the indices of the non nan values along y axis
  for k in range(size[2]):
    for i in range(size[0]):
      j = find_last_false_idx(mask_3d[k,:,i])
      idx[k,i] = k*size[1]*size[0] + j*size[0] + i
  
  # flatten the array of non nan indices
  top_idx = np.ravel(idx)

  fcopy[nan_indices] = fcopy[top_idx]

  return fcopy

def symtensor_get_unique_components(tensor:np.ndarray) -> np.ndarray:
  """ 
  symtensor_get_unique_components(tensor)
  Get the unique components of a symmetric tensor
  Order of the outputed components: (0,4,8,1,2,5) -> (xx,yy,zz,xy,xz,yz)
  Parameters:
    tensor: the components of the symmetric tensor of shape (npoints,9)
  
  Output:
    the unique components of the symmetric tensor of shape (npoints,6)
  """
  return tensor[:,[0,4,8,1,2,5]]

def compute_xi(mesh:pvs.StructuredGrid,e2_key:str,field_name:str="xi") -> None:
  mesh[field_name] = np.exp( np.log10(mesh[e2_key]) - np.min( np.log10(mesh[e2_key]) ) )
  return

def replace_extension(fname:str,new_extension:str) -> str:
  """
  Replace the extension of a file name
  """
  f,_ = os.path.splitext(fname)
  return f + new_extension