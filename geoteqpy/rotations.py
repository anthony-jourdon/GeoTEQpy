#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  GeoTEQpy
#  filename: rotations.py
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

import numpy as np

class Rotation:
  """
  Class to perform 3D rotations.

  :param float angle: Angle of rotation
  :param numpy.ndarray axis: Axis of rotation, expected shape ``(3,)``. The vector is expected to be **normalized**.
  
  :Attributes:

  .. py:attribute:: angle
    :type: float

    Angle of rotation

  .. py:attribute:: axis
    :type: numpy.ndarray

    Axis of rotation

  :Methods:
  
  """
  def __init__(self,angle:float,axis:np.ndarray) -> None:
    self.angle = angle
    self.axis  = axis
    return

  def rotation_matrix(self) -> np.ndarray:
    """
    Construct the 3D rotation matrix for an angle :math:`\\theta` and a rotation axis :math:`\\mathbf r`.
    
    :return: The rotation matrix. Shape ``(3,3)``
    :ret type: numpy.ndarray
    """
    theta = self.angle
    r     = self.axis
    R = np.array([ [np.cos(theta)+(1.0-np.cos(theta))*r[0]**2, 
                    r[0]*r[1]*(1.0-np.cos(theta))-r[2]*np.sin(theta), 
                    r[0]*r[2]*(1.0-np.cos(theta))+r[1]*np.sin(theta)],
                  [r[0]*r[1]*(1.0-np.cos(theta))+r[2]*np.sin(theta),
                    np.cos(theta)+(1.0-np.cos(theta))*r[1]**2,
                    r[1]*r[2]*(1.0-np.cos(theta))-r[0]*np.sin(theta)],
                  [r[0]*r[2]*(1.0-np.cos(theta))-r[1]*np.sin(theta), 
                    r[1]*r[2]*(1.0-np.cos(theta))+r[0]*np.sin(theta), 
                    np.cos(theta)+(1.0-np.cos(theta))*r[2]**2]])
    return R
  
  def rotate_vector(self,R:np.ndarray,u:np.ndarray,ccw:bool=True) -> np.ndarray:
    """
    Rotate vector(s) :math:`\\mathbf u` given the rotation matrix :math:`\\boldsymbol R`.

    .. warning:: 

      This is **not** a rotation of the vector field, but a rotation of the vectors themselves.

    :param np.ndarray R: rotation matrix of the shape ``(dim,dim)``
    :param np.ndarray u: vector(s) to be rotated of the shape ``(npoints,dim)``
    :param bool ccw: rotate counter-clockwise (default is True)

    :return: Rotated vector of the shape ``(npoints,dim)``
    """
    # u is expected to be in the form (npoints,dim)
    if ccw: # rotate couter-clockwise
      u_R = np.matmul(R,u.T).T
    else: # rotate clockwise
      u_R = np.matmul(R.T,u.T).T
    return u_R
  
  def rotate_tensor(self,R:np.ndarray,tensor:np.ndarray,ccw:bool=True) -> np.ndarray:
    """
    Rotate a tensor such that

    .. math::

      \\boldsymbol{\\tau}_R = \\boldsymbol{R} \\boldsymbol{\\tau} \\boldsymbol{R}^T

    :param np.ndarray R: rotation matrix of the shape ``(dim,dim)``
    :param np.ndarray tensor: tensor(s) to be rotated of the shape ``(npoints,dim,dim)``
    :param bool ccw: rotate counter-clockwise (default is True)

    :return: Rotated tensor of the shape ``(npoints,dim,dim)``
    """
    if tensor.ndim != 3:
      raise ValueError(f"Expected tensor to be of shape (npoints,3,3), found: {tensor.shape}")
    if ccw:
      tensor_r = (np.matmul(R.T,(np.matmul(tensor,R)).T)).T
    else:
      tensor_r = (np.matmul(R,(np.matmul(tensor,R.T)).T)).T
    return tensor_r

  def rotate_referential(self,coor:np.ndarray,O:np.ndarray,L:np.ndarray,ccw:bool=True) -> np.ndarray:
    """
    rotate_referential(self,coor,O,L,ccw=True)
    Rotate the referential of the coordinates :math:`\\mathbf{x}` given the rotation matrix :math:`\\boldsymbol R`.
    The referential is first translated to be centred on :math:`\\mathbf{0}`, then rotated and finally translated back to its original position.
    
    .. math:: 
      \\mathbf x_T &= \\mathbf x - \\frac{1}{2}(\\mathbf L + \\mathbf O) \\\\
      \\mathbf x_{TR} &= \\boldsymbol R \\mathbf x_T \\\\
      \\mathbf x_R &= \\mathbf x_{TR} + \\frac{1}{2}(\\mathbf L + \\mathbf O)

    :param np.ndarray coor: coordinates to be rotated of the shape ``(npoints,dim)``
    :param np.ndarray O: origin of the referential of the shape ``(dim,)``
    :param np.ndarray L: maximum coordinates of the referential of the shape ``(dim,)``
    :param bool ccw: rotate counter-clockwise (default is True)

    :return: Rotated coordinates of the shape ``(npoints,dim)``
    """
    # coor is expected to be in the form (npoints,dim)
    # get the rotation matrix
    R = self.rotation_matrix()
    # translate referential to get centred on 0
    coorT = coor - 0.5*(L + O)
    # rotate
    coorTR = self.rotate_vector(R,coorT,ccw)
    # translate back
    coorR = coorTR + 0.5*(L + O)
    return coorR

def test():
  r = Rotation(np.pi/2.0, np.array([1,0,0],dtype=np.float32))
  R = r.rotation_matrix()
  print(R)
  # test rotation of vector
  u = np.array([1,0,0],dtype=np.float32)
  u_R = r.rotate_vector(R,u)
  print(u_R)
  # test rotation of tensor
  tau9 = np.array([1,2,3,2,4,5,3,5,6],dtype=np.float32)
  tau  = np.zeros(shape=(2,9), dtype=np.float32)
  tau[0,:] = tau9
  tau[1,:] = tau9*2
  tau3x3 = np.reshape(tau, (-1,3,3))
  print(tau3x3)
  tt = r.rotate_tensor(R,tau3x3)
  print(tt)
  return

if __name__ == "__main__":
  test()
