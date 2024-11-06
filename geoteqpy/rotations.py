import numpy as np

class Rotation:
  def __init__(self,angle:float,axis:np.ndarray) -> None:
    self.angle = angle
    self.axis  = axis
    return

  def rotation_matrix(self) -> np.ndarray:
    """
    rotation_matrix(theta, r)
    Construct the 3D rotation matrix for an angle theta and a rotation axis r

    Parameters:
      theta: angle of rotation in radians
      r:     axis of rotation. A 3D unit vector is expected (i.e., ||r|| = 1).

    Output:
      R: the 3D rotation matrix of the shape (3,3) 
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
  
  def rotate_vector(self,R:np.ndarray,u:np.ndarray,ccw:bool=True):
    """
    rotate_vector(self,R,u,ccw=True)
    Rotate vector(s) :math:`\\mathbf u` given the rotation matrix :math:`\\boldsymbol R`.

    .. warning:: 

      This is **not** a rotation of the vector field, but a rotation of the vectors themselves.
      To rotate the vector field, have a look at how it is done in
      :py:meth:`evaluate_velocity_symbolic() <genepy.VelocityLinear.evaluate_velocity_symbolic>`.

    :param np.ndarray R: rotation matrix of the shape ``(dim,dim)``
    :param np.ndarray u: vector(s) to be rotated of the shape ``(npoints,dim)``
    :param bool ccw: rotate counter-clockwise (default is True)

    :return: **u_R** rotated vector(s) of the shape ``(npoints,dim)``
    """
    # u is expected to be in the form (npoints,dim)
    if ccw: # rotate couter-clockwise
      u_R = np.matmul(R,u.T).T
    else: # rotate clockwise
      u_R = np.matmul(R.T,u.T).T
    return u_R
  
  def rotate_tensor(self,R:np.ndarray,tensor:np.ndarray,ccw:bool=True) -> np.ndarray:
    if tensor.ndim != 3:
      raise ValueError(f"Expected tensor to be of shape (npoints,3,3), found: {tensor.shape}")
    if ccw:
      tensor_r = (np.matmul(R.T,(np.matmul(tensor,R)).T)).T
    else:
      tensor_r = (np.matmul(R,(np.matmul(tensor,R.T)).T)).T
    return tensor_r

  def rotate_referential(self,coor:np.ndarray,O:np.ndarray,L:np.ndarray,ccw:bool=True):
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

    :return: **coorR** rotated coordinates of the shape ``(npoints,dim)``
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