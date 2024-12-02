#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  GeoTEQpy
#  filename: asagi.py
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
import netCDF4 as nc

class Asagi:
  """
  Class to write and read data in the `ASAGI`_ format.
  `ASAGI`_ uses the `netCDF4 <https://unidata.github.io/netcdf4-python/>`_ format.

  :param np.ndarray x:
    1D array of the x-coordinates. Shape ``(nx,)``.
  :param np.ndarray y:
    1D array of the y-coordinates. Shape ``(ny,)``.
  :param np.ndarray z:
    1D array of the z-coordinates. Shape ``(nz,)``.
  :param dict fields:
    Dictionary of the fields to write. 
    The keys are the names of the fields and the values are the data arrays.
    Data arrays must be of shape ``(nz,ny,nx)``.
  :param str fname:
    Name of the netCDF file to read.

  :Attributes:

  .. py:attribute:: x
    :type: np.ndarray
  
    1D array of the x-coordinates. Shape ``(nx,)``. Default: ``None``.

  .. py:attribute:: y
    :type: np.ndarray

    1D array of the y-coordinates. Shape ``(ny,)``. Default: ``None``.

  .. py:attribute:: z
    :type: np.ndarray
      
    1D array of the z-coordinates. Shape ``(nz,)``. Default: ``None``.
  
  .. py:attribute:: fields
    :type: dict

    Dictionary of the fields to write. 
    The keys are the names of the fields and the values are the data arrays.
    Data arrays must be of shape ``(nz,ny,nx)``. Default: ``None``.

  .. py:attribute:: shape
    :type: tuple

    Shape of the data arrays. Default: ``None``, if a file is passed, the shape is set from the coordinates in the file.

  .. py:attribute:: fname
    :type: str

    Name of the netCDF file to read. Default: ``None``.

    
  :Example:

  Read data from a file, coordinates and fields are stored in the class attributes:

  >>> asagi = Asagi(fname='test_netcdf.nc')

  Create an instance of the class with coordinates and fields:

  .. code:: python

    O = np.array([0,-1,0],dtype=np.float32)
    L = np.array([3,1,2],dtype=np.float32)
    n = np.array([3,4,5],dtype=np.int32)
    
    cds = []
    for i in range(3):
      cds.append(np.linspace(O[i],L[i],n[i]))
    
    density = np.ones(shape=(n[2],n[1],n[0]), dtype=np.float32)*3000.0
    fields = {"density": density}
    asagi = Asagi(cds[0],cds[1],cds[2],fields)

  Write the data to a netCDF file for `ASAGI`_:

  >>> asagi.write_netcdf('test_asagi.nc')
    
  Write the data to a netCDF file for `ParaView <https://www.paraview.org/>`_:

  >>> asagi.write_paraview('test_paraview.nc')

  """
  def __init__(self,
               x:np.ndarray|None=None, 
               y:np.ndarray|None=None, 
               z:np.ndarray|None=None, 
               fields:dict|None=None,
               fname:str|None=None,) -> None:
    # coordinates
    self.x = x
    self.y = y
    self.z = z
    # variables
    self.fields = fields
    # grid size
    self.shape = None
    if x is not None and y is not None and z is not None:
      self.shape = (z.shape[0],y.shape[0],x.shape[0])
    # if a filename is provided, read the netcdf file
    self.fname = fname
    if self.fname is not None:
      print(f"Reading data from file: {self.fname}")
      self.read_netcdf(self.fname)
    if self.fields is not None:
      for key in self.fields:
        if not self.check_size(self.fields[key]):
          serr = f"Field '{key}' is expected to be of shape: {self.shape}, found: {self.fields[key].shape}"
          raise ValueError(serr)
    return
  
  def __str__(self) -> str:
    s = f"class {self.__class__.__name__}\n"
    if self.fname is not None:
      s += f"Data from file: {self.fname}\n"
    s += f"  size:\n"
    s += f"    [ x: {self.shape[2]}, y:{self.shape[1]}, z: {self.shape[0]} ]\n"
    s += f"  shape: {self.shape}\n"
    if self.fields is not None:
      s += f"  fields:\n"
      for key in self.fields:
        s += f"    {key}, {self.fields[key].shape}, {self.fields[key].dtype}\n"
    return s
  
  def check_size(self,field:np.ndarray) -> bool:
    """
    Check if the given field has the same shape as the coordinates.

    :param np.ndarray field:
      Field to check.
    :return:
      ``True`` if the field has the same shape as the coordinates, ``False`` otherwise.
    :rtype: bool
    """
    if field.shape == self.shape:
      return True
    return False

  def get_fields_dtype(self) -> np.dtype:
    """
    Get the data type of the fields and store it in a numpy dtype object with the following format:
    ``[ ("field1",dtype1), ("field2",dtype2), ... ]`` 
    where ``"field1"`` are the keys of the :py:attr:`field <geoteqpy.MedialAxis.medial_axis>` dictionary 
    and ``dtype1`` are the data type of the corresponding keys.

    :return:
      Numpy data type object.
    :rtype: np.dtype
    """
    data_type = []
    for key,value in self.fields.items():
      data_type.append((key,value.dtype))
    return np.dtype(data_type)
  
  def fields_to_tuple(self, dtype:np.dtype) -> np.ndarray:
    """
    Reshape the fields to a tuple of the given data type.
  
    :param np.dtype dtype:
      Numpy data type object. 
      Use the :py:meth:`get_fields_dtype <geoteqpy.Asagi.get_fields_dtype>` method to get the data type.
    :return:
      Numpy array of the reshaped fields.
    :rtype: np.ndarray
    """
    reshaped_fields = np.empty(shape=self.shape, dtype=dtype)
    for field in self.fields:
      if field in dtype.names:
        reshaped_fields[:,:,:][field] = self.fields[field]
    return reshaped_fields

  def write_netcdf(self,fname:str,dimensions=('z','y','x')) -> None:
    """
    Write the data to a netCDF file for `ASAGI`_.

    :param str fname:
      Name of the netCDF file to write.
    :param tuple dimensions:
      Tuple of the dimensions names. Default is ``('z','y','x')``.
    """
    fp = nc.Dataset(fname,'w',format='NETCDF4')
    dims = {0: {'name': 'x', 'size': self.x.shape[0], 'data': self.x},
            1: {'name': 'y', 'size': self.y.shape[0], 'data': self.y},
            2: {'name': 'z', 'size': self.z.shape[0], 'data': self.z}}
    # coordinates
    for i in range(3):
      dim = dims[i]
      fp.createDimension(dim['name'],dim['size'])
      var = fp.createVariable(dim['name'],dim['data'].dtype,(dim['name'],))
      var[:] = dim['data']
    # fields
    material_types = self.get_fields_dtype()
    nc_dtypes = fp.createCompoundType(material_types,'material')
    nc_fields = fp.createVariable('data',nc_dtypes,dimensions)
    nc_fields[:,:,:] = self.fields_to_tuple(material_types)
    fp.close()
    return
  
  def write_paraview(self,fname:str) -> None:
    """
    Write the data to a netCDF file for `ParaView <https://www.paraview.org/>`_.

    :param str fname:
      Name of the netCDF file to write.
    """
    fp = nc.Dataset(fname,'w',format='NETCDF4')
    dims = {0: {'name': 'x', 'size': self.x.shape[0], 'data': self.x},
            1: {'name': 'y', 'size': self.y.shape[0], 'data': self.y},
            2: {'name': 'z', 'size': self.z.shape[0], 'data': self.z}}
    # coordinates
    for i in range(3):
      dim = dims[i]
      fp.createDimension(dim['name'],dim['size'])
      var = fp.createVariable(dim['name'],dim['data'].dtype,(dim['name'],))
      var[:] = dim['data']
    # fields
    for field in self.fields:
      var = fp.createVariable(field,self.fields[field].dtype,('z','y','x'))
      var[:,:,:] = self.fields[field]
    fp.close()
    return
  
  def read_netcdf(self,fname:str) -> None:
    """
    Read the data from a netCDF file.
    Set the class attributes values from the read netCDF file.

    :param str fname:
      Name of the netCDF file to read.
    """
    fp = nc.Dataset(fname,'r')
    self.x = fp.variables['x'][:]
    self.y = fp.variables['y'][:]
    self.z = fp.variables['z'][:]
    self.shape = (self.z.shape[0],self.y.shape[0],self.x.shape[0])
    self.fields = {}
    for key in fp.variables['data'].dtype.names:
      self.fields[key] = fp.variables['data'][:,:,:][key]
    fp.close()
    return

def test():
  O = np.array([0,-1,0],dtype=np.float32)
  L = np.array([3,1,2],dtype=np.float32)
  n = np.array([3,4,5],dtype=np.int32)
  
  cds = []
  for i in range(3):
    cds.append(np.linspace(O[i],L[i],n[i]))

  coor = np.zeros(shape=(n[2]*n[1]*n[0],3),dtype=np.float32)
  for k in range(n[2]):
    for j in range(n[1]):
      for i in range(n[0]):
        nidx = k*n[1]*n[0] + j*n[0] + i
        coor[nidx,0] = cds[0][i]
        coor[nidx,1] = cds[1][j]
        coor[nidx,2] = cds[2][k]

  density = np.ones(shape=(n[2]*n[1]*n[0]),dtype=np.float32)*3000.0
  regions = np.zeros(shape=(n[2]*n[1]*n[0]),dtype=np.int32)
  sxx     = np.ones(shape=(n[2]*n[1]*n[0]),dtype=np.float32)

  density[ coor[:,1] > 0.0 ] = 2700.0
  regions[ coor[:,1] > 0.0 ] = 1
  sxx[ coor[:,1] > 0.0 ]     = 2.0
  stress = { "xx": sxx, "yy": 2*sxx, "zz": 3*sxx, "xy": 4*sxx, "xz": 5*sxx, "yz": 6*sxx }
  fields = {"density": density, 
            "regions": regions, 
            "s_xx": stress['xx'],
            "s_yy": stress['yy'],
            "s_zz": stress['zz'],
            "s_xy": stress['xy'],
            "s_xz": stress['xz'],
            "s_yz": stress['yz'] }

  for field in fields:
    fields[field] = np.reshape(fields[field],newshape=(n[2],n[1],n[0]))

  nc_fname = 'test_netcdf.nc'
  asagi = Asagi(cds[0],cds[1],cds[2],fields)
  print(asagi)
  asagi.write_netcdf(nc_fname)
  asagi_r = Asagi(fname=nc_fname)
  print(asagi_r)
  asagi_r.write_paraview('test_paraview.nc')
  return

if __name__ == '__main__':
  test()
