#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  GeoTEQpy
#  filename: mesh_extrusion.py
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
import pyvista as pvs

class Extrusion:
  # TODO: add ne, face, and dx as class attributes
  def __init__(self, mesh:pvs.StructuredGrid) -> None:
    self.mesh = mesh
    # number of nodes in each direction
    self.nodes = np.array([mesh.dimensions[0], mesh.dimensions[1], mesh.dimensions[2]], dtype=np.int32)
    # number of cells in each direction (Q1 mesh is assumed)
    self.cells = np.array([mesh.dimensions[0]-1, mesh.dimensions[1]-1, mesh.dimensions[2]-1], dtype=np.int32)
    # min coords in each direction
    self.O = np.array([mesh.bounds[0], mesh.bounds[2], mesh.bounds[4]], dtype=np.float32)
    # max coords in each direction
    self.L = np.array([mesh.bounds[1], mesh.bounds[3], mesh.bounds[5]], dtype=np.float32)
    self.faces_info = self.faces_indicators()
    return

  def faces_indicators(self) -> dict:
    """
    Create a dictionnary with the following face informations:

      - ``axis``: the axis of the face (x:0, y:1, z:2)
      - ``node``: the node index of the face in its normal direction
      - ``cell``: the cell index of the face in its normal direction
      - ``scale``: ``+1`` if the face is a max face, ``-1`` if the face is a min face
    
    :return: dictionnary with face informations
    :rtype: dict
    """
    directions = ['x','y','z']
    limits     = ['min','max']
    faces = {}

    for d in directions:
      for l in limits:
        idx  = directions.index(d)
        face = d+l
        axis = idx
        if l == 'min': 
          node = 0
          cell = 0
          scale = -1
        else: 
          node = self.nodes[idx]-1
          cell = node-1
          scale = 1
        faces[face] = {
          'axis': axis,
          'node': node,
          'cell': cell,
          'scale': scale
        }
    return faces
  
  def extrude_coordinates(self,ne:int,dx:float,face:str) -> list[np.ndarray]:
    """
    
    """
    face_info = self.faces_info[face]
    axis      = face_info['axis']
    
    # newshape for the extruded mesh, will be modified depending on the face
    shape  = [self.nodes[0], self.nodes[1], self.nodes[2]]
    # slices for the extruded mesh, will be modified depending on the face
    slices = [slice(None), slice(None), slice(None)]

    layers = []
    for d in range(3):
      # get the coordinates in the direction d
      x = np.reshape(self.mesh.points[:,d], newshape=(self.nodes[0], self.nodes[1], self.nodes[2]), order='F')
      shape[axis] = 1
      lx = np.reshape(np.take(x,face_info['node'],axis=axis), newshape=tuple(shape))
      shape[axis] = ne
      layer_x = np.zeros(shape=tuple(shape), dtype=np.float32)
      for j in range(ne):
        slices[axis] = j # [:,:,j] or [:,j,:] or [j,:,:]
        sl0 = tuple(slices)
        if 'max' in face:
          slices[axis] = j # [:,:,j] or [:,j,:] or [j,:,:]
          sl1 = tuple(slices)
          slices[axis] = j-1 # [:,:,j-1] or [:,j-1,:] or [j-1,:,:]
          sl2 = tuple(slices)
        else:
          slices[axis] = ne - j - 1 # [:,:,ne - j - 1] or [:,ne - j - 1,:] or [ne - j - 1,:,:]
          sl1 = tuple(slices)
          slices[axis] = ne - j # [:,:,ne - j] or [:,ne - j,:] or [ne - j,:,:]
          sl2 = tuple(slices)
        slices[axis] = 0 # [:,:,0] or [:,0,:] or [0,:,:]
        sl3 = tuple(slices)
        if d == axis: # if the direction is the direction of extrusion
          if j == 0: layer_x[sl1] = lx[sl0] + dx * face_info['scale'] # first layer, take the coords from original mesh
          else:      layer_x[sl1] = layer_x[sl2] + dx * face_info['scale'] # following layers, add to the previous layer
        else:
          layer_x[sl0] = lx[sl3] # other directions, copy the layer
      layers.append(layer_x)
    return layers

  def extrude_cellfield(self,field_key:str,nc:int,face:str,bs:int=1) -> np.ndarray:
    face_info = self.faces_info[face]
    axis      = face_info['axis']
    # newshape for the extruded field, will be modified depending on the face
    shape  = [self.cells[0], self.cells[1], self.cells[2]]
    if bs != 1:
      shape = [self.cells[0], self.cells[1], self.cells[2], bs]
    # slices for the extruded field, will be modified depending on the face
    slices = [slice(None), slice(None), slice(None)]

    field_i = np.reshape(self.mesh.cell_data[field_key], newshape=tuple(shape), order='F')
    field = np.take(field_i, face_info['cell'], axis=axis)

    shape[axis] = nc
    layer_field = np.zeros(shape=tuple(shape), dtype=field.dtype)
    for j in range(nc):
      slices[axis] = j
      layer_field[tuple(slices)] = field
    if 'max' in face: extruded_field = np.concatenate([field_i, layer_field], axis=axis)
    else:             extruded_field = np.concatenate([layer_field, field_i], axis=axis)
    return extruded_field
  
  def merge_cellfield(self,field_key:str,extruded_mesh:pvs.StructuredGrid,ne:int,face:str) -> None:
    # default to scalar field
    block_size = 1
    shape      = [extruded_mesh.n_cells]
    sl0 = [slice(None)]
    sl1 = [slice(None), slice(None), slice(None)]
    # vector field
    if self.mesh.cell_data[field_key].ndim == 2:
      block_size = self.mesh.cell_data[field_key].shape[1]
      shape.append(block_size) # shape = [n_cells, block_size]
      sl0.append(slice(None))
      sl1.append(slice(None))
    # extrude field
    extruded_field = self.extrude_cellfield(field_key,ne,face,block_size)

    extruded_mesh.cell_data[field_key] = np.zeros(shape=tuple(shape), dtype=self.mesh.cell_data[field_key].dtype)
    for dof in range(block_size):
      if block_size != 1: 
        sl0[-1] = dof
        sl1[-1] = dof
      extruded_mesh.cell_data[field_key][tuple(sl0)] = np.reshape(
        extruded_field[tuple(sl1)], 
        newshape=(
          (extruded_mesh.dimensions[0] - 1) * 
          (extruded_mesh.dimensions[1] - 1) * 
          (extruded_mesh.dimensions[2] - 1)
        ), 
        order='F'
      )
    return

  def extrude_pointfield(self,field_key:str,ne:int,face:str,bs:int=1) -> np.ndarray:
    face_info = self.faces_info[face]
    axis      = face_info['axis']

    # newshape for the extruded field, will be modified depending on the face
    shape  = [self.nodes[0], self.nodes[1], self.nodes[2]]
    if bs != 1:
      shape = [self.nodes[0], self.nodes[1], self.nodes[2], bs]
    # slices for the extruded field, will be modified depending on the face
    slices = [slice(None), slice(None), slice(None)]

    field_i = np.reshape(self.mesh.point_data[field_key], newshape=tuple(shape), order='F')
    field = np.take(field_i, face_info['node'], axis=axis)

    shape[axis] = ne
    layer_field = np.zeros(shape=tuple(shape), dtype=field.dtype)
    for j in range(ne):
      slices[axis] = j
      layer_field[tuple(slices)] = field
    if 'max' in face: extruded_field = np.concatenate([field_i, layer_field], axis=axis)
    else:             extruded_field = np.concatenate([layer_field, field_i], axis=axis)   
    return extruded_field
  
  def merge_pointfield(self,field_key:str,extruded_mesh:pvs.StructuredGrid,ne:int,face:str) -> None:
    # default to scalar field
    block_size = 1
    shape      = [extruded_mesh.n_points]
    sl0 = [slice(None)]
    sl1 = [slice(None), slice(None), slice(None)]
    # vector field
    if self.mesh.point_data[field_key].ndim == 2:
      block_size = self.mesh.point_data[field_key].shape[1]
      shape.append(block_size) # shape = [n_points, block_size]
      sl0.append(slice(None))
      sl1.append(slice(None))
    # extrude field
    extruded_field = self.extrude_pointfield(field_key,ne,face,block_size)
    
    extruded_mesh.point_data[field_key] = np.zeros(shape=tuple(shape), dtype=self.mesh.point_data[field_key].dtype)
    for dof in range(block_size):
      if block_size != 1: 
        sl0[-1] = dof
        sl1[-1] = dof
      extruded_mesh.point_data[field_key][tuple(sl0)] = np.reshape(
        extruded_field[tuple(sl1)], 
        newshape=(
          extruded_mesh.dimensions[0] * 
          extruded_mesh.dimensions[1] * 
          extruded_mesh.dimensions[2]
        ), 
        order='F'
      )
    return

  def extrude_mesh(self,ne:int,dx:float,face:str,fields:list[str]|None=None) -> pvs.StructuredGrid:
    extruded_coords = self.extrude_coordinates(ne,dx,face)
    dirmap = {0: self.mesh.x, 1: self.mesh.y, 2: self.mesh.z}

    axis = self.faces_info[face]['axis']
    coor = []
    for d in range(3):
      if 'max' in face: order = [dirmap[d], extruded_coords[d]]
      else:             order = [extruded_coords[d], dirmap[d]] 
      coor.append(np.concatenate([*order], axis=axis))
    extruded_mesh = pvs.StructuredGrid(*coor)

    if fields is not None:
      for field in fields:
        if field in self.mesh.cell_data:
          self.merge_cellfield(field,extruded_mesh,ne,face)
        if field in self.mesh.point_data:
          self.merge_pointfield(field,extruded_mesh,ne,face)
    return extruded_mesh

def test_extrude_coor():
  x_ = np.linspace(-2,2,6, dtype=np.float32)
  y_ = np.linspace(-3,3,8, dtype=np.float32)
  z_ = np.linspace(-1,1,4, dtype=np.float32)
  x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')
  mesh = pvs.StructuredGrid(x, y, z)

  np.random.seed(0)
  mesh['e2'] = np.random.rand(mesh.n_cells)
  mesh['e3'] = np.zeros(mesh.n_points, dtype=np.float32)
  mesh['e3'][mesh.points[:,2] >= 0] = 1.0
  mesh['e3'][mesh.points[:,0] <= 0] = 2.0

  mesh['e4'] = np.zeros(shape=(mesh.n_points, 3), dtype=np.float32)
  mesh['e4'][mesh.points[:,2] >= 0, 0] = 1.0
  mesh['e4'][mesh.points[:,0] <= 0, 2] = 2.0

  extrusion = Extrusion(mesh)
  extruded_mesh = extrusion.extrude_mesh(5,0.2,'ymin',['e2','e3','e4'])
  
  extrusion = Extrusion(extruded_mesh)
  extruded_mesh = extrusion.extrude_mesh(3,0.1,'zmax',['e2','e3','e4'])
  extrusion = Extrusion(extruded_mesh)
  extruded_mesh = extrusion.extrude_mesh(4,0.1,'xmax',['e2','e3','e4'])
  extrusion = Extrusion(extruded_mesh)
  extruded_mesh = extrusion.extrude_mesh(4,0.1,'xmin',['e2','e3','e4'])

  p = pvs.Plotter()
  p.add_mesh(extruded_mesh, show_edges=True, edge_color='black', color='green', scalars="e4",)
  p.add_axes()
  p.show()
  return

if __name__ == '__main__':
  test_extrude_coor()
