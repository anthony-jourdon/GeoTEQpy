from geoteqpy import ModelData
from geoteqpy import rotations
from geoteqpy import Asagi
import geoteqpy.field_utils as gtfu
import numpy as np

class pTatin2SeisSol:
  def __init__(self,model:ModelData,uniform_mesh_size:np.ndarray,**fields_key) -> None:
    self.model      = model
    self.fields_key = fields_key
    self.rotation   = rotations.Rotation(angle=np.pi/2.0,axis=np.array([1,0,0],dtype=np.float32))
    self.R          = self.rotation.rotation_matrix()

    self.model.check_ptatin_mesh()
    self.size = uniform_mesh_size
    return
  
  def uniform_reshape(self,field):
    # reshape for asagi, note the swapaxes to intervert y and z
    return np.reshape(field,newshape=(self.size[2]+1,self.size[1]+1,self.size[0]+1)).swapaxes(0,1)

  def get_1d_coordinates(self,coor):
    x = self.uniform_reshape(coor[:,0])
    y = self.uniform_reshape(coor[:,1])
    z = self.uniform_reshape(coor[:,2])

    x_1d = x[0,0,:]
    y_1d = y[0,:,0]
    z_1d = z[:,0,0]
    return x_1d,y_1d,z_1d

  def rotate_stress(self,stress:np.ndarray):
    # shape for rotation
    stress = np.reshape(stress,(stress.shape[0],3,3))
    rotated_stress = self.rotation.rotate_tensor(R=self.R,tensor=stress)
    # shape for output
    rotated_stress = np.reshape(rotated_stress,(rotated_stress.shape[0],9))
    return rotated_stress
  
  def export_to_asagi(self,fname,fields_to_export:list[str]|None=None,pvoutput:str|None=None):
    # 1: we create the uniform mesh
    self.model.uniform_mesh = self.model.mesh_create_from_options(self.model.O,self.model.L,self.size)
    # 2: we compute stress
    self.model.ptatin_mesh.cell_data[self.fields_key['total_stress_key']] = self.model.get_total_stress(**self.fields_key)
    # 3: we interpolate the fields from ptatin to the newly created mesh
    if fields_to_export is None:
      fields_to_export = []
      for field in self.fields_key:
        fields_to_export.append(self.fields_key[field])
    self.model.interpolate_to_uniform_mesh(fields=fields_to_export)
    # 4: remove nan that could have been introduced by the interpolation
    for field in fields_to_export:
      self.model.uniform_mesh.point_data[field] = gtfu.fill_nan_j(self.model.uniform_mesh.point_data[field],self.size+1)
    # 5: we rotate the fields needing rotation
    coor_r   = self.rotation.rotate_referential(coor=self.model.uniform_mesh.coor,O=self.model.O,L=self.model.L)
    stress_r = self.rotate_stress(self.model.uniform_mesh.point_data[self.fields_key['total_stress_key']])
    # 6: reshape the fields
    asagi_fields = {}
    for field in fields_to_export:
      # special case for stress
      if field == self.fields_key['total_stress_key']:
        asagi_fields["s_xx"] = self.uniform_reshape(stress_r[:,0])
        asagi_fields["s_yy"] = self.uniform_reshape(stress_r[:,4])
        asagi_fields["s_zz"] = self.uniform_reshape(stress_r[:,8])
        asagi_fields["s_xy"] = self.uniform_reshape(stress_r[:,1])
        asagi_fields["s_xz"] = self.uniform_reshape(stress_r[:,2])
        asagi_fields["s_yz"] = self.uniform_reshape(stress_r[:,5])
      else:
        asagi_fields[field] = self.uniform_reshape(self.model.uniform_mesh.point_data[field])
    # 7: get 1d coordinates
    x,y,z = self.get_1d_coordinates(coor_r)
    # 8: export to asagi
    asagi = Asagi(x=x,y=y,z=z,fields=asagi_fields)
    print(asagi)
    asagi.write_netcdf(fname)
    # 9: write paraview readable netcdf file if required
    if pvoutput is not None:
      asagi.write_paraview(pvoutput)
    return
  
def test():
  import os
  import pyvista as pvs

  mdir = os.path.join(os.environ['M3D_DIR'],'StrikeSlip3D/MiddleRes/EastAnatolia')
  vts  = os.path.join(mdir,'step000825-cellfields.vts')
  fault = os.path.join(mdir,'Postproc/FaultSurface/step825/Main_SZ_surface_0.vtp')
  out_netcdf = os.path.join(mdir,'step000825-cellfields-asagi.nc')
  out_ncpv   = os.path.join(mdir,'step000825-cellfields-paraview.nc')

  ptatin = ModelData()
  ptatin.ptatin_mesh  = ptatin.mesh_create_from_vts(vts)

  p2s = pTatin2SeisSol(
    ptatin,
    uniform_mesh_size=np.array([256,64,128],dtype=np.int32),
    total_stress_key='total_stress',
    deviatoric_stress_key='dev_stress',
    viscosity_key='viscosity',
    strainrate_key='StrainRate',
    velocity_key='u',
    pressure_key='pressure_p',
    density_key='density',
    temperature_key='T',
  )
  p2s.export_to_asagi(out_netcdf,fields_to_export=['total_stress','density','T'],pvoutput=out_ncpv)

  fault_mesh = pvs.read(fault)
  fault_mesh.points = p2s.rotation.rotate_referential(coor=fault_mesh.points,O=ptatin.O,L=ptatin.L)
  fault_mesh.point_data['Normals'] = p2s.rotation.rotate_vector(p2s.R,fault_mesh.point_data['Normals'])
  fault_mesh.save(os.path.join(mdir,'Main_SZ_surface_0_rot.vtp'))

  return

if __name__ == "__main__":
  test()