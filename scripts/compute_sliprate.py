import os
import numpy as np
import pyvista as pvs
import pyviztools as pvt
import argparse
import yaml

from mesh_extrude import extrusion

def compute_sliprate(fault_mesh:pvs.PolyData,fe_mesh:pvt.FEMeshQ1,velocity:np.ndarray,dx:float):
  # shift fault mesh by dx in its normal direction
  f_plus  = fault_mesh.points + dx * fault_mesh.point_data['Normals']
  f_minus = fault_mesh.points - dx * fault_mesh.point_data['Normals']

  # locate shifted fault mesh in the FE mesh
  el_plus,xi_plus   = pvt.utils.femesh_point_location(fe_mesh, np.float32(f_plus))
  el_minus,xi_minus = pvt.utils.femesh_point_location(fe_mesh, np.float32(f_minus))

  # interpolate velocity on shifted fault mesh
  v_plus  = fe_mesh.interp_vertex_field(velocity,el_plus,xi_plus)
  v_minus = fe_mesh.interp_vertex_field(velocity,el_minus,xi_minus)

  # compute slip rate
  sliprate = v_plus - v_minus
  fault_mesh.point_data['sliprate'] = sliprate
  return

def process_faults(ptatin:pvs.StructuredGrid,fault_dir:str,faults:list[str],dx:float,velocity_key:str):
  # create Q1 FE mesh
  mesh = pvt.CFEMeshQ1(3,ptatin.dimensions[0]-1,ptatin.dimensions[1]-1,ptatin.dimensions[2]-1)
  # create the element to vertex connectivity
  mesh.create_e2v()
  # attach coordinates to the mesh
  mesh.coor = ptatin.points

  for fault in faults:
    fname = os.path.join(fault_dir,fault)
    fault_mesh = pvs.read(fname)
    compute_sliprate(fault_mesh,mesh,ptatin.point_data[velocity_key],dx)
    fault_mesh.save(fname)
  return


def main():
  description = 'Compute slip rate on 3D fault meshes.\n'
  description += 'Fault meshes must contain normal vectors.\n'
  parser = argparse.ArgumentParser(
    prog='compute_sliprate.py',
    description=description,
    formatter_class=argparse.RawTextHelpFormatter
  )
  helptxt = 'Path to the yaml configuration file.\n'
  helptxt += 'The configuration file must contain the following keys:\n'
  helptxt += '  model_fname: Path to the model file.\n'
  helptxt += '  fault_dir: Path to the directory containing the fault meshes.\n'
  helptxt += '  faults_fname: List of fault mesh filenames.\n'
  helptxt += '  dx: Distance to shift the fault mesh in the normal direction.\n'
  helptxt += '  velocity_key: Key to access the velocity field in the model file.\n'
  parser.add_argument('-f', '--file', type=str, dest='file', help=helptxt, required=True)

  args = parser.parse_args()
  with open(args.file,'r') as f:
    config:dict = yaml.load(f,Loader=yaml.FullLoader)

  if config.get('model_fname', None) is None:
    serr = f'model_fname not found in the provided configuration file: {args.file}\n'
    serr += 'Use -h option to see the help message.\n'
    raise ValueError(serr)
  
  if config.get('extrusion', None) is not None:
    print("Extruding mesh...")
    mesh_input = pvs.read(config['model_fname'])
    # extrude mesh
    mesh = extrusion(
      mesh_input,
      config['extrusion'],
      fields=config.get('fields',None),
      e2_key=config.get('e2_key',None)
    )
  else:
    # Load the mesh
    mesh = pvs.read(config['model_fname'])

  process_faults(
    mesh,
    config['fault_dir'],
    config['faults_fname'],
    float(config['dx']),
    config.get('velocity_key','velocity')
  )
  return


if __name__ == '__main__':
  main()