import os
import yaml
import argparse 
import pyvista as pvs
import geoteqpy as gte
from mesh_extrude import extrusion

def get_contour(mesh:pvs.StructuredGrid,isovalue:float=20.0,field_name:str="xi",flip_normals:bool=False) -> pvs.PolyData:
  point_mesh = mesh.cell_data_to_point_data(pass_cell_data=True)
  contour_mesh = point_mesh.contour(isosurfaces=[isovalue],scalars=field_name,compute_normals=False)
  contour_mesh.compute_normals(flip_normals=flip_normals,inplace=True)
  return contour_mesh

def generate_contour_mesh(data:dict, file_basename:str, output_dir:str) -> pvs.PolyData:
  mesh = pvs.read(data['model']['file'])
  if "e2_key" not in data['model']:
    print("The key for the strain-rate norm field is not provided, trying the default value: \"e2\".")
    data['model']['e2_key'] = "e2"
  if data['model']['e2_key'] not in mesh.array_names:
    raise ValueError(f"Key for the strain-rate norm field not found. Given key: {data['model']['e2_key']}")

  fields = data['model'].get('fields',None)
  # extrude mesh if the extrusion key is present
  if "extrusion" in data:
    print("Extruding mesh...")
    # extrude
    mesh = extrusion(
      mesh,
      data['extrusion'],
      fields=fields,
      e2_key=data['model']['e2_key']
    )
    output = gte.replace_extension(os.path.join(output_dir,file_basename),"-extruded.vts")
    mesh.save(output)
  else:
    # If no extrusion, load the mesh and compute the xi field
    mesh = pvs.read(data['model']['file'])
    gte.compute_xi(mesh,field_name="xi",e2_key=data['model']['e2_key'])
  
  if 'contour' not in data:
    data['contour'] = {}
  contour_mesh = get_contour(
    mesh,
    isovalue=data['contour'].get('isovalue',0.75),
    field_name=data['contour'].get('field_name',"xi"),
    flip_normals=data['contour'].get('flip_normals',False)
  )
  
  # add normals to the fields to output
  if fields is None:
    fields = ["Normals"]
  else:
    fields.append("Normals")

  # remove all fields except the ones in the list, reduces the size of the output file
  for f in contour_mesh.array_names:
    # remove all cell data, it does not make sense for the contour mesh
    if f in contour_mesh.cell_data:
      contour_mesh.cell_data.remove(f)
    if f not in fields:
      if f in contour_mesh.point_data:
        contour_mesh.point_data.remove(f)
  return contour_mesh

def main():
  description  = 'Extracts the isosurface of the required field.\n'
  description += 'If required, an extrusion of the mesh is performed in the direction normal to the face(s) provided.\n'
  parser = argparse.ArgumentParser(prog='mesh_extrude.py',description=description,formatter_class=argparse.RawTextHelpFormatter)
  
  help_str = 'Path to the yaml file containing the info required to extrude the mesh.\n'
  help_str += 'The yaml file should contain the following keys:\n'
  help_str += '  model:\n'
  help_str += '    file: path to the vts file containing ptatin model data\n'
  help_str += '    output: path to the output directory\n'
  help_str += '    e2_key: (optional) key for the e2 field. Default: \"e2\"\n'
  help_str += '    fields: (optional) list of fields to extrude, if e2_key is provided at least \"xi\" field will be extruded\n'
  help_str += '  contour: (optional). If not provided all default values listed below are used.\n'
  help_str += '    isovalue: (optional) value of the isosurface to extract the contour. Default: 0.75\n'
  help_str += '    field_name: (optional) name of the field to extract the contour. Default: \"xi\"\n'
  help_str += '  extrusion: (optional). If provided, extrudes the mesh and saves it as a .vts file.\n'
  help_str += '  - name: name of the face to extrude\n'
  help_str += '    nsteps: number of steps to extrude\n'
  help_str += '    dx: step size for the extrusion\n'
  parser.add_argument('-f','--yaml_file',type=str,help=help_str,dest='yaml_file',required=True)

  args = parser.parse_args()
  yaml_file = args.yaml_file
  # read the yaml file
  with open(yaml_file,'r') as f:
    data = yaml.load(f,Loader=yaml.FullLoader)
  
  # extract the file name without the path
  file_basename = os.path.basename(data['model']['file'])
  # output directory, use the one provided or the directory in which the model file is otherwise
  if "output" in data['model']: output_dir = data['model']['output']
  else:                         output_dir = os.path.dirname(data['model']['file'])

  contour_mesh = generate_contour_mesh(data,file_basename,output_dir)

  output = gte.replace_extension(os.path.join(output_dir,file_basename),"-contour.vtk")
  contour_mesh.save(output,recompute_normals=False)
  return

if __name__ == "__main__":
  main()