#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  GeoTEQpy
#  filename: get_medial_axis.py
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
import yaml
import argparse
import pyvista as pvs
import geoteqpy as gte
from extract_contour import generate_contour_mesh

def compute_medial_axis(mesh:pvs.PolyData,radius_ma:float,get_eigv:bool,radius_cov:float) -> pvs.PolyData:
  ma = gte.MedialAxis(mesh=mesh,radius_ma=radius_ma,radius_cov=radius_cov)
  medial_axis = ma.get_medial_axis_mesh(get_eigv=get_eigv)
  return medial_axis

def main():
  description = "Compute the medial axis of a contour mesh.\n"
  parser = argparse.ArgumentParser(prog='get_medial_axis.py',description=description,formatter_class=argparse.RawTextHelpFormatter)
  help_str = 'Path to the yaml file containing the info required to compute the medial axis.\n'
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
  help_str += '  medial_axis:\n'
  help_str += '    radius_ma: (optional) radius of the medial axis. Default: 1e5\n'
  help_str += '    get_eigv_cov: (optional) if True, the eigenvectors and covariance are computed. Default: False\n'
  help_str += '    radius_cov: (optional) radius of the covariance. Default: 1e4\n'
  parser.add_argument('-f','--yaml_file',type=str,help=help_str,dest='yaml_file',required=True)

  args = parser.parse_args()
  yaml_file = args.yaml_file
  # read the yaml file
  with open(yaml_file,'r') as f:
    data = yaml.load(f,Loader=yaml.FullLoader)

  mesh = None
  if "model" not in data:
    contour_file = data.get("contour_file",None)
    if contour_file is None:
      raise ValueError("File containing the envelope for which the medial axis needs to be computed is missing. Use the key contour_file: \"path/to/file\".")
    file_basename = os.path.basename(contour_file)
    if "output" in data: output_dir = data["output"]
    else: output_dir = os.path.dirname(contour_file)
    # get the contour mesh
    contour_mesh = pvs.read(contour_file)
  else:
    # extract the file name without the path
    file_basename = os.path.basename(data['model']['file'])
    # output directory, use the one provided or the directory in which the model file is otherwise
    if "output" in data['model']: output_dir = data['model']['output']
    else:                         output_dir = os.path.dirname(data['model']['file'])
    # extract the contour mesh from the model
    contour_mesh = generate_contour_mesh(data,file_basename,output_dir)
    output = gte.replace_extension(os.path.join(output_dir,file_basename),"-contour.vtk")
    contour_mesh.save(output,recompute_normals=False)

  medial_axis = compute_medial_axis(
    contour_mesh,
    float(data["medial_axis"].get("radius_ma",1e5)),
    data["medial_axis"].get("get_eigv_cov",False),
    float(data["medial_axis"].get("radius_cov",1e4))
  )
  output = gte.replace_extension(os.path.join(output_dir,file_basename),"-ma.vtp")
  medial_axis.save(output)
  return

if __name__ == "__main__":
  main()
