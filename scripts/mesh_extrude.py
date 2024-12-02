#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  GeoTEQpy
#  filename: mesh_extrude.py
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

import yaml
import argparse 
import pyvista as pvs
import geoteqpy as gte

def extrusion(
    vts_file:str,
    output:str,
    extrusion_params:dict,
    fields:list[str]|None=None,
    e2_key:str|None=None) -> None:
  
  mesh = pvs.read(vts_file)
  if e2_key is not None:
    gte.compute_xi(mesh,e2_key)
    if fields is not None:
      fields.append("xi")
    else:
      fields = ["xi"]

  for face in extrusion_params:
    extrusion = gte.Extrusion(mesh=mesh)
    mesh = extrusion.extrude_mesh(
      int(face["nsteps"]),
      float(face["dx"]),
      str(face["name"]),
      fields=fields
    )
  mesh.save(output)
  return

def main():

  description  = 'Extrude a mesh in the direction normal to the face(s) provided.\n'
  parser = argparse.ArgumentParser(prog='mesh_extrude.py',description=description,formatter_class=argparse.RawTextHelpFormatter)
  
  help_str = 'Path to the yaml file containing the info required to extrude the mesh.\n'
  help_str += 'The yaml file should contain the following keys:\n'
  help_str += '  model:\n'
  help_str += '    file: path to the vts file containing ptatin model data\n'
  help_str += '    output: path to the output vtu file\n'
  help_str += '    e2_key: (optional) key for the e2 field\n'
  help_str += '    fields: (optional) list of fields to extrude, if e2_key is provided at least \"xi\" field will be extruded\n'
  help_str += '  extrusion: \n'
  help_str += '  - name: name of the face to extrude\n'
  help_str += '    nsteps: number of steps to extrude\n'
  help_str += '    dx: step size for the extrusion\n'
  parser.add_argument('-f','--yaml_file',type=str,help=help_str,dest='yaml_file',required=True)

  args = parser.parse_args()
  yaml_file = args.yaml_file
  # read the yaml file
  with open(yaml_file,'r') as f:
    data = yaml.load(f,Loader=yaml.FullLoader)
  
  extrusion(
    data['model']['file'],
    data['model']['output'],
    data['extrusion'],
    fields=data['model'].get('fields',None),
    e2_key=data['model'].get('e2_key',None)
  )
  return

if __name__ == "__main__":
  main()
