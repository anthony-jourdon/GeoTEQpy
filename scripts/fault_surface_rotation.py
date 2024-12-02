#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  GeoTEQpy
#  filename: fault_surface_rotation.py
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
import geoteqpy as gte
import pyvista as pvs
import argparse 

def fault_rotation(vts_file:str,fault_files:list[str]):
  ptatin = gte.ModelData(vts=vts_file)
  p2s    = gte.pTatin2SeisSol(ptatin,None)

  for fault_file in fault_files:
    fault_mesh = pvs.read(fault_file)
    fault_mesh.points = p2s.rotation.rotate_referential(coor=fault_mesh.points,O=ptatin.O,L=ptatin.L)
    if 'Normals' in fault_mesh.point_data:
      fault_mesh.point_data['Normals'] = p2s.rotation.rotate_vector(p2s.R,fault_mesh.point_data['Normals'])
    f,ext = os.path.splitext(fault_file)
    outfname = f + '_rot' + ext
    fault_mesh.save(outfname)
  return

def main():
  description  = 'Rotate faults surface from ptatin3d coord system xyz (y vertical), to xzy (z vertical).'
  parser = argparse.ArgumentParser(prog='fault_surface_rotation.py',description=description,formatter_class=argparse.RawTextHelpFormatter)
  
  help_str = 'Path to the yaml file containing the info required to rotate fault(s).\n'
  help_str += 'The yaml file should contain the following keys:\n'
  help_str += '  model:\n'
  help_str += '    file: path to the vts file containing ptatin model data\n'
  help_str += '  faults:\n'
  help_str += '    file: list of paths to the fault files to rotate\n'
  parser.add_argument('-f','--yaml_file',type=str,help=help_str,dest='yaml_file',required=True)

  args = parser.parse_args()
  yaml_file = args.yaml_file
  # read the yaml file
  with open(yaml_file,'r') as f:
    data = yaml.load(f,Loader=yaml.FullLoader)
  
  fault_rotation(data['model']['file'],data['faults']['file'])
  
  return

if __name__ == "__main__":
  main()
