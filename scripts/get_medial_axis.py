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

import yaml
import argparse
import pyvista as pvs
import geoteqpy as gte

def compute_medial_axis(contour_file:str,radius_ma:float,radius_cov:float,output:str|None=None) -> None:
  mesh = pvs.read(contour_file)
  ma = gte.MedialAxis(mesh=mesh,radius_ma=radius_ma,radius_cov=radius_cov)
  medial_axis = ma.get_medial_axis_mesh()
  if output is None:
    output = gte.replace_extension(contour_file,"_medial_axis.vtp")
  medial_axis.save(output)
  return

def main():
  description = "Compute the medial axis of a contour mesh.\n"
  parser = argparse.ArgumentParser(prog='get_medial_axis.py',description=description,formatter_class=argparse.RawTextHelpFormatter)
  help_str = 'Path to the yaml file containing the info required to compute the medial axis.\n'
  help_str += 'The yaml file should contain the following keys:\n'
  help_str += '  contour_file: path to the contour mesh file\n'
  help_str += '  radius_ma: radius of research for the medial axis\n'
  help_str += '  radius_cov: radius of the covariance matrix\n'
  help_str += '  output: (optional) path to the output medial axis file\n'
  parser.add_argument('-f','--yaml_file',type=str,help=help_str,dest='yaml_file',required=True)

  args = parser.parse_args()
  yaml_file = args.yaml_file
  # read the yaml file
  with open(yaml_file,'r') as f:
    data = yaml.load(f,Loader=yaml.FullLoader)

  compute_medial_axis(
    data["contour_file"],
    float(data["radius_ma"]),
    float(data["radius_cov"]),
    data.get("output",None)
  )
  return

if __name__ == "__main__":
  main()
