import yaml
import geoteqpy as gte
import numpy as np
import argparse 

def tatin2asagi(
    uniform_mesh_size:np.ndarray, # size of the uniform mesh for asagi
    vts_file:str, # path to the vts file containing ptatin model data
    netcdf_file:str, # path to the output netcdf file
    fields_to_export:list[str], # list of fields to export into the netcdf file
    pvout:bool, # if True, a paraview readable netcdf file will be created
    **fields_keys # keys for the fields of ptatin model
  ):
  # ensure the extension is .nc
  if not netcdf_file.endswith('.nc'): netcdf_file += '.nc'
  ptatin = gte.ModelData(vts=vts_file)
  p2s = gte.pTatin2SeisSol(ptatin,uniform_mesh_size=uniform_mesh_size,**fields_keys)
  # if paraview output is required make a new file name
  if pvout: pv_file = netcdf_file.replace('.nc','-paraview.nc')
  p2s.export_to_asagi(netcdf_file,fields_to_export=fields_to_export,pvoutput=pv_file)
  return

def main():
  description  = '[1]: Extract data from pTatin3d model result stored in a .vts file.\n'
  description += '[2]: Create a uniform grid to interpolate pTatin3d data.\n'
  description += '[3]: Export the interpolated fields to a netcdf file compatible with ASAGI for SeisSol.'
  parser = argparse.ArgumentParser(prog='ptatin2asagi.py',description=description,formatter_class=argparse.RawTextHelpFormatter)
  
  help_str = 'Path to the yaml file containing the info required to export from ptatin to netcdf.\n'
  help_str += 'The yaml file should contain the following keys:\n'
  help_str += '  model:\n'
  help_str += '    file: path to the vts file containing ptatin model data\n'
  help_str += '    fields_keys: keys for the fields of ptatin model\n'
  help_str += '  asagi:\n'
  help_str += '    file: path to the output netcdf file\n'
  help_str += '    fields_keys: list of fields to export into the netcdf file\n'
  help_str += '    mesh_size: size of the uniform mesh for asagi\n'
  help_str += '    paraview: (optionnal) if True, a paraview readable netcdf file will be created\n'
  parser.add_argument('-f','--yaml_file',type=str,help=help_str,dest='yaml_file',required=True)

  args = parser.parse_args()
  yaml_file = args.yaml_file
  # read the yaml file
  with open(yaml_file,'r') as f:
    data = yaml.load(f,Loader=yaml.FullLoader)
  
  ptatin:dict = data['model']
  asagi:dict  = data['asagi']

  mandatory_keys = ['file','fields_keys']
  for key in mandatory_keys:
    if key not in asagi:
      raise ValueError(f"Parameter {key} is missing from the 'asagi' section in the provided file: {yaml_file}.")
    if key not in ptatin:
      raise ValueError(f"Parameter {key} is missing from the 'model' section in the provided file: {yaml_file}.")
  if 'mesh_size' not in asagi:
    raise ValueError(f"Parameter mesh_size is missing from the 'asagi' section in the provided file: {yaml_file}.")

  tatin2asagi(
    np.asarray(asagi['mesh_size'], dtype=np.int32), # size of the uniform mesh for asagi
    ptatin['file'], # path to the vts file containing ptatin model data
    asagi['file'], # path to the output netcdf file
    asagi['fields_keys'], # list of fields to export into the netcdf file
    asagi.get('paraview', False), # if True, a paraview readable netcdf file will be created
    **ptatin['fields_keys'] # keys for the fields of ptatin model
  )

if __name__ == "__main__":
  main()
