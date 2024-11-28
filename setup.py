from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import subprocess
import os
import shutil

class CustomBuildExtCommand(build_ext):
    """Custom build_ext command to run `make all`."""
    def run(self):
        # Run the original build_ext command
        build_ext.run(self)
        # Run `make all` to compile the C files
        subprocess.check_call(['make', 'all'], cwd='geoteqpy/c')
        # Verify that the shared library was created
        source_path = os.path.join('geoteqpy', 'c', 'lib', 'faulttools.so')
        if not os.path.exists(source_path):
            raise FileNotFoundError(f"Expected shared library not found: {source_path}")
        # Move the compiled shared library to the desired location
        build_lib = self.build_lib
        target_dir = os.path.join(build_lib, 'geoteqpy', 'c', 'lib')
        os.makedirs(target_dir, exist_ok=True)
        destination_path = os.path.join(target_dir, 'faulttools.so')
        if os.path.exists(destination_path):
            os.remove(destination_path)
        shutil.move(source_path, destination_path)
        print(f"Moved {source_path} to {destination_path}")

# Generate the list of source files automatically
src_dir = 'geoteqpy/c/src'
sources = [os.path.join(src_dir, f) for f in os.listdir(src_dir) if f.endswith('.c')]

# Define the extension module
extension = Extension(
    name='geoteqpy.c.faulttools',  # Name of the extension
    sources=sources, # Source files  
    include_dirs=[],  # Include directories
    libraries=[],  # Libraries to link against
    library_dirs=[],  # Directories to search for libraries
    extra_compile_args=[],  # Extra arguments to pass to the compiler
    extra_link_args=[]  # Extra arguments to pass to the linker
)

setup(
    name='geoteqpy',
    version='1.0.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        # List your package dependencies here
        "numpy",
        "pyvista",
        "netCDF4",
        "ptatin3d-pyviztools @ git+https://bitbucket.org/ptatin/ptatin3d-pyviztools.git@anthony_jourdon/post-proc-script",
    ],
    entry_points={
        'console_scripts': [
            # Define command-line scripts here
        ],
    },
    cmdclass={
        'build_ext': CustomBuildExtCommand,
    },
    ext_modules=[extension],
    author='Anthony Jourdon',
    author_email='',
    description='Python package to extrude structured mesh, compute the medial axis of a 3D shape and export data from long-term geodynamic models to ASAGI using NetCDF.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/anthony-jourdon/GeoTEQpy',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12',
)