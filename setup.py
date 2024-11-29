import os
from setuptools import setup, find_packages, Extension

# Generate the list of source files automatically
src_dir = 'geoteqpy/c/src'
sources = [os.path.join(src_dir, f) for f in os.listdir(src_dir) if f.endswith('.c')]

# Define the extension module
extension = Extension(
    name='geoteqpy.c.lib.faulttools',  # Name of the extension
    sources=sources, # Source files  
    include_dirs=['geoteqpy/c/src'],  # Include directories
    libraries=['m'],  # Libraries to link against
    library_dirs=[],  # Directories to search for libraries
    extra_compile_args=['-march=native'],  # Extra arguments to pass to the compiler
    extra_link_args=['-shared', '-fPIC']  # Extra arguments to pass to the linker
)

setup(
    name='geoteqpy',
    version='1.0.0',
    ext_modules=[extension],
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