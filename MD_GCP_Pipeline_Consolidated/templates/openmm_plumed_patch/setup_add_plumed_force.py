from setuptools import setup, Extension
import os
import sys

# Get paths from environment or use defaults
openmm_dir = os.environ.get('OPENMM_DIR', '/Users/eric/micromamba/envs/omm')
plumed_dir = os.environ.get('PLUMED_DIR', '/Users/eric/miniconda3/envs/afmi-tutorial')

# Find include directories
openmm_include = os.path.join(openmm_dir, 'include')
plumed_include = os.path.join(plumed_dir, 'include')

# Find library directories
openmm_lib = os.path.join(openmm_dir, 'lib')
plumed_lib = os.path.join(plumed_dir, 'lib')

# Find the PlumedForce header
plumed_force_include = os.path.join(os.path.dirname(__file__), '..', 'openmmapi', 'include')

module = Extension(
    'add_plumed_force',
    sources=['add_plumed_force.cpp'],
    include_dirs=[
        openmm_include,
        plumed_include,
        plumed_force_include,
        '/Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/OpenMM-Plumed-MPI/openmmapi/include',
    ],
    library_dirs=[
        openmm_lib,
        plumed_lib,
    ],
    libraries=['OpenMM', 'plumed'],
    language='c++',
    extra_compile_args=['-std=c++11'],
)

setup(
    name='add_plumed_force',
    version='1.0',
    description='C extension to add PlumedForce to System',
    ext_modules=[module],
)

