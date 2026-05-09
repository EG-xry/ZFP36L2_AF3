#!/usr/bin/env python
"""Test PlumedForce integration"""
import os
os.environ['DYLD_LIBRARY_PATH'] = '/Users/eric/micromamba/envs/omm/lib'

from openmm.app import *
from openmm import *
from openmmplumed import PlumedForce
from mpi4py import MPI

print("Testing PlumedForce integration...")

# Load system
pdb = PDBFile('input_af.pdb')
top = pdb.topology
forcefield = ForceField('forcefield.xml')
system = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic)
system.removeForce(3)

# Create PlumedForce
comm1 = MPI.COMM_SELF
comm2 = MPI.COMM_WORLD

with open('plumed.dat', 'r') as f:
    script = f.read()

print("Creating PlumedForce...")
plumed_force = PlumedForce(script, comm1, comm2)
print("PlumedForce created successfully!")

print("Adding PlumedForce to system...")
system.addForce(plumed_force)
print("SUCCESS: PlumedForce added to system without errors!")
print("Python wrapper rebuild was successful!")

