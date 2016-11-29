from __future__ import print_function
# Copyright (C) 2008 CSC - Scientific Computing Ltd.
"""
This module defines an ASE interface to ORCA.

"""

import os
import sys
import re
import warnings
from .general import Calculator
from os.path import join, isfile, islink

import numpy as np

import ase
import ase.io
from ase.utils import devnull

from ase.calculators.singlepoint import SinglePointCalculator


block_keys = [
    'BASIS',      # Basis sets are specified
    'CASSCF',     # Control of CASSCF/NEVPT2 and DMRG calculations
    'CIS',        # Control of CIS and TD-DFT calculations (synonym is TDDFT)
    'COORDS',     # Input of atomic coordinates
    'COSMO',      # Control of the conductor like screening model
    'ELPROP',     # Control of electric property calculations
    'EPRNMR',     # Control of SCF level EPR and NMR calculations
    'FREQ',       # Control of frequency calculations
    'GEOM',       # Control of geometry optimization
    'MD',         # Control of molecular dynamics simulation
    'LOC',        # Localization of orbitals
    'MDCI',       # Controls single reference correlation methods
    'METHOD',     # Here a computation method is specified
    'MP2',        # Controls the details of the MP2 calculation
    'MRCI',       # Control of MRCI calculations
    'OUTPUT',     # Control of output
    'PAL',        # Control of parallel jobs
    'PARAS',      # Input of semiempirial parameters
    'PLOTS',      # Control of plot generation
    'REL',        # Control of relativistic options
    'RR',         # Control of resonance Raman and abs/fluor bandshape calcs
    'SCF'         # Control of the SCF procedure
]

class Orca(Calculator):
    name = 'Orca'

    # Parameters corresponding to 'xc' settings.  This may be modified
    # by the user in-between loading calculators.vasp submodule and
    # instantiating the calculator object with calculators.vasp.Vasp()
    default_params = {
        }

    def __init__(self, restart=None,
                 **kwargs):
        self.block_params = {}
        self.params = {}
        for key in block_keys:
            self.block_params[key] = {}
        self.atoms = None
        self.positions = None
        self.run_counts = 0
        self.set(**kwargs)

    def set(self, **kwargs):
        for key in kwargs:
            self.params[key] = kwargs[key]
        self.set_block_params(**kwargs)
    
    def set_block_params(self, **kwargs):
        for key in kwargs:
            splitkey = key.split('_')
            self.block_params[splitkey[0]][splitkey[1]] = kwargs[key]

    def initialize(self, atoms):
        """
        Initialize an ORCA calculation

        """

        p = self.input_params

        self.all_symbols = atoms.get_chemical_symbols()
        self.natoms = len(atoms)
        self.spinpol = atoms.get_initial_magnetic_moments().any()
        atomtypes = atoms.get_chemical_symbols()

        # Determine the number of atoms of each atomic species
        # sorted after atomic species
        symbols = []
        symbolcount = {}

        for m, atom in enumerate(atoms):
            symbol = atom.symbol
            if m in special_setups:
                pass
            else:
                if symbol not in symbols:
                    symbols.append(symbol)
                    symbolcount[symbol] = 1
                else:
                    symbolcount[symbol] += 1


        self.converged = None
        self.setups_changed = None

    def write_inp(self, atoms, **kwargs):
        inp = open('.inp','w')

        for key, val in self.block_params.items():
            if val:
                inp.write('    %%%s\n' % (key.upper()))
                for innerKey, innerVal in val.items():
                    if val is not None: # need to go through types of vals: int,float,etc
                        pass
                inp.write(' end\n\n')


        
    def print_atoms(self, atoms):
        out=''
        for atom in atoms:
            atom_str= '{0}\t{1}\t{2}\t{3}\n'.format(atom.symbol, atom.position[0], atom.position[1], atom.position[2])
            out=out+atom_str
        return out

    def calculate(self, atoms):
        """Generate necessary files in the working directory and run ORCA.

        The method first write ORCA input files, then calls the method
        which executes ORCA. When the ORCA run is finished energy, forces,
        etc. are read from the ORCA output.
        """

        # Initialize calculations
        self.initialize(atoms)

        # Write input

        # Execute 
        self.run()
        # Read output

    def set_results(self, atoms):
        pass

    def run(self):
        """Method which explicitely runs ORCA."""
        pass
    
    def clean(self):
        pass

    def set_atoms(self, atoms):
        if (atoms != self.atoms):
            self.converged = None
        self.atoms = atoms.copy()

    def get_atoms(self):
        atoms = self.atoms.copy()
        atoms.set_calculator(self)
        return atoms

from ase.atoms import Atoms
from ase.io.extxyz import read_extxyz as read_xyz, write_extxyz as write_xyz

__all__ = ['read_xyz', 'write_xyz']

def read_orca(fileobj, index):
    lines = fileobj.readlines()
    natoms = int(lines[0])
    nimages = len(lines) // (natoms + 2)
    for i in range(*index.indices(nimages)):
        symbols = []
        positions = []
        n = i * (natoms + 2) + 2
        for line in lines[n:n + natoms]:
            symbol, x, y, z = line.split()[:4]
            symbol = symbol.lower().capitalize()
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
        yield Atoms(symbols=symbols, positions=positions)


def write_orca(fileobj, images, comment=''):
    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    for atoms in images:
        fileobj.write('%d\n%s\n' % (natoms, comment))
        for s, (x, y, z) in zip(symbols, atoms.positions):
            fileobj.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))

