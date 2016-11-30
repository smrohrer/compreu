from __future__ import print_function
# Copyright (C) 2008 CSC - Scientific Computing Ltd.
"""
This module defines an ASE interface to ORCA.

"""

import os
import sys
import re
import warnings
from ase.calculators.general import Calculator
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

    default_params = {
        'METHOD_METHOD' :   'HF',
        'METHOD_RUNTYP' :   'ENERGY',
        'BASIS_BASIS'   :   'STO_3G',
        'COORDS_CHARGE' : 0,
        'COORDS_MULT' : 1,
        'COORDS_CTYP' : 'XYZ'
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
        self.set(**Orca.default_params)
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
            if symbol not in symbols:
                symbols.append(symbol)
                symbolcount[symbol] = 1
            else:
                symbolcount[symbol] += 1
        self.converged = None
        self.setups_changed = None

    def write_inp(self, atoms, projectName="orca", **kwargs):
        inp = open('%s.inp' % (projectName),'w')
        for key, val in self.block_params.items():
            if val:
                inp.write('%%%s\n' % (key.upper()))
                for innerKey, innerVal in val.items():
                    self.write_line(inp, innerKey, innerVal)
                if key == 'COORDS':
                    self.write_atoms(inp, atoms)
                inp.write('end\n\n')
        inp.close()

    def write_atoms(self, inp, atoms):
        inp.write('    COORDS\n')
        for atom in atoms:
            s = atom.symbol
            [x, y, z] = atom.position
            inp.write('        %s\t%5.6f\t%5.6f\t%5.6f\n' % (s,x,y,z))
        inp.write('    end\n')

    def write_line(self, open_file, key, val):
        # Line non-block line with these possible types:
        #   int
        #   float
        #   exponential form (not sure how to do this one, 
        #       currenty just a redundent elif)
        #   string
        #   Bool: should just be having the keyword or not
        if type(val) == int:
            open_file.write('    %s %d\n' % (key, val))
        elif type(val) == float:
            open_file.write('    %s %5.6f\n' % (key, val))
        elif type(val) == float:
            open_file.write('    %s %5.2e\n' % (key, val))
        elif type(val) == str:
            open_file.write('    %s %s\n' % (key, val))
        elif type(val) == bool:
            open_file.write('    %s\n' % (key))
        else:
            pass

    def calculate(self, atoms):
        """
        Generate necessary files in the working directory and run ORCA.
        """

        # Initialize calculations
        self.initialize(atoms)

        # Write input
        self.write_inp(atoms)

        # Execute 
        self.run()
        # Read output

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
