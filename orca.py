from __future__ import print_function
# Copyright (C) 2008 CSC - Scientific Computing Ltd.
"""This module defines an ASE interface to VASP.

Developed on the basis of modules by Jussi Enkovaara and John
Kitchin.  The path of the directory containing the pseudopotential
directories (potpaw,potpaw_GGA, potpaw_PBE, ...) should be set
by the environmental flag $VASP_PP_PATH.

The user should also set the environmental flag $VASP_SCRIPT pointing
to a python script looking something like::

   import os
   exitcode = os.system('vasp')


http://cms.mpi.univie.ac.at/vasp/
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

# Parameters that can be set in inp file. The values which are None
# are not written and default parameters of ORCA are used for them.

float_keys = [
]

exp_keys = [
]

string_keys = [
    'METHOD',
    'RUNTYP',
    'CITYPP',
    'FROZENCORE',
    'HFTYP',
    'ALLOWRHF',
    'RI',
    'KMATRIX',
    'JMATRIX',
    'SCFMODE',
    'MAXITER',
]

int_keys = [
]

bool_keys = [
]

list_keys = [
]

dict_keys = [
]

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
    'RR',         # Control of resonance Raman and absorption/fluorescence bandshape calculations
    'SCF'         # Control of the SCF procedure
]

class Vasp(Calculator):
    name = 'Vasp'

    # Parameters corresponding to 'xc' settings.  This may be modified
    # by the user in-between loading calculators.vasp submodule and
    # instantiating the calculator object with calculators.vasp.Vasp()
    defaults = {
        }

    def __init__(self, restart=None,
                 output_template='vasp',
                 track_output=False,
                 **kwargs):
        self.float_params = {}
        self.exp_params = {}
        self.string_params = {}
        self.int_params = {}
        self.bool_params = {}
        self.list_params = {}
        for key in float_keys:
            self.float_params[key] = None
        for key in exp_keys:
            self.exp_params[key] = None
        for key in string_keys:
            self.string_params[key] = None
        for key in int_keys:
            self.int_params[key] = None
        for key in bool_keys:
            self.bool_params[key] = None
        for key in list_keys:
            self.list_params[key] = None
        for key in special_keys:
            self.special_params[key] = None
        for key in dict_keys:
            self.dict_params[key] = None

        # Initialize internal dictionary of input parameters which are
        # not regular VASP keys
        self.input_params = {
            'xc': None,  # Exchange-correlation recipe (e.g. 'B3LYP')
            'pp': None,  # Pseudopotential file (e.g. 'PW91')
            'setups': None,  # Special setups (e.g pv, sv, ...)
            'txt': '-',  # Where to send information
            'kpts': (1, 1, 1),  # k-points
            # Option to use gamma-sampling instead of Monkhorst-Pack:
            'gamma': False,
            # number of points between points in band structures:
            'kpts_nintersections': None,
            # Option to write explicit k-points in units
            # of reciprocal lattice vectors:
            'reciprocal': False}

        self.restart = restart
        self.track_output = track_output
        self.output_template = output_template
        if restart:
            self.restart_load()
            return


        self.nbands = self.int_params['nbands']
        self.atoms = None
        self.positions = None
        self.run_counts = 0
        self.set(**kwargs)

    def set(self, **kwargs):
        if 'xc' in kwargs:
            self.set_xc_params(kwargs['xc'])
        for key in kwargs:
            if key in self.float_params:
                self.float_params[key] = kwargs[key]
            elif key in self.exp_params:
                self.exp_params[key] = kwargs[key]
            elif key in self.string_params:
                self.string_params[key] = kwargs[key]
            elif key in self.int_params:
                self.int_params[key] = kwargs[key]
            elif key in self.bool_params:
                self.bool_params[key] = kwargs[key]
            elif key in self.list_params:
                self.list_params[key] = kwargs[key]
            elif key in self.special_params:
                self.special_params[key] = kwargs[key]
            elif key in self.dict_params:
                self.dict_params[key] = kwargs[key]
            elif key in self.input_params:
                self.input_params[key] = kwargs[key]
            else:
                raise TypeError('Parameter not defined: ' + key)

    def update(self, atoms):
        if self.calculation_required(atoms, ['energy']):
            if (((self.atoms is None) or
                 (self.atoms.positions.shape != atoms.positions.shape)
                 )):
                # Completely new calculation just reusing the same
                # calculator, so delete any old VASP files found.
                self.clean()
            self.calculate(atoms)

    def initialize(self, atoms):
        """Initialize a VASP calculation

        """

        p = self.input_params

        self.all_symbols = atoms.get_chemical_symbols()
        self.natoms = len(atoms)
        self.spinpol = atoms.get_initial_magnetic_moments().any()
        atomtypes = atoms.get_chemical_symbols()

        # Determine the number of atoms of each atomic species
        # sorted after atomic species
        special_setups = []
        symbols = []
        symbolcount = {}
        if self.input_params['setups']:
            for m in self.input_params['setups']:
                try:
                    special_setups.append(int(m))
                except ValueError:
                    continue

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

        # Build the sorting list
        self.sort = []
        self.sort.extend(special_setups)

        for symbol in symbols:
            for m, atom in enumerate(atoms):
                if m in special_setups:
                    pass
                else:
                    if atom.symbol == symbol:
                        self.sort.append(m)
        self.resort = list(range(len(self.sort)))
        for n in range(len(self.resort)):
            self.resort[self.sort[n]] = n
        self.atoms_sorted = atoms[self.sort]

        # Check if the necessary POTCAR files exists and
        # create a list of their paths.
        self.symbol_count = []
        for m in special_setups:
            self.symbol_count.append([atomtypes[m], 1])
        for m in symbols:
            self.symbol_count.append([m, symbolcount[m]])

        sys.stdout.flush()

        # Potpaw folders may be identified by an alias or full name
        for pp_alias, pp_folder in (('lda', 'potpaw'),
                                    ('pw91', 'potpaw_GGA'),
                                    ('pbe', 'potpaw_PBE')):
            if p['pp'].lower() == pp_alias:
                break
        else:
            pp_folder = p['pp']

        if 'VASP_PP_PATH' in os.environ:
            pppaths = os.environ['VASP_PP_PATH'].split(':')
        else:
            pppaths = []
        self.ppp_list = []
        # Setting the pseudopotentials, first special setups and
        # then according to symbols
        for m in special_setups:
            if m in p['setups']:
                special_setup_index = m
            elif str(m) in p['setups']:
                special_setup_index = str(m)
            else:
                raise Exception("Having trouble with special setup index {0}."
                                " Please use an int.".format(m))
            potcar = join(pp_folder,
                          p['setups'][special_setup_index],
                          'POTCAR')
            for path in pppaths:
                filename = join(path, potcar)

                if isfile(filename) or islink(filename):
                    self.ppp_list.append(filename)
                    break
                elif isfile(filename + '.Z') or islink(filename + '.Z'):
                    self.ppp_list.append(filename + '.Z')
                    break
            else:
                print('Looking for %s' % potcar)
                raise RuntimeError('No pseudopotential for %s!' % symbol)

        for symbol in symbols:
            try:
                potcar = join(pp_folder, symbol + p['setups'][symbol],
                              'POTCAR')
            except (TypeError, KeyError):
                potcar = join(pp_folder, symbol, 'POTCAR')
            for path in pppaths:
                filename = join(path, potcar)

                if isfile(filename) or islink(filename):
                    self.ppp_list.append(filename)
                    break
                elif isfile(filename + '.Z') or islink(filename + '.Z'):
                    self.ppp_list.append(filename + '.Z')
                    break
            else:
                print('''Looking for %s
                The pseudopotentials are expected to be in:
                LDA:  $VASP_PP_PATH/potpaw/
                PBE:  $VASP_PP_PATH/potpaw_PBE/
                PW91: $VASP_PP_PATH/potpaw_GGA/''' % potcar)
                raise RuntimeError('No pseudopotential for %s!' % symbol)
        self.converged = None
        self.setups_changed = None

    def calculate(self, atoms):
        """Generate necessary files in the working directory and run ORCA.

        The method first write ORCA input files, then calls the method
        which executes ORCA. When the ORCA run is finished energy, forces,
        etc. are read from the ORCA output.
        """

        # Initialize calculations
        self.initialize(atoms)

        # Write input

        # Execute VASP
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

