#get eigenvalues in ev, fillings (for homo/lumo), total energy, orbital energy

"""https://gitlab.com/ase/ase/blob/master/ase/calculators/dftb.py
http://www.dftb-plus.info/fileadmin/DFTB-Plus/public/recipes/html/basics/firstcalc.html"""


import os
import numpy as np
from ase.calculators.calculator import FileIOCalculator, kpts2mp
class Dftb(FileIOCalculator):
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='dftb', atoms=None, kpts=None,
                 run_manyDftb_steps=False,
                 **kwargs):
        # I stole this from a site and didn't touch their init function
        """Construct a DFTB+ calculator.
        run_manyDftb_steps:  Logical
            True: many steps are run by DFTB+,
            False:a single force&energy calculation at given positions
        ---------
        Additional object (to be set by function embed)
        pcpot: PointCharge object
            An external point charge potential (only in qmmm)
        """
        from ase.dft.kpoints import monkhorst_pack
        if 'DFTB_PREFIX' in os.environ:
            slako_dir = os.environ['DFTB_PREFIX']
        else:
            slako_dir = './'
        # to run Dftb as energy and force calculator use
        # Driver_MaxSteps=0,
        if run_manyDftb_steps:
            # minimisation of molecular dynamics is run by native DFTB+
            self.default_parameters = dict(
                Hamiltonian_='DFTB',
                Hamiltonian_SlaterKosterFiles_='Type2FileNames',
                Hamiltonian_SlaterKosterFiles_Prefix=slako_dir,
                Hamiltonian_SlaterKosterFiles_Separator='"-"',
                Hamiltonian_SlaterKosterFiles_Suffix='".skf"',
                Hamiltonian_MaxAngularMomentum_='')
        else:
            # using ase to get forces and energy only
            # (single point calculation)
            self.default_parameters = dict(
                Hamiltonian_='DFTB',
                Driver_='ConjugateGradient',
                Driver_MaxForceComponent='1E-4',
                Driver_MaxSteps=0,
                Hamiltonian_SlaterKosterFiles_='Type2FileNames',
                Hamiltonian_SlaterKosterFiles_Prefix=slako_dir,
                Hamiltonian_SlaterKosterFiles_Separator='"-"',
                Hamiltonian_SlaterKosterFiles_Suffix='".skf"',
                Hamiltonian_MaxAngularMomentum_='')
        self.pcpot = None
        self.lines = None
        self.atoms = None
        self.atoms_input = None
        self.outfilename = 'dftb.out'
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms,
                                  **kwargs)
        self.kpts = kpts
        # kpoint stuff by ase
        if self.kpts is not None:
            mpgrid = kpts2mp(atoms, self.kpts)
            mp = monkhorst_pack(mpgrid)
            initkey = 'Hamiltonian_KPointsAndWeights'
            self.parameters[initkey + '_'] = ''
            for i, imp in enumerate(mp):
                key = initkey + '_empty' + str(i)
                self.parameters[key] = str(mp[i]).strip('[]') + ' 1.0'
    def read_results(self):
        """ all results are read from results.tag file
            It will be destroyed after it is read to avoid
            reading it once again after some runtime error """
        # from ase.units import Hartree, Bohr
        # myfile = open('results.tag'), 'r')
        # self.lines = myfile.readlines()
        # myfile.close()
        
        # energy = 0.0
        # energy = self.read_energy()
        
        # self.results['energy'] = energy

        """these results are pulled from the "detailed.out" file. 
        So far I've pulled eigenvalues and fillings, which can be used for homo/lumo 
        Note that these are spin unpolarized calculations"""

        myfile2 = open("detailed.out", "r")
        self.lines2 = myfile2.readlines() #changed to lines2
        myfile2.close()
     
        eigenvalues = 0
        fillings = 0
        eigenvalues = self.read_eigenvalues()
        fillings = self.read_fillings()
        self.results['eigenvalues'] = eigenvalues
        self.results['fillings'] = fillings

        print(self.results)
        return self.results #why won't it return RIP


    def read_energy(self):
        """Read Energy from dftb output file (results.tag)."""
        from ase.units import Hartree
        # Energy line index
        for iline, line in enumerate(self.lines):
            estring = 'total_energy'
            if line.find(estring) >= 0:
                index_energy = iline + 1
                break
        try:
            return float(self.lines[index_energy].split()[0]) * Hartree
        except:
            raise RuntimeError('Problem in reading energy')

    def read_eigenvalues(self):
        
        # eigenvalues line index
        for iline, line in enumerate(self.lines2):
            evstring = 'Eigenvalues /eV'
            if line.find(evstring) >= 0:
                index_eigenvalues_start = iline + 1
                break
        try:
            eigenvals = []
            i = index_eigenvalues_start
            #print(self.lines2)
            try:
                while self.lines2[i].split()[0] != "\n":
                    eigenvals.append(self.lines2[i].split()[0])
                    i += 1
            except: 
                return eigenvals

        except:
            raise RuntimeError('Problem in reading eigenvalues')

    def read_fillings(self):
        # eigenvalues line index
        for iline, line in enumerate(self.lines2):
            fstring = 'Fillings'
            if line.find(fstring) >= 0:
                index_fillings_start = iline + 1
                break
        try:
            fillings = []
            i = index_fillings_start
            #print(self.lines2)
            try:
                while self.lines2[i].split()[0] != "\n":
                    fillings.append(self.lines2[i].split()[0])
                    i += 1
            except: 
                return fillings
        except:
            raise RuntimeError('Problem in reading fillings')


fileToTest = Dftb()
fileToTest.read_results()
    
  