# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from ase import Atoms
from ase.visualize import view
import numpy as np
from cclib.parser import ccopen
from cllib.parser import ORCA
import logging


np.set_printoptions(precision=3,suppress=True)

def build_sheet(nx, nz='nx'):
    nx=nx+1
    basic_cell= Atoms('C4', 
        positions=[[3.11488, 2.50000, 0.71000], 
               [4.34463, 2.50000, 1.42000], 
               [4.34463, 2.50000, 2.84000], 
               [3.11488, 2.50000, 3.55000]],
          cell=[[2.45951, 0, 0], 
                [0, 1, 0],
                [0, 0, 4.26]])
    atoms=basic_cell.repeat((nx,1,nz))
    atoms.pop(basic_cell.get_number_of_atoms()*nz-1)
    atoms.pop(0)
    return atoms

def nitrogenate(sheet, position):
    symbols = atoms.get_chemical_symbols()
    symbols[position] = 'N'
    sheet.set_chemical_symbols(symbols)
    
def print_atoms(atoms):
    out=''
    for atom in atoms:
        atom_str= '{0}\t{1}\t{2}\t{3}\n'.format(atom.symbol, atom.position[0], atom.position[1], atom.position[2])
        out=out+atom_str
    return out

def orca_parameters(charge="0", multiplicity="2"):
    # if sum(atoms.get_atomic_numbers()) % 2 == 1:
    #        spin= "! HF"
    # elif sum(atoms.get_atomic_numbers()) % 2 == 0:
    #        spin= "! UHF"
    out=''
    parameters0= '{0}\t{1}\t{2}\n{3}\n'.format("%method", "method", "am1", "end")
    parameters1='{0}\t{1}\t{2}\t{3}'.format("*", "xyz", charge, multiplicity)
    out=out+parameters0+parameters1
    return out
    
def parse(filename, cclib_attribute):
    filename= "filename.txt"
    myfile= ccopen(filename)
    data = myfile.parse()
    data.cclib_attribute

    
atoms=build_sheet(1,1)
nitrogenate(atoms, 0)
print orca_parameters()
print print_atoms(atoms)
#view(atoms, viewer='avogadro')

# <codecell>

