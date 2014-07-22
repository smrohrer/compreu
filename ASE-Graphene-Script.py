# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from ase import Atoms
from ase.visualize import view
import numpy as np

np.set_printoptions(precision=3,suppress=True)

def build_sheet(nx, nz):
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

def calc_parameters(calculation_type, basis_set, kwargs, charge, multiplicity):
    print string.rjust("!", 4), string.rjust(calculation_type, 4), string.rjust(basis_set, 4), string.rjust(kwargs, 4)
    print string.rjust("* xyz", 4), string.rjust(charge, 4), string.rjust(multiplicity, 4)
    

    
atoms=build_sheet(1,1)
nitrogenate(atoms, 0)
print print_atoms(atoms)
view(atoms, viewer='avogadro')

# <codecell>

