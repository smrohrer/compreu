# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from ase import Atoms, Atom
from ase.visualize import view
import numpy as np
from cclib.parser import *
import logging
from scipy.spatial import KDTree
from collections import Counter
import subprocess


np.set_printoptions(precision=3,suppress=True)

def build_sheet(nx, nz='nx'):
    nx=nx+1
    nz=(nz+1)/2
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
    return position
    
def print_atoms(atoms):
    out=''
    for atom in atoms:
        atom_str= '{0}\t{1}\t{2}\t{3}\n'.format(atom.symbol, atom.position[0], atom.position[1], atom.position[2])
        out=out+atom_str
    return out

def make_orca(atoms, filename="filename.inp", charge="0", multiplicity="1", method="am1"):
    # if sum(atoms.get_atomic_numbers()) % 2 == 1:
    #        spin= "! HF"
    # elif sum(atoms.get_atomic_numbers()) % 2 == 0:
    #        spin= "! UHF"
    out=''
    parameters0= '{0}\t{1}\t{2}\n{3}\n'.format("%method", "method", method, "end")
    parameters1='{0}\t{1}\t{2}\t{3}\n'.format("*", "xyz", charge, multiplicity)
    out=out+parameters0+parameters1
    end_of_atom_coordinates="*"
    with open(filename, 'w') as f:
        f.write(out)
        f.write(print_atoms(atoms))
        f.write(end_of_atom_coordinates)
    subprocess.call("/home/matthew/orca/orca "+ filename + " > temp.out", shell=True)
    return parse("temp.out")

def parse(filename):
    myfile= ccopen(filename)
    data = myfile.parse()
    return data

def saturate(atoms = 'atoms'):
    __bond3__, __bond2__, __saturate__ = (np.array([]) for i in range(3))
    
    tree = KDTree(atoms.get_positions())
    list_tree = list(tree.query_pairs(1.430))
    print list_tree
    dictionary_count=Counter(elem[0] for elem in list_tree) + Counter(elem[1] for elem in list_tree)
    for k, v in dictionary_count.iteritems():
        if v == 3:
            __bond3__ = np.append(__bond3__, k)
        elif v == 2:
            __bond2__ = np.append(__bond2__, k)
        else:
            pass
    print __bond3__
    print __bond2__
    for i in range(__bond2__.size):
        __saturate__ = np.append(__saturate__, (atoms.get_positions()[__bond2__.item(i),:]))
    __saturate__ = np.reshape(__saturate__, (__saturate__.size/3, 3))
    print __saturate__
    

def daves_super_saturate(atoms):
    pos = atoms.get_positions()
    tree = KDTree(atoms.get_positions())
    list_tree = list(tree.query_pairs(1.430))
    bondedTo = [ [] for i in xrange(len(atoms))] 

    for bond in list_tree:
        bondedTo[bond[0]].append(bond[1])
        bondedTo[bond[1]].append(bond[0])

    Zs = atoms. get_atomic_numbers()
# figure out what needs a hydrogen atom
    for iatom in xrange(len(atoms)):
        nbonds = len( bondedTo[iatom] )
        Z = Zs[iatom]
        if (Z,nbonds) == (6,2):
            print "we should add H to atom ", iatom
            

            r0 = pos[iatom, :]
            bond1 = pos[ bondedTo[iatom][0] , : ] - r0
            bond2 = pos[ bondedTo[iatom][1],   :]  -r0
            rH = -(bond1 + bond2)
            rH = 1.09 * rH / np.linalg.norm(rH)
            atoms.append(Atom('H',  r0+rH ))
            
def find_edge_atoms(atoms):
    edge_atoms = []
    pos = atoms.get_positions()
    tree = KDTree(atoms.get_positions())
    list_tree = list(tree.query_pairs(1.430))
    bondedTo = [ [] for i in xrange(len(atoms))] 

    for bond in list_tree:
        bondedTo[bond[0]].append(bond[1])
        bondedTo[bond[1]].append(bond[0])

    Zs = atoms. get_atomic_numbers()
    # figure out what needs a hydrogen atom
    for iatom in xrange(len(atoms)):
        nbonds = len( bondedTo[iatom] )
        Z = Zs[iatom]
        if (Z,nbonds) == (6,2):
            edge_atoms.append(iatom)
    return edge_atoms


    
atoms = build_sheet(3,3)
daves_super_saturate(atoms)
data = make_orca(atoms, filename="3x3graphene.inp")

with open("results.txt", 'w') as r:
    r.write("3x3 Graphene Sheet -- No Nitrogens\n")
    r.write("Total SCF energy in eV:\t")
    r.write(str(data.scfenergies))
    r.write("\nMolecular orbital energy of HOMO in eV:\t")
    moenergies_array = data.moenergies[0]
    r.write(str(moenergies_array[data.homos]))


#view(atoms, viewer="avogadro")
#print data.atomcharges

edge_carbon_index = [6, 14, 22, 13, 21, 29]
for index_number in edge_carbon_index:
    atoms = build_sheet(3,3)
    nitrogenate(atoms, index_number)
    daves_super_saturate(atoms)
    data = make_orca(atoms, filename = "3x3sheetN%d" % index_number)
    with open("results.txt", 'a+') as r:
        r.write("\n\n3x3sheetN%d\n" % index_number)
        r.write("Total SCF energy in eV:\t")
        r.write(str(data.scfenergies))
        r.write("\nMolecular orbital energy of HOMO in eV:\t")
        moenergies_array = data.moenergies[0]
        r.write(str(moenergies_array[data.homos]))