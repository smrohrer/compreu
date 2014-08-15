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
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm

np.set_printoptions(precision=3,suppress=True)

def build_sheet(nx, nz):
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


def calc_edge_nitrogens(nx="1", nz="1"):
    atoms = build_sheet(nx, nz)
    daves_super_saturate(atoms)
    data = make_orca(atoms, filename="%dx%dgraphene.inp" % (nx, nz), multiplicity="1")
    os.popen("mkdir /home/matthew/compreu/%dx%dsheet" % (nx, nz))
    os.chdir("/home/matthew/compreu/%dx%dsheet" % (nx, nz))
    moenergies_array = data.moenergies[0]
    almost_LUMO_index = data.homos

##Writes carbon only sheet energy values
    with open("%dx%d_edge_results.txt" % (nx, nz), 'w') as r:
        r.write("%dx%d Graphene Sheet -- No Nitrogens\n" % (nx, nz))
        r.write("Total SCF energy in eV:\t")
        r.write(str(data.scfenergies))
        r.write("\nMolecular orbital energy of HOMO in eV:\t")
        r.write(str(moenergies_array[data.homos]))
        r.write("\nMolecular orbital energy of LUMO in eV:\t")
        r.write(str(moenergies_array[almost_LUMO_index+1]))

##Pattern for describing zig-zag edge atom positions.  Creates edge_carbon_index array based on nx & nz parameters
    addition = [0]
    multiplication = []
    edge_carbon_index =[]

    for number in xrange(0, 2*nx):
        addition.append(number)
    for number in xrange(2, 2*(nx+2), 2):
        multiplication.append(number)
        multiplication.append(number)
    print multiplication
    print addition
    for value in xrange(0, len(addition)):
        edge_carbon_index.append(nz*multiplication[value]+addition[value])
    edge_carbon_index.pop(1)

##setup arrays for color map that will be utilized in draw_colormap function later    
    x_pos, y_pos, scf_energy, HOMO_energy, LUMO_energy = (np.array([]) for i in range(5))
    x_pos = np.append(x_pos, atoms.get_positions()[:,0])
    y_pos = np.append(y_pos, atoms.get_positions()[:,2])
    scf_energy = np.append(scf_energy, data.scfenergies)
    scf_energy = np.repeat(scf_energy, atoms.get_number_of_atoms())
    HOMO_energy = np.append(HOMO_energy, moenergies_array[data.homos])
    HOMO_energy = np.repeat(HOMO_energy, atoms.get_number_of_atoms())
    LUMO_energy = np.append(LUMO_energy, moenergies_array[almost_LUMO_index+1])
    LUMO_energy = np.repeat(LUMO_energy, atoms.get_number_of_atoms())

##Writes energy values associated with single nitrogen substitions into energy arrays
    for index_number in edge_carbon_index:
        atoms = build_sheet(nx, nz)
        nitrogenate(atoms, index_number)
        daves_super_saturate(atoms)
        #view(atoms, viewer="avogadro")
        data = make_orca(atoms, filename = "%dx%dsheetN%d" % (nx, nz, index_number), multiplicity="1")
        scf_energy[index_number] = data.scfenergies
        HOMO_energy[index_number] = moenergies_array[data.homos]
        LUMO_energy[index_number] = moenergies_array[almost_LUMO_index+1]

##Writes results text file
        with open("%dx%d_edge_results.txt" % (nx, nz), 'a+') as e:
            e.write("\n\n%dx%dsheetN%d\n" % (nx, nz, index_number))
            e.write("Total SCF energy in eV:\t")
            e.write(str(data.scfenergies))
            e.write("\nMolecular orbital energy of HOMO in eV:\t")
            e.write(str(moenergies_array[data.homos]))
            e.write("\nMolecular orbital energy of LUMO in eV:\t")
            e.write(str(moenergies_array[almost_LUMO_index+1]))

    ##Creates colormaps
    cm = plt.get_cmap("hot")
    title_list =["SCF Energy Map", "HOMO Energy Map", "LUMO Energy Map"]
    energy_list = [scf_energy, HOMO_energy, LUMO_energy]
    plt.xlabel("Atom X Position on Sheet")
    plt.ylabel("Atom Y Position on Sheet")

    for i in xrange(3):
        plt.title(title_list[i])
        COLOR = (energy_list[i]*energy_list[i].min())/energy_list[i].max()
        plt.scatter(x_pos, y_pos, c=COLOR, s=100, marker='o', edgecolors='none')
        plt.savefig(title_list[i]+".png")





atoms = build_sheet(3, 3)
#edge_carbon_list(atoms, 3)


#view(atoms, viewer="avogadro")
calc_edge_nitrogens(3, 3)
#print data.atomcharges