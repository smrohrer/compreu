# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from ase import Atoms, Atom
from ase.visualize import write
import numpy as np
from cclib.parser import *
import logging
from scipy.spatial import KDTree
from collections import Counter
import subprocess
import math
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import re
import random
# from Tkinter import *

EV_TO_KCAL = 23.0605
H2_TOTAL_ENERGY = -634.997248712
H_ENERGY = -13.6

np.set_printoptions(precision=3,suppress=True)



    

def daves_super_saturate(atoms):
    pos = atoms.get_positions()
    tree = KDTree(atoms.get_positions())
    list_tree = list(tree.query_pairs(1.430))
    bondedTo = [ [] for i in range(len(atoms))] 

    for bond in list_tree:
        bondedTo[bond[0]].append(bond[1])
        bondedTo[bond[1]].append(bond[0])

    Zs = atoms. get_atomic_numbers()
# figure out what needs a hydrogen atom
    for iatom in range(len(atoms)):
        nbonds = len( bondedTo[iatom] )
        Z = Zs[iatom]
        if (Z,nbonds) == (6,2):
            print("we should add H to atom ", iatom)
            

            r0 = pos[iatom, :]
            bond1 = pos[ bondedTo[iatom][0] , : ] - r0
            bond2 = pos[ bondedTo[iatom][1],   :]  -r0
            rH = -(bond1 + bond2)
            rH = 1.09 * rH / np.linalg.norm(rH)
            atoms.append(Atom('H',  r0+rH ))


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def random_structure(rings=1, pyrroles=0, nitrogens=1, alcohols=0, nCOOH=1):

    # start with a single benzene ring
    atoms = Atoms('C6',
                  positions=[[3.11488, 2.50000, 0.71000],
                             [4.34463, 2.50000, 1.42000],
                             [4.34463, 2.50000, 2.84000],
                             [3.11488, 2.50000, 3.55000],
                             [1.88513, 2.50000, 1.42000],
                             [1.88513, 2.50000, 2.84000]])

    for iring in range(0, rings-1):
        pos = atoms.get_positions()
        tree = KDTree(pos)
        list_tree = list(tree.query_pairs(1.430))
        bondedTo = [[] for i in range(len(atoms))]
        for bond in list_tree:
            bondedTo[bond[0]].append(bond[1])
            bondedTo[bond[1]].append(bond[0])

        # find a random pair of edge atoms where a new ring will be fused
        edge_atoms = []
        for iatom, neighbors in enumerate(bondedTo):
            if len(neighbors) == 2:
                edge_atoms.append(iatom)
        C1 = -1
        C2 = -1
        anchor1 = []
        while C1 < 0 or C2 < 0:
            C1 = random.choice(edge_atoms)
            if len(bondedTo[bondedTo[C1][0]]) == 2:
                C2 = bondedTo[C1][0]
                anchor1 = pos[C1] - pos[bondedTo[C1][1]]
            elif len(bondedTo[bondedTo[C1][1]]) == 2:
                C2 = bondedTo[C1][1]
                anchor1 = pos[C1] - pos[bondedTo[C1][0]]
        bond = pos[C2] - pos[C1]
        anchor2 = - np.dot(2 * np.dot(anchor1, bond) / np.dot(bond, bond), bond) + anchor1
        atoms.append(Atom('C', pos[C1] + anchor2))
        C3 = len(atoms) - 1
        atoms.append(Atom('C', pos[C2] + anchor1))
        C4 = len(atoms) - 1
        pos = atoms.get_positions()
        atoms.append(Atom('C', pos[C3] + anchor1))
        atoms.append(Atom('C', pos[C4] + anchor2))

    for ipyrrole in range(pyrroles):
        pos = atoms.get_positions()
        tree = KDTree(pos)
        list_tree = list(tree.query_pairs(1.430))
        bondedTo = [[] for i in range(len(atoms))]
        for bond in list_tree:
            bondedTo[bond[0]].append(bond[1])
            bondedTo[bond[1]].append(bond[0])

        # find a random pair of edge atoms to be replaced with a single nitrogen
        Zs = atoms.get_atomic_numbers()
        edge_atoms = []
        for iatom, neighbors in enumerate(bondedTo):
            if Zs[iatom] == 6 and len(neighbors) == 2:
                if Zs[neighbors[0]] == 6 and Zs[neighbors[1]] == 6:
                    edge_atoms.append(iatom)
        C1 = -1
        C2 = -1
        while C1 < 0 or C2 < 0:
            C1 = random.choice(edge_atoms)
            if len(bondedTo[bondedTo[C1][0]]) == 2:
                C2 = bondedTo[C1][0]
                C3 = bondedTo[C1][1]
            elif len(bondedTo[bondedTo[C1][1]]) == 2:
                C2 = bondedTo[C1][1]
                C3 = bondedTo[C1][0]
        if bondedTo[C2][0]==C1:
            C4 = bondedTo[C2][1]
        else:
            C4 = bondedTo[C2][0]
        N_pos = (pos[C1] + pos[C2] + pos[C3] + pos[C4]) / 4
        map = [False for i in range(len(atoms))]
        map[C1] = True
        map[C2] = True
        atoms.pop(map.index(True))
        map.remove(True)
        atoms.pop(map.index(True))
        atoms.append(Atom('N', N_pos))

    for iN in range(nitrogens):
        pos = atoms.get_positions()
        tree = KDTree(pos)
        list_tree = list(tree.query_pairs(1.430))
        bondedTo = [[] for i in range(len(atoms))]
        for bond in list_tree:
            bondedTo[bond[0]].append(bond[1])
            bondedTo[bond[1]].append(bond[0])

        # find a edge atoms
        Zs = atoms.get_atomic_numbers()
        edge_atoms = []
        for iatom, neighbors in enumerate(bondedTo):
            if Zs[iatom] == 6 and len(neighbors) == 2:
                if Zs[neighbors[0]] == 6 and Zs[neighbors[1]] == 6:
                    edge_atoms.append(iatom)

        # choose a random edge atom to replace with nitrogen
        atom = random.choice(edge_atoms)
        Zs[atom] = 7
        atoms.set_atomic_numbers(Zs)

    coBondLength = 1.36
    ohBondLength = 0.96
    ohBondAngle = 109.5 * np.pi / 180
    for iO in range(alcohols):
        pos = atoms.get_positions()
        tree = KDTree(pos)
        list_tree = list(tree.query_pairs(1.430))
        bondedTo = [[] for i in range(len(atoms))]
        for bond in list_tree:
            bondedTo[bond[0]].append(bond[1])
            bondedTo[bond[1]].append(bond[0])

        # find a edge atoms
        Zs = atoms.get_atomic_numbers()
        edge_atoms = []
        for iatom, neighbors in enumerate(bondedTo):
            if Zs[iatom] == 6 and len(neighbors) == 2:
                if Zs[neighbors[0]] == 6 and Zs[neighbors[1]] == 6:
                    edge_atoms.append(iatom)

        # choose a random edge atom to attach an alcohol
        atom = random.choice(edge_atoms)

        r0 = pos[atom]
        bond1 = pos[bondedTo[atom][0]] - r0
        bond2 = pos[bondedTo[atom][1]] - r0
        CO_bond = -(bond1 + bond2)
        CO_bond = coBondLength * CO_bond / np.linalg.norm(CO_bond)
        atoms.append(Atom('O', r0 + CO_bond))
        rotation = [[np.cos(ohBondAngle), 0, np.sin(ohBondAngle)],
                    [0, 1, 0],
                    [-np.sin(ohBondAngle), 0, np.cos(ohBondAngle)]]
        OH_bond = np.dot(-CO_bond, rotation)
        OH_bond = ohBondLength * OH_bond / np.linalg.norm(OH_bond)
        atoms.append(Atom('H', r0 + CO_bond + OH_bond))

    for iO in range(nCOOH):
        atom = random.choice(edge_atoms)

        r0 = pos[atom]
        bond1 = pos[bondedTo[atom][0]] - r0
        bond2 = pos[bondedTo[atom][1]] - r0
        CC_bond = -(bond1 + bond2)
        CC_bond = 1.4 * CC_bond / np.linalg.norm(CC_bond)
        atoms.append(Atom('C', r0 + CC_bond))
        
        dihedralRotationUp = rotation_matrix(CC_bond, np.pi/2)
        dihedralRotationDown = rotation_matrix(CC_bond, -1 * np.pi/2)

        angle1 = 120 * 2.0 * np.pi / 360
        rotation1 = [[np.cos(angle1), 0, np.sin(angle1)],
                     [0, 1, 0],
                     [-np.sin(angle1), 0, np.cos(angle1)]]
        CO_double_bond = np.dot(-CC_bond, rotation1)
        CO_double_bond = np.dot(CO_double_bond, dihedralRotationUp)
        CO_double_bond = 1.2 * CO_double_bond / np.linalg.norm(CO_double_bond)
        atoms.append(Atom('O', r0 + CC_bond + CO_double_bond))
        angle2 = -120 * 2.0 * np.pi / 360
        rotation2 = [[np.cos(angle2), 0, np.sin(angle2)],
                     [0, 1, 0],
                     [-np.sin(angle2), 0, np.cos(angle2)]]
        CO_bond = np.dot(-CC_bond, rotation2)
        CO_bond = np.dot(CO_bond, dihedralRotationUp)
        CO_bond = 1.4 * CO_bond / np.linalg.norm(CO_bond)
        atoms.append(Atom('O', r0 + CC_bond + CO_bond))
        angle3 = 109 * 2 * np.pi / 360
        OH_rotation_axis = np.cross(CO_bond,CO_double_bond)
        rotation3 = rotation_matrix(OH_rotation_axis,angle3)
        OH_bond = np.dot(-CO_bond, rotation3)
        OH_bond = 0.96 * OH_bond / np.linalg.norm(OH_bond)
        atoms.append(Atom('H', r0 + CC_bond + CO_bond + OH_bond))
        

    daves_super_saturate(atoms)

    atoms1 = atoms.copy()
    atoms1.rotate("x", np.pi / 2.0)
    write("random_structure.png", atoms1)

    xyz = open("test.xyz", "w")
    xyz.write(str(len(atoms1))+"\n")
    xyz.write("comment\n")
    symbols_temp = atoms1.get_chemical_symbols()
    pos_temp = atoms1.get_positions()
    for iatom in range(0,len(atoms1)):
        xyz.write(symbols_temp[iatom]+"\t"+str(pos_temp[iatom][0])+"\t"+str(pos_temp[iatom][1])+"\t"+str(pos_temp[iatom][2])+"\n")
    xyz.close()

    return atoms1


random_structure(rings=20, pyrroles=0, nitrogens=5, alcohols=0, nCOOH=1)

