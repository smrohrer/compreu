# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from ase import Atoms, Atom
from ase.visualize import write
from ase.io import write
import numpy as np
from scipy.spatial import KDTree
import math
import os
import random

EV_TO_KCAL = 23.0605
H2_TOTAL_ENERGY = -634.997248712
H_ENERGY = -13.6

coBondLength = 1.36
chBondLength = 1.09
ohBondLength = 0.96
ohBondAngle = 109.5 * np.pi / 180

np.set_printoptions(precision=3,suppress=True)

def base_ring():
    return Atoms('C6',
                  positions=[[3.11488, 2.50000, 0.71000],
                             [4.34463, 2.50000, 1.42000],
                             [4.34463, 2.50000, 2.84000],
                             [3.11488, 2.50000, 3.55000],
                             [1.88513, 2.50000, 1.42000],
                             [1.88513, 2.50000, 2.84000]])

def daves_super_saturate(atoms):
    pos = atoms.get_positions()
    bondedTo = build_bonded_to(pos)

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
            rH = chBondLength * rH / np.linalg.norm(rH)
            atoms.append(Atom('H',  r0+rH ))

def build_bonded_to(pos):
    tree = KDTree(pos)
    list_tree = list(tree.query_pairs(1.430))
    bondedTo = [[] for i in range(len(pos))]
    for bond in list_tree:
        bondedTo[bond[0]].append(bond[1])
        bondedTo[bond[1]].append(bond[0])
    return bondedTo

def get_edge_atoms(atoms,bondedTo):
    Zs = atoms.get_atomic_numbers()
    edge_atoms = []
    for iatom, neighbors in enumerate(bondedTo):
        if Zs[iatom] == 6 and len(neighbors) == 2:
            if Zs[neighbors[0]] == 6 and Zs[neighbors[1]] == 6:
                edge_atoms.append(iatom)
    return edge_atoms

def add_alcohol(atoms):
    pos = atoms.get_positions()
    bondedTo = build_bonded_to(pos)
    edge_atoms = get_edge_atoms(atoms,bondedTo)

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
    return atoms

def add_COOH(atoms):
    pos = atoms.get_positions()
    bondedTo = build_bonded_to(pos)
    edge_atoms = get_edge_atoms(atoms,bondedTo)

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
    CO_bond = coBondLength * CO_bond / np.linalg.norm(CO_bond)
    atoms.append(Atom('O', r0 + CC_bond + CO_bond))
    angle3 = 109 * 2 * np.pi / 360
    OH_rotation_axis = np.cross(CO_bond,CO_double_bond)
    rotation3 = rotation_matrix(OH_rotation_axis,angle3)
    OH_bond = np.dot(-CO_bond, rotation3)
    OH_bond = ohBondLength * OH_bond / np.linalg.norm(OH_bond)
    atoms.append(Atom('H', r0 + CC_bond + CO_bond + OH_bond))
    return atoms

def add_pyridinic(atoms):
    pos = atoms.get_positions()
    bondedTo = build_bonded_to(pos)

    edge_atoms = get_edge_atoms(atoms,bondedTo)

    # choose a random edge atom to replace with nitrogen
    atom = random.choice(edge_atoms)
    Zs = atoms.get_atomic_numbers()
    Zs[atom] = 7
    atoms.set_atomic_numbers(Zs)
    return atoms

def add_pyrrollic(atoms):
    pos = atoms.get_positions()
    bondedTo = build_bonded_to(pos)

    # find a random pair of edge atoms to be replaced with a single nitrogen
    edge_atoms = get_edge_atoms(atoms,bondedTo)
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
    return atoms

def add_epoxide(atoms,basal):
    pos = atoms.get_positions()
    bonded_to = build_bonded_to(pos)
    
    bond_length = 1.43
    epoxide_dist = 1.24

    # Find a random pair of adjacent atoms to get an epoxide
    # Only if pair is a fully exposed edge (eg armchair) place in plane

    C1 = random.choice(basal)
    C2 = random.choice(bonded_to[C1])
    middle_pos = (pos[C1] + pos[C2]) / 2
    # If the selected carbons are both edge (ie have only 2 bonds)
    if len(bonded_to[C1]) == 2 and len(bonded_to[C2]) == 2:
        # Place epoxide in plane of sheet
        carbs=[C1,C2,bonded_to[C1][0],bonded_to[C1][1],bonded_to[C2][0],bonded_to[C2][1]]
        carbs=np.unique(carbs)
        avPos = 0
        for carb in carbs:
            avPos += pos[carb]
        avPos/=len(carbs)
        bond_vect = middle_pos - avPos
        bond_vect = epoxide_dist * bond_vect / np.linalg.norm(bond_vect)
        o_pos = middle_pos + bond_vect
    else:
        # Place above or below sheet
        o_pos = middle_pos
        o_pos[1] = o_pos[1] + random.choice([-1, 1]) * epoxide_dist
    atoms.append(Atom('O', o_pos))
    return atoms

def build_random_ring_structure(rings):
    atoms = base_ring()
    for iring in range(0, rings-1):
        pos = atoms.get_positions()
        bondedTo = build_bonded_to(pos)

        edge_atoms = get_edge_atoms(atoms,bondedTo)
        
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

    return atoms

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


def random_structure(rings=1, pyrroles=0, nitrogens=0, alcohols=0, COOH=0, epoxide=0):
    atoms = build_random_ring_structure(rings)
    for ipyrrole in range(pyrroles):
        atoms = add_pyrrollic(atoms)
    for iN in range(nitrogens): # pyridinic N
        atoms = add_pyridinic(atoms)
    basal_Cs = range(len(atoms))
    to_pop = []
    # Find all basal carbons. Remove from index list from last to first
    for i,j in enumerate(atoms.get_chemical_symbols()):
        if j == 'N':
            to_pop.append(i)
    to_pop.reverse()
    for i in to_pop:
        basal_Cs.pop(i)
    for iOH in range(alcohols):
        atoms = add_alcohol(atoms)
    for iCOOH in range(COOH):
        atoms = add_COOH(atoms) 
    for iEpoxi in range(epoxide):
        atoms = add_epoxide(atoms,basal_Cs)

    daves_super_saturate(atoms)
    atoms1 = atoms.copy()
    return atoms1

def name_structure(atoms, rings, pyrroles, nitrogens, alcohols, COOH, epoxide):
    chem_str = atoms.get_chemical_formula()
    struct_input = '_r{0}'.format(rings)
    if pyrroles != 0:
        struct_input += '_5py{0}'.format(pyrroles)
    if nitrogens != 0:
        struct_input += '_6py{0}'.format(nitrogens)
    if alcohols != 0:
        struct_input += '_OH{0}'.format(alcohols) 
    if COOH != 0:
        struct_input += '_COOH{0}'.format(COOH)
    if epoxide != 0:
        struct_input += '_ep{0}'.format(epoxide)
    final_str = chem_str + struct_input
    return final_str

