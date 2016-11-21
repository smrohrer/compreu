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
ORCA_filepath = os.getcwd()

def build_sheet(nx, nz, symmetry=0):
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

    if symmetry == 1:
        print("Making graphene sheet symmetric")
        global to_be_removed
        to_be_removed = []
        sym_fix_int = nz-1
        for i in range(2, nz+sym_fix_int, 2):
            to_be_removed.append(2*i-2)
            to_be_removed.append(2*i-1)
        
        to_be_removed = to_be_removed[::-1]

        for entry in to_be_removed:
            atoms.pop(entry)
    else:
        print("Molecule was not made symmetric")
    return atoms


def nitrogenate(sheet, position):
    symbols = sheet.get_chemical_symbols()
    symbols[position] = 'N'
    sheet.set_chemical_symbols(symbols)
    return position
    
def print_atoms(atoms):
    out=''
    for atom in atoms:
        atom_str= '{0}\t{1}\t{2}\t{3}\n'.format(atom.symbol, atom.position[0], atom.position[1], atom.position[2])
        out=out+atom_str
    return out

def make_orca(atoms, filename="filename.inp", charge="0", multiplicity="1", method="am1", geometry_opt=0,
              output="temp.out", compute_cores=3, convergence_tolerance="default"):
    parameters = '%method\n\tMethod AM1\n'
    if geometry_opt==1:
        parameters = parameters + '\tRunTyp Opt\n'
    parameters = parameters + 'end\n' \
                              '%scf\n' \
                              '\tMaxIter 2000\n' \
                              'end\n' \
                              '* xyz '+charge+' '+multiplicity+'\n'

    # if method == "am1":
    #     if geometry_opt == 0:
    #         parameters1= '{0}\t{1}\t{2}\n\t{3}'.format("%method", "method", method, "end")
    #     elif geometry_opt == 1:
    #         if convergence_tolerance == "default":
    #             parameters1= '{0}\n\t{1}\n\t{2}\n\n{3}\n{4}\n\n{5}\n{6}\n\n'.format("%method method am1", "method OPT", "end", "%SCF MAXITER 300", "     end", "%geom MaxIter 250", "      end")
    #         else:
    #             parameters1= '{0}\n\t{1}\n\t{2}\n\n{3}\n{4}\n{5}\n\n{6}\n{7}\n\n'.format("%method method am1", "method OPT", "end", "%SCF MAXITER 300", "     Convergence " + convergence_tolerance, "     end", "%geom MaxIter 250", "      end")
    #
    # elif method == "DFT":
    #     if geometry_opt == 0:
    #         parameters1 = '{0}\t{1}\n'.format("!", "DFT-Energy")
    #     elif geometry_opt == 1:
    #         parameters1= '{0}\t{1}\n'.format("!", "Good-Opt")
    #
    # out = ''
    # parameters0 = '{0}\n{1}\n\n'.format("%pal nprocs " + str(compute_cores), "     end")
    # parameters0 = ''
    # parameters2 ='{0}\t{1}\t{2}\t{3}\n'.format("*", "xyz", charge, multiplicity)
    # out=out+parameters0+parameters1+parameters2
    # end_of_atom_coordinates="*\n" \
    #                         "%Output\n" \
    #                         "\tPrint[P_MOs]\t1\n" \
    #                         "\tPrint[P_Overlap]\t1\n" \
    #                         "end\n"
    end_of_atom_coordinates = "*"
    with open(filename, 'w') as f:
        f.write(parameters)
        f.write(print_atoms(atoms))
        f.write(end_of_atom_coordinates)
    subprocess.call("$HOME/Library/Orca/orca "+ filename + " > " + output, shell=True)
    return parse(output)

def parse(filename):
    myfile = ccopen(filename)
    data = myfile.parse()
    return data

def make_gaussian(atoms, method="am1", geometry_opt=False, output="temp.chk"):
    if method == "am1" and geometry_opt == False:
        out = ''
        parameters0 = '{0}\n{1}\n{2}'.format("%mem=48MB", "%chk=/scratch/", output)

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
            
def find_edge_atoms(atoms):
    edge_atoms = []
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
            edge_atoms.append(iatom)
    return edge_atoms


def calc_edge_nitrogens(nx="1", nz="1", method="am1", optimize_geometry=0, make_symmetric=0, sub_all_zigzag=0):
    if make_symmetric == 1:
        print("atoms sheet at start of calc is symmetric")
        atoms = build_sheet(nx, nz, symmetry=1)
        symmetry_folder_string = "_symmetric"
        ##set parameter for multiplicity of 2 here##########
    else:
        atoms = build_sheet(nx, nz, symmetry=0)
        symmetry_folder_string = "_asymmetric"

    no_hydrogen = atoms.get_positions()
    no_hydrogen_count = atoms.get_number_of_atoms()
    daves_super_saturate(atoms)
    os.popen("mkdir " + ORCA_filepath + "/%dx%dsheet%s" % (nx, nz, symmetry_folder_string))
    os.chdir(ORCA_filepath + "/%dx%dsheet%s" % (nx, nz, symmetry_folder_string))
    data = make_orca(atoms, filename="%dx%dgraphene.inp" % (nx, nz), multiplicity="1", method=method, geometry_opt=optimize_geometry, output= ORCA_filepath + "/%dx%dsheet%s/orca_%dx%dsheet.out" % (nx, nz, symmetry_folder_string, nx, nz))
    moenergies_array = data.moenergies[0]

##Writes carbon only sheet energy values
    with open("%dx%d_edge_results.txt" % (nx, nz), 'w') as r:
        r.write("%dx%d Graphene Sheet -- No Nitrogens\n" % (nx, nz))
        r.write("Total SCF energy in eV:\t")
        r.write(str(data.scfenergies))
        r.write("\nMolecular orbital energy of HOMO in eV:\t")
        r.write(str(moenergies_array[data.homos]))
        r.write("\nMolecular orbital energy of LUMO in eV:\t")
        r.write(str(moenergies_array[data.homos+1]))


##Pattern for describing zig-zag edge atom positions.  Creates edge_carbon_index array based on nx & nz parameters
    addition = [0]
    multiplication = []
    edge_carbon_index =[]

    for number in range(0, 2*nx):
        addition.append(number)
    for number in range(2, 2*(nx+2), 2):
        multiplication.append(number)
        multiplication.append(number)
    for value in range(0, len(addition)):
        edge_carbon_index.append(nz*multiplication[value]+addition[value])
    edge_carbon_index.pop(1)

    if make_symmetric == 1:
        edge_carbon_index[:] = [x-len(to_be_removed) for x in edge_carbon_index]
    else:
        pass

##setup arrays for color map that will be utilized in draw_colormap function later    
    x_pos, y_pos, scf_energy, HOMO_energy, LUMO_energy = (np.array([]) for i in range(5))
    x_pos = np.append(x_pos, no_hydrogen[:,0])
    y_pos = np.append(y_pos, no_hydrogen[:,2])
    all_carbon_scf = data.scfenergies

##Writes energy values associated with single nitrogen substitions into energy arrays
    nitrogenated_x_pos, nitrogenated_y_pos, nitrogenated_scf, nitrogenated_HOMO, nitrogenated_LUMO = (np.array([]) for i in range(5))
    for index_number in edge_carbon_index:
        if make_symmetric == 1:
            print("nitrogenate make atoms sheet is symmetric")
            atoms = build_sheet(nx, nz, symmetry=1)
        else:
            atoms = build_sheet(nx, nz, symmetry=0)
        nitrogenate(atoms, index_number)
        daves_super_saturate(atoms)
        data = make_orca(atoms, filename = "%dx%dsheetN%d" % (nx, nz, index_number), multiplicity="1", method=method, geometry_opt=optimize_geometry, output=ORCA_filepath + "/%dx%dsheet%s/orca_%dx%dsheet.out" % (nx, nz, symmetry_folder_string, nx, nz))
        moenergies_array = data.moenergies[0]
    ##Segregate all carbon energies from substituted nitrogen energies within all pertaining arrays
        
        ##X & Y position array segregation
        nitrogenated_x_pos = np.append(nitrogenated_x_pos, x_pos[index_number])
        nitrogenated_y_pos = np.append(nitrogenated_y_pos, y_pos[index_number])

        ##Energies arrays segregation
        nitrogenated_scf = np.append(nitrogenated_scf, (data.scfenergies-all_carbon_scf))
        nitrogenated_HOMO = np.append(nitrogenated_HOMO, moenergies_array[data.homos])
        nitrogenated_LUMO = np.append(nitrogenated_LUMO, moenergies_array[data.homos+1])

##Writes nitrogen substition energies into results text file
        with open("%dx%d_edge_results.txt" % (nx, nz), 'a+') as e:
            moenergies_array = data.moenergies[0]
            e.write("\n\n%dx%dsheetN%d\n" % (nx, nz, index_number))
            e.write("Total SCF energy in eV:\t")
            e.write(str(data.scfenergies))
            e.write("\nMolecular orbital energy of HOMO in eV:\t")
            e.write(str(moenergies_array[data.homos]))
            e.write("\nMolecular orbital energy of LUMO in eV:\t")
            e.write(str(moenergies_array[data.homos+1]))

    ##Creates colormaps
    cm = plt.get_cmap("hot")
    title_list =["SCF_Energy_Map", "HOMO_Energy_Map", "LUMO_Energy_Map"]
    energy_list = [nitrogenated_scf, nitrogenated_HOMO, nitrogenated_LUMO]
    plt.xlabel("Atom X Position on Sheet")
    plt.ylabel("Atom Z Position on Sheet")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(3):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.xlabel("Atom X Position on Sheet")
        plt.ylabel("Atom Z Position on Sheet")
        plt.title(title_list[i])
        plt.axis('equal')
        #COLOR = (energy_list[i]-energy_list[i].min())/(energy_list[i].max()-energy_list[i].min())
        COLOR = energy_list[i]
        ax.scatter(x_pos, y_pos, c="0.5", s=100, marker='o', edgecolors='none')
        p = ax.scatter(nitrogenated_x_pos, nitrogenated_y_pos, c=COLOR, s=100, marker='o', edgecolors='none', label="something")
        cbar = plt.colorbar(p)
        cbar.set_label("Energy in eV")
        plt.savefig(title_list[i]+".png")
        plt.clf()

def get_outer_zigzag_atoms(nx, nz):
    addition = [0]
    multiplication = []
    edge_carbon_index = []

    for number in range(0, 2 * nx):
        addition.append(number)
    for number in range(2, 2 * (nx + 2), 2):
        multiplication.append(number)
        multiplication.append(number)
    for value in range(0, len(addition)):
        edge_carbon_index.append(nz * multiplication[value] + addition[value])
    edge_carbon_index.pop(1)
    edge_carbon_index[:] = [x - len(to_be_removed) for x in edge_carbon_index]
    return edge_carbon_index

def nitrogenate_all_zig_zag(nx_min, nx_max, nz_min, nz_max, method="am1", optimize_geometry=0, make_symmetric=0, saturate_nitrogens=0):
    all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list = (np.array([]) for dummy_var in range(3))
    pltylabel_list = "SCF Energy", "HOMO Energy", "LUMO Energy"
    sheet_dimensions, sheet_count = ([] for dummy_var in range(2))
    for nx in range(nx_min, nx_max+1, 2):
        for nz in range(nz_min, nz_max+1, 2):

            if make_symmetric == 1:
                atoms = build_sheet(nx, nz, symmetry=1)
                symmetry_folder_string = "_symmetric"
                edge_carbon_index =get_outer_zigzag_atoms(nx,nz)

            elif make_symmetric == 0:
                atoms = build_sheet(nx, nz, symmetry=0)
                symmetry_folder_string = "_asymmetric"

            if saturate_nitrogens == 1:
                for edge_carbon in edge_carbon_index:
                    symbols = atoms.get_chemical_symbols()
                    symbols[edge_carbon] = 'N'
                    atoms.set_chemical_symbols(symbols)
                pos = atoms.get_positions()
                tree = KDTree(atoms.get_positions())
                list_tree = list(tree.query_pairs(1.430))
                bondedTo = [[] for i in range(len(atoms))]
                for bond in list_tree:
                    bondedTo[bond[0]].append(bond[1])
                    bondedTo[bond[1]].append(bond[0])
                Zs = atoms.get_atomic_numbers()
                for iatom in range(len(atoms)):
                    nbonds = len(bondedTo[iatom])
                    Z = Zs[iatom]
                    if (Z,nbonds) == (7,2):
                        r0 = pos[iatom]
                        bond1 = pos[ bondedTo[iatom][0]] - r0
                        bond2 = pos[ bondedTo[iatom][1]]  -r0
                        rH = -(bond1 + bond2)
                        rH = 1.09 * rH / np.linalg.norm(rH)
                        atoms.append(Atom('H',  r0+rH ))
                        daves_super_saturate(atoms)

            else:
                daves_super_saturate(atoms)

            os.popen("mkdir " + ORCA_filepath + "/all_N_zigzag")
            os.chdir(ORCA_filepath + "/all_N_zigzag")
            write("nitrogenated_%dx%dgraphene.png" % (nx, nz), atoms)
            if (nx>4 and nx<6) and (nz>4 and nz<6):
                data = make_orca(atoms, filename="nitrogenated_%dx%dgraphene.inp" % (nx, nz), multiplicity="1", method=method, geometry_opt=optimize_geometry, output= ORCA_filepath + "/all_N_zigzag/orca%s_nitrogenated_%dx%dsheet.out" % (symmetry_folder_string, nx, nz), convergence_tolerance="Medium")
            if nx>6:
                data = make_orca(atoms, filename="nitrogenated_%dx%dgraphene.inp" % (nx, nz), multiplicity="1", method=method, geometry_opt=optimize_geometry, output= ORCA_filepath + "/all_N_zigzag/orca%s_nitrogenated_%dx%dsheet.out" % (symmetry_folder_string, nx, nz), convergence_tolerance="Loose")
            else:
                data = make_orca(atoms, filename="nitrogenated_%dx%dgraphene.inp" % (nx, nz), multiplicity="1", method=method, geometry_opt=optimize_geometry, output= ORCA_filepath + "/all_N_zigzag/orca%s_nitrogenated_%dx%dsheet.out" % (symmetry_folder_string, nx, nz))
            moenergies_array = data.moenergies[0]
            all_N_SCF_energy_list = np.append(all_N_SCF_energy_list, int(data.scfenergies))
            all_N_HOMO_energy_list = np.append(all_N_HOMO_energy_list, moenergies_array[data.homos])
            all_N_LUMO_energy_list = np.append(all_N_LUMO_energy_list, moenergies_array[data.homos+1])
            sheet_dimensions.append("%sx%s" % (nx, nz))

            with open("nitrogenated_edge_results.txt", 'a+') as e:
                e.write("\n##########################\n")
                e.write("Nitrogenated %dx%d sheet\n" % (nx, nz))
                e.write("##########################\n")
                e.write("Total SCF energy in eV:\t")
                e.write(str(data.scfenergies))
                e.write("\nMolecular orbital energy of HOMO in eV:\t")
                e.write(str(moenergies_array[data.homos]))
                e.write("\nMolecular orbital energy of LUMO in eV:\t")
                e.write(str(moenergies_array[data.homos+1]))
    for dimension in range(0, len(sheet_dimensions)):
        sheet_count.append(dimension)
    sheet_count = [x+0.2 for x in sheet_count]
    for y in range(0, 3):
        pltylist = all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list
        rectangles = plt.bar(np.arange(all_N_SCF_energy_list.size), pltylist[y], 0.4, alpha=0.4, color='b')
        plt.title("%s vs Sheet Dimensions" % pltylabel_list[y])
        plt.ylabel(pltylabel_list[y])
        plt.xticks(sheet_count, sheet_dimensions)
        plt.xlabel("Sheet Dimension")
        plt.savefig("%s.png" % pltylabel_list[y])
        plt.clf()

def prepare_nitrogenated_sheet(nx=5, nz=3):

    atoms = build_sheet(nx, nz, symmetry=1)
    edge_carbon_index = get_outer_zigzag_atoms(nx,nz)

    for edge_carbon in edge_carbon_index:
        symbols = atoms.get_chemical_symbols()
        symbols[edge_carbon] = 'N'
        atoms.set_chemical_symbols(symbols)

    return atoms

def selectively_nitrogenate(unit_cells = 1, nx=5, nz=3, method="am1", optimize_geometry=0):

    cell = 0
    all_N_SCF_energy_list = [[]]
    all_N_HOMO_energy_list = [[]]
    all_N_LUMO_energy_list = [[]]
    descriptions = [[]]

    # determine the order of positions to replace C with N
    order = []
    for i in range(0, nx/2):
        order.append(nx/2-i-1)
        order.append(nx-i-1)

    atoms = build_sheet(nx, nz, symmetry=1)
    edge_carbon_index = get_outer_zigzag_atoms(nx,nz)

    atoms_temp = atoms.copy()
    daves_super_saturate(atoms_temp)
    atoms1 = atoms_temp.copy()
    atoms1.rotate("x", np.pi/2.0)
    description = "nz"+str(nz)+"cell"+str(cell+1)+"step0"
    os.popen("mkdir " + ORCA_filepath + "/selective_nitrogenation")
    os.chdir(ORCA_filepath + "/selective_nitrogenation")
    write(description+".png", atoms1)
    data = make_orca(atoms_temp, filename=description + ".inp", multiplicity="1", method=method,
                     geometry_opt=optimize_geometry,
                     output=ORCA_filepath + "/selective_nitrogenation/" + description + ".out")
    moenergies_array = data.moenergies[0]
    all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
    all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos]))
    all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos + 1]))
    descriptions.append(description)

    # get bond structure data
    pos = atoms.get_positions()
    tree = KDTree(atoms.get_positions())
    # find pairs of atoms that are bonded; they are bonded if their distance apart is less than 1.430
    list_tree = list(tree.query_pairs(1.430))
    bondedTo = [[] for i in range(len(atoms))]
    for bond in list_tree:
        bondedTo[bond[0]].append(bond[1])
        bondedTo[bond[1]].append(bond[0])
    Cs = edge_carbon_index
    xC = [pos[i][0] for i in edge_carbon_index]
    uxC = np.unique(xC)
    # sort into pairs of atoms at the same x position
    sorted_Cs = [[] for i in range(len(uxC))]
    for C in Cs:
        sorted_Cs[list(uxC).index(pos[C][0])].append(C)
    # replace with N
    for step,position in enumerate(order):
        for k in range(0,cell+1):
            for m in range(0,2):
                iatom = sorted_Cs[position+k*cell*nx][m]
                symbols = atoms.get_chemical_symbols()
                symbols[iatom] = 'N'
                atoms.set_chemical_symbols(symbols)
                atoms_temp = atoms.copy()
            daves_super_saturate(atoms_temp)
            atoms1 = atoms_temp.copy()
            atoms1.rotate("x", np.pi/2.0)
            description = "nz"+str(nz)+"cell"+str(cell+1)+"step"+str(step+1)
            write(description+".png", atoms1)
            data = make_orca(atoms_temp, filename=description+".inp", multiplicity="1", method=method, geometry_opt=optimize_geometry, output=ORCA_filepath+"/selective_nitrogenation/"+description+".out")
            moenergies_array = data.moenergies[0]
            all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
            all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos]))
            all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos + 1]))
            descriptions.append(description)

            with open("results.txt", "a+") as e:
                e.write("\n##########################\n")
                e.write(description)
                e.write("\n##########################\n")
                e.write("Total SCF energy in eV:\t")
                e.write(str(data.scfenergies))
                e.write("\nMolecular orbital energy of HOMO in eV:\t")
                e.write(str(moenergies_array[data.homos]))
                e.write("\nMolecular orbital energy of LUMO in eV:\t")
                e.write(str(moenergies_array[data.homos + 1]))

    return descriptions, all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list


def selectively_hydrogenate_nitrogens(unit_cells=1, nx=5, nz=3, method="am1", optimize_geometry=0, oxide=0):
    all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list = ([] for dummy_var in range(3))
    descriptions, sheet_count = ([] for dummy_var in range(2))
    # determine the order of N positions to add H
    order = []
    for i in range(0, nx/2):
        order.append(nx/2-i-1)
        order.append(nx-i-1)
    for cell in range(0,unit_cells):
        all_N_SCF_energy_list.append([])
        all_N_HOMO_energy_list.append([])
        all_N_LUMO_energy_list.append([])

        atoms = prepare_nitrogenated_sheet(nx*(cell+1) if cell%2==0 else nx*(cell+1)+1, nz)
        edge_carbon_index = get_outer_zigzag_atoms(nx,nz)

        if oxide == 1:
            atoms = form_alcohol(atoms=atoms)
            atoms = form_epoxide(atoms=atoms)

        # saturate C's on armchair edges
        daves_super_saturate(atoms)
        # calculate energy of initial sheet
        atoms1 = atoms.copy()
        atoms1.rotate("x", np.pi / 2.0)
        description = "nz" + str(nz) + "cell" + str(cell + 1) + "step0"
        os.popen("mkdir " + ORCA_filepath + "/selective_hydrogenation")
        os.chdir(ORCA_filepath + "/selective_hydrogenation")
        write(description + ".png", atoms1)
        data = make_orca(atoms, filename=description + ".inp", multiplicity="1", method=method,
                         geometry_opt=optimize_geometry,
                         output=ORCA_filepath + "/selective_hydrogenation/" + description + ".out")
        moenergies_array = data.moenergies[0]
        all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies)-1])
        all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos]))
        all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos + 1]))
        descriptions.append(description)
        # get bond structure data
        pos = atoms.get_positions()
        tree = KDTree(atoms.get_positions())
        # find pairs of atoms that are bonded; they are bonded if their distance apart is less than 1.430
        list_tree = list(tree.query_pairs(1.430))
        bondedTo = [[] for i in range(len(atoms))]
        for bond in list_tree:
            bondedTo[bond[0]].append(bond[1])
            bondedTo[bond[1]].append(bond[0])
        atomic_numbers = atoms.get_atomic_numbers()
        # find N atoms and their x and z positions
        Ns = edge_carbon_index
        xN = [pos[i][0] for i in edge_carbon_index]
        uxN = np.unique(xN)
        # sort into pairs of atoms at the same x position
        sorted_Ns = [[] for i in range(len(uxN))]
        for N in Ns:
            sorted_Ns[list(uxN).index(pos[N][0])].append(N)
        # add H and do calculation
        for step,position in enumerate(order):
            changed_Ns = []
            for k in range(0,cell+1):
                for m in range(0,2):
                    iatom = sorted_Ns[position+k*cell*nx][m]
                    changed_Ns.append(iatom)
                    r0 = pos[iatom]
                    bond1 = pos[bondedTo[iatom][0]] - r0
                    bond2 = pos[bondedTo[iatom][1]] - r0
                    rH = -(bond1+bond2)
                    rH = 1.09 * rH / np.linalg.norm(rH)
                    atoms.append(Atom('H', r0+rH))
            atoms1 = atoms.copy()
            atoms1.rotate("x", np.pi / 2.0)
            description = "nz"+str(nz)+"cell"+str(cell+1)+"step"+str(step+1)
            write(description+".png", atoms1)
            data = make_orca(atoms, filename=description+".inp", multiplicity="1", method=method, geometry_opt=optimize_geometry, output=ORCA_filepath+"/selective_hydrogenation/"+description+".out")
            moenergies_array = data.moenergies[0]
            all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies)-1])
            all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos]))
            all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos + 1]))
            descriptions.append(description)

            charges = parse_mulliken(ORCA_filepath+"/selective_hydrogenation/"+description+".out")
            for i in range(0,2*(cell+1)):
                charges.pop()
            for iatom in changed_Ns:
                charges.pop(iatom)
            qsum = np.sum(charges)

            with open("results.txt", "a+") as e:
                e.write("\n##########################\n")
                e.write(description)
                e.write("\n##########################\n")
                e.write("Total SCF energy in eV:\t")
                e.write(str(data.scfenergies))
                e.write("\nMolecular orbital energy of HOMO in eV:\t")
                e.write(str(moenergies_array[data.homos]))
                e.write("\nMolecular orbital energy of LUMO in eV:\t")
                e.write(str(moenergies_array[data.homos + 1]))
                e.write("\nqsum:\t")
                e.write(str(qsum))

    cleanup(ORCA_filepath+"/selective_hydrogenation/")

    # plot_orbitals(ORCA_filepath+"/selective_hydrogenation/")
    plot_mulliken(ORCA_filepath+"/selective_hydrogenation/")
    plot_energy(all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list, unit_cells, nx, nz)

    return descriptions, all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list

def hydrogenate_inner_carbon(unit_cells=1, nx=5, nz=3, method="am1", optimize_geometry=0):
    N_SCF, N_HOMO, N_LUMO = ([] for dummy_var in range(3))
    C_SCF, C_HOMO, C_LUMO = ([] for dummy_var in range(3))
    descriptions, sheet_count = ([] for dummy_var in range(2))
    # determine the order of N positions to add H
    order = []
    for i in range(0, nx/2):
        order.append(nx/2-i-1)
        order.append(nx-i-1)
    for cell in range(0,unit_cells):
        N_SCF.append([])
        N_HOMO.append([])
        N_LUMO.append([])
        C_SCF.append([])
        C_HOMO.append([])
        C_LUMO.append([])

        atoms = prepare_nitrogenated_sheet(nx*(cell+1) if cell%2==0 else nx*(cell+1)+1, nz)
        edge_carbon_index = get_outer_zigzag_atoms(nx,nz)

        # if oxide == 1:
        #     atoms = form_alcohol(atoms=atoms)
        #     atoms = form_epoxide(atoms=atoms)

        # saturate C's on armchair edges
        daves_super_saturate(atoms)
        # calculate energy of initial sheet
        atoms1 = atoms.copy()
        atoms1.rotate("x", np.pi / 2.0)
        description = "nz" + str(nz) + "cell" + str(cell + 1) + "step0"
        os.popen("mkdir " + ORCA_filepath + "/hydrogenate_inner_carbon")
        os.chdir(ORCA_filepath + "/hydrogenate_inner_carbon")
        write(description + ".png", atoms1)
        data = make_orca(atoms, filename=description + ".inp", multiplicity="1", method=method,
                         geometry_opt=optimize_geometry,
                         output=ORCA_filepath + "/hydrogenate_inner_carbon/" + description + ".out")
        moenergies_array = data.moenergies[0]
        N_SCF[cell].append(data.scfenergies[len(data.scfenergies)-1])
        N_HOMO[cell].append(float(moenergies_array[data.homos]))
        N_LUMO[cell].append(float(moenergies_array[data.homos + 1]))
        descriptions.append(description)
        # get bond structure data
        pos = atoms.get_positions()
        tree = KDTree(atoms.get_positions())
        # find pairs of atoms that are bonded; they are bonded if their distance apart is less than 1.430
        list_tree = list(tree.query_pairs(1.430))
        bondedTo = [[] for i in range(len(atoms))]
        for bond in list_tree:
            bondedTo[bond[0]].append(bond[1])
            bondedTo[bond[1]].append(bond[0])
        Zs = atoms.get_atomic_numbers()
        # find N atoms and their x and z positions
        Ns = edge_carbon_index
        xN = [pos[i][0] for i in edge_carbon_index]
        uxN = np.unique(xN)
        # sort into pairs of atoms at the same x position
        sorted_Ns = [[] for i in range(len(uxN))]
        for N in Ns:
            sorted_Ns[list(uxN).index(pos[N][0])].append(N)

        for step,position in enumerate(order):
            for k in range(0,cell+1):
                # hydrogenate one nitrogen
                iatom = sorted_Ns[position+k*cell*nx][0]
                r0 = pos[iatom]
                bond1 = pos[bondedTo[iatom][0]] - r0
                bond2 = pos[bondedTo[iatom][1]] - r0
                rH = -(bond1+bond2)
                rH = 1.09 * rH / np.linalg.norm(rH)
                atoms.append(Atom('H', r0 + rH))

                atoms_temp = atoms.copy()

                # find all C atoms
                Cs = []
                for iZ, Z in enumerate(Zs):
                    if Z == 6:
                        Cs.append(iZ)

                # pick one at random and add a hydrogen
                C = np.random.choice(Cs)
                CH_bond = [0,1.09,0]
                atoms_temp.append(Atom('H', pos[C]+CH_bond))

                atoms1 = atoms_temp.copy()
                atoms1.rotate("x", np.pi / 2.0)
                description = "nz" + str(nz) + "cell" + str(cell + 1) + "step" + str(step + 1) + "C"
                write(description + ".png", atoms1)
                data = make_orca(atoms_temp, filename=description + ".inp", multiplicity="1", method=method,
                                 geometry_opt=optimize_geometry,
                                 output=ORCA_filepath + "/hydrogenate_inner_carbon/" + description + ".out")
                moenergies_array = data.moenergies[0]
                C_SCF[cell].append(data.scfenergies[len(data.scfenergies) - 1])
                C_HOMO[cell].append(float(moenergies_array[data.homos]))
                C_LUMO[cell].append(float(moenergies_array[data.homos + 1]))
                descriptions.append(description)

                with open("results.txt", "a+") as e:
                    e.write("\n##########################\n")
                    e.write(description)
                    e.write("\n##########################\n")
                    e.write("Total SCF energy in eV:\t")
                    e.write(str(data.scfenergies))
                    e.write("\nMolecular orbital energy of HOMO in eV:\t")
                    e.write(str(moenergies_array[data.homos]))
                    e.write("\nMolecular orbital energy of LUMO in eV:\t")
                    e.write(str(moenergies_array[data.homos + 1]))
                    e.flush()

                iatom = sorted_Ns[position+k*cell*nx][1]
                r0 = pos[iatom]
                bond1 = pos[bondedTo[iatom][0]] - r0
                bond2 = pos[bondedTo[iatom][1]] - r0
                rH = -(bond1 + bond2)
                rH = 1.09 * rH / np.linalg.norm(rH)
                atoms.append(Atom('H', r0 + rH))

                atoms1 = atoms.copy()
                atoms1.rotate("x", np.pi / 2.0)
                description = "nz" + str(nz) + "cell" + str(cell + 1) + "step" + str(step + 1) + "N"
                write(description + ".png", atoms1)
                data = make_orca(atoms, filename=description + ".inp", multiplicity="1", method=method,
                                 geometry_opt=optimize_geometry,
                                 output=ORCA_filepath + "/hydrogenate_inner_carbon/" + description + ".out")
                moenergies_array = data.moenergies[0]
                N_SCF[cell].append(data.scfenergies[len(data.scfenergies) - 1])
                N_HOMO[cell].append(float(moenergies_array[data.homos]))
                N_LUMO[cell].append(float(moenergies_array[data.homos + 1]))
                descriptions.append(description)

                with open("results.txt", "a+") as e:
                    e.write("\n##########################\n")
                    e.write(description)
                    e.write("\n##########################\n")
                    e.write("Total SCF energy in eV:\t")
                    e.write(str(data.scfenergies))
                    e.write("\nMolecular orbital energy of HOMO in eV:\t")
                    e.write(str(moenergies_array[data.homos]))
                    e.write("\nMolecular orbital energy of LUMO in eV:\t")
                    e.write(str(moenergies_array[data.homos + 1]))
                    e.flush()

    cleanup(ORCA_filepath + "/hydrogenate_inner_carbon/")

    N_SCF_kcal = np.array(N_SCF) * EV_TO_KCAL
    C_SCF_kcal = np.array(C_SCF) * EV_TO_KCAL
    delta_E_N = []
    delta_E_C = []

    for cell in range(0,unit_cells):
        delta_E_N.append([])
        delta_E_C.append([])
        delta_E_N[cell] = [N_SCF_kcal[cell][i+1]-(N_SCF_kcal[cell][i]+H2_TOTAL_ENERGY) for i in range(0,len(N_SCF_kcal[cell])-1)]
        delta_E_C[cell] = [C_SCF_kcal[cell][i]-(N_SCF_kcal[cell][i]+H2_TOTAL_ENERGY) for i in range(0,len(C_SCF_kcal[cell]))]

        fig = plt.figure()
        plt.plot([(i+1)*(cell+1) for i in range(0,len(delta_E_N[cell]))], delta_E_N[cell], 'bo-')

        plt.plot([1], [delta_E_C[cell][0]], 'ro-')
        for i in range(0,len(C_SCF_kcal[cell])-1):
            plt.plot([(i+1)*(cell+1),(i+2)*(cell+1)], [delta_E_N[cell][i],delta_E_C[cell][i+1]], 'ro-')

        plt.xticks(np.arange(cell+1,(cell+1)*len(delta_E_N[cell])+1,cell+1))
        plt.savefig("rxn energy.jpg")

def hydrogenate_random(nx=5, nz=3, method="am1", optimize_geometry=0, iter=5, directory="/hydrogenate_random"):
    os.popen("mkdir " + ORCA_filepath + directory)
    os.chdir(ORCA_filepath + directory)

    i = 0
    while i < iter:
        atoms = prepare_nitrogenated_sheet(nx,nz)
        Ns = get_outer_zigzag_atoms(nx,nz)
        daves_super_saturate(atoms)

        # get bond structure data
        pos = atoms.get_positions()
        tree = KDTree(atoms.get_positions())
        # find pairs of atoms that are bonded; they are bonded if their distance apart is less than 1.430
        list_tree = list(tree.query_pairs(1.430))
        bondedTo = [[] for k in range(len(atoms))]
        for bond in list_tree:
            bondedTo[bond[0]].append(bond[1])
            bondedTo[bond[1]].append(bond[0])
        Zs = atoms.get_atomic_numbers()

        # iterate over N's, randomly decide to add H
        n_nH = 0
        for iatom in Ns:
            if np.random.random() > 0.5:
                r0 = pos[iatom]
                bond1 = pos[bondedTo[iatom][0]] - r0
                bond2 = pos[bondedTo[iatom][1]] - r0
                rH = -(bond1 + bond2)
                rH = 1.09 * rH / np.linalg.norm(rH)
                atoms.append(Atom('H', r0 + rH))
                n_nH += 1

        # add an additional H to a random C if an odd number was added to N's
        if n_nH%2 == 1:
            # find all C atoms
            Cs = []
            for iZ, Z in enumerate(Zs):
                if Z == 6:
                    Cs.append(iZ)

            # pick one at random and add a hydrogen
            C = np.random.choice(Cs)
            CH_bond = [0, 1.09, 0]
            atoms.append(Atom('H', pos[C] + CH_bond))

            # data = make_orca(atoms, filename="sheet"+str(i+1)+".inp", multiplicity="1", method=method,
            #                  geometry_opt=optimize_geometry,
            #                  output=ORCA_filepath+"/hydrogenate_random/sheet"+str(i+1)+".out")
            # cleanup(ORCA_filepath + "/hydrogenate_random/")
            #
            # if optimize_geometry == 0:
            #     xyz = open("step" + str(i + 1) + ".xyz", "w")
            #     xyz.write(str(len(atoms)) + "\n")
            #     xyz.write("comment\n")
            #     symbols_temp = atoms.get_chemical_symbols()
            #     pos_temp = atoms.get_positions()
            #     for iatom in range(0, len(atoms)):
            #         xyz.write(
            #             symbols_temp[iatom] + "\t" + str(pos_temp[iatom][0]) + "\t" + str(pos_temp[iatom][1]) + "\t" + str(
            #                 pos_temp[iatom][2]) + "\n")
            #     xyz.close()

            i += 1

        # create an image of the sheet
        atoms1 = atoms.copy()
        atoms1.rotate("x", np.pi / 2.0)
        write(str(i) + ".png", atoms1)


def sh_vary_width(unit_cells=2, nx=5, nz_min=3, nz_max=9, method="am1", optimize_geometry=0, make_symmetric=0):

    data = []

    for nz in range(nz_min,nz_max+1,2):

        descriptions, scf_list, homo_list, lumo_list = selectively_hydrogenate_nitrogens(unit_cells,nx,nz,method,optimize_geometry,make_symmetric)
        scf = np.array(scf_list)
        scf_kcal = scf * EV_TO_KCAL
        sheet_nH2 = []
        ref_to_0 = []
        for cell in range(0, unit_cells):
            sheet_nH2.append([])
            ref_to_0.append([])
            sheet_nH2[cell] = [scf_kcal[cell][i] + (nx - i - 1) * (cell + 1) * H2_TOTAL_ENERGY for i in
                               range(0, len(scf_kcal[cell]))]
            ref_to_0[cell] = [sheet_nH2[cell][i] - sheet_nH2[cell][0] for i in range(0, len(sheet_nH2[cell]))]

        data.append({'nz':nz, 'ref_to_0':ref_to_0})

    colors = ['r','y','g','b']
    for cell in range(0, unit_cells):
        fig = plt.figure()
        for i in range(0, len(data)):
            plt.plot([k*(cell+1) for k in range(0,len(data[i]['ref_to_0'][cell]))], data[i]['ref_to_0'][cell], colors[i]+'o-', label="Width = "+str(data[i]['nz']))
        plt.xticks(np.arange(0,(cell+1)*len(data[0]['ref_to_0'][cell]),cell+1))
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel("n H2 adsorbed")
        ax.set_ylabel("Energy relative to sheet before hydrogenation (kcal/mol)")
        plt.legend(loc="best")
        plt.savefig("cell"+str(cell+1)+" ref to 0.jpg")

def plot_energy(all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list, unit_cells, nx, nz):

    SCF_kcal = np.array(all_N_SCF_energy_list) * EV_TO_KCAL
    delta_E_sheet = []
    delta_E_rxn = []
    sheet_nH2 = []
    ref_to_0 = []
    E_plus2_species = []
    H_attachment = []
    Hplus_attachment = []
    for cell in range(0,unit_cells):
        delta_E_sheet.append([])
        delta_E_rxn.append([])
        sheet_nH2.append([])
        ref_to_0.append([])
        E_plus2_species.append([])
        H_attachment.append([])
        Hplus_attachment.append([])
        delta_E_sheet[cell] = [SCF_kcal[cell][i+1]-SCF_kcal[cell][i] for i in range(0,len(SCF_kcal[cell])-1)]
        delta_E_rxn[cell] = [SCF_kcal[cell][i+1]-(SCF_kcal[cell][i]+H2_TOTAL_ENERGY) for i in range(0,len(SCF_kcal[cell])-1)]
        sheet_nH2[cell] = [SCF_kcal[cell][i]+(nx-i-1)*(cell+1)*H2_TOTAL_ENERGY for i in range(0,len(SCF_kcal[cell]))]
        ref_to_0[cell] = [sheet_nH2[cell][i]-sheet_nH2[cell][0] for i in range(0,len(sheet_nH2[cell]))]
        E_plus2_species[cell] = [all_N_SCF_energy_list[cell][i]-2*cell*all_N_HOMO_energy_list[cell][i] for i in range(0,len(all_N_SCF_energy_list[cell]))]
        H_attachment[cell] = [(all_N_SCF_energy_list[cell][i+1]-all_N_SCF_energy_list[cell][i])/(2*(cell+1)) for i in range(0,len(all_N_SCF_energy_list[cell])-1)]
        Hplus_attachment[cell] = [(E_plus2_species[cell][i+1]-E_plus2_species[cell][i])/(2*(cell+1)) for i in range(0,len(E_plus2_species[cell])-1)]

        fig1 = plt.figure()
        plt.plot([(i+1)*(cell+1) for i in range(0,len(delta_E_rxn[cell]))], delta_E_rxn[cell], 'bo-')
        plt.xticks(np.arange(cell+1,(cell+1)*len(delta_E_rxn[cell])+1,cell+1))
        ax1 = fig1.add_subplot(1,1,1)
        ax1.set_xlabel("n H2 adsorbed")
        ax1.set_ylabel("Reaction energy (kcal/mol)")
        plt.savefig("cell"+str(cell+1)+" rxn energy.jpg")

        fig2 = plt.figure()
        plt.plot([i*(cell+1) for i in range(0,len(ref_to_0[cell]))], ref_to_0[cell], 'bo-')
        plt.xticks(np.arange(0,(cell+1)*len(ref_to_0[cell]),cell+1))
        ax2 = fig2.add_subplot(1,1,1)
        ax2.set_xlabel("n H2 adsorbed")
        ax2.set_ylabel("Energy relative to sheet before hydrogenation (kcal/mol)")
        plt.savefig("cell"+str(cell+1)+" ref to 0.jpg")

        fig3 = plt.figure()
        plt.plot([i*2*(cell+1) for i in range(0,len(all_N_HOMO_energy_list[cell]))], all_N_HOMO_energy_list[cell], "bo-", label="HOMO")
        plt.plot([i*2*(cell+1) for i in range(0,len(all_N_LUMO_energy_list[cell]))], all_N_LUMO_energy_list[cell], "ro-", label="LUMO")
        plt.plot([(i+1)*2*(cell+1) for i in range(0,len(H_attachment[cell]))], H_attachment[cell], "go-", label="H attachment")
        plt.xticks(np.arange(2*(cell+1), 2*(cell+1)*len(all_N_HOMO_energy_list[cell]),2*(cell+1)))
        ax3 = fig3.add_subplot(1,1,1)
        ax3.set_xlabel("n H atoms adsorbed")
        ax3.set_ylabel("eV")
        plt.legend(loc="best")
        plt.savefig("cell"+str(cell+1)+" HOMO LUMO.jpg")

def NH_combinations(nx=5, nz=3, method="am1", optimize_geometry=0, make_symmetric=1):
    binaries = []
    scf = []
    homo = []
    lumo = []
    os.popen("mkdir " + ORCA_filepath + "/NH_combinations")
    os.chdir(ORCA_filepath + "/NH_combinations")
    results = open("results.txt", "w")
    results.truncate()
    # find possible combinations, expressed in binary
    # 0 = no H, 1 = H added
    for i in range(0,2**nx):
        binary = bin(i)[2:]
        while len(binary) < nx:
            binary = '0'+binary
        binaries.append(binary)

    combinations = []
    for binary1 in binaries:
        for binary2 in binaries:
            combo = [binary1,binary2]
            mirror = [combo[0][::-1],combo[1][::-1]]
            if (not mirror in combinations) and (binary1.count('1')+binary2.count('1'))%2 == 0:
                combinations.append(combo)

    # iterate over combinations
    for combo in combinations:

        atoms = prepare_nitrogenated_sheet(nx,nz)
        edge_carbon_index = get_outer_zigzag_atoms(nx,nz)

        front_Ns = []
        back_Ns = []
        pos = atoms.get_positions()
        z = [pos[k][2] for k in range(0,len(pos))]
        zfront = np.amax(z)
        zback = np.amin(z)

        # iterate over outer C's on zigzag
        for edge_carbon in edge_carbon_index:
            # replace with N
            symbols = atoms.get_chemical_symbols()
            symbols[edge_carbon] = 'N'
            atoms.set_chemical_symbols(symbols)
            # determine which edge
            if abs(z[edge_carbon]-zfront) < abs(z[edge_carbon]-zback):
                front_Ns.append(edge_carbon)
            else:
                back_Ns.append(edge_carbon)

        # sort N's by x position
        front_x = np.sort([pos[N][0] for N in front_Ns])
        sorted_front_Ns = [0 for k in range(0,len(front_Ns))]
        for N in front_Ns:
            sorted_front_Ns[list(front_x).index(pos[N][0])] = N
        back_x = np.sort([pos[N][0] for N in back_Ns])
        sorted_back_Ns = [0 for k in range(0,len(back_Ns))]
        for N in back_Ns:
            sorted_back_Ns[list(back_x).index(pos[N][0])] = N

        # get bond structure data
        tree = KDTree(pos)
        # find bonded pairs
        list_tree = list(tree.query_pairs(1.430))
        bondedTo = [[] for k in range(len(atoms))]
        for bond in list_tree:
            bondedTo[bond[0]].append(bond[1])
            bondedTo[bond[1]].append(bond[0])

        nH = 0
        # iterate over N's
        for iN,N in enumerate(sorted_front_Ns):
            if combo[0][iN] == '1':
                # add H
                r0 = pos[N]
                bond1 = pos[bondedTo[N][0]] - r0
                bond2 = pos[bondedTo[N][1]] - r0
                rH = -(bond1+bond2)
                rH = 1.09 * rH / np.linalg.norm(rH)
                atoms.append(Atom('H', r0+rH))
                nH = nH+1
        for iN,N in enumerate(sorted_back_Ns):
            if combo[1][iN] == '1':
                r0 = pos[N]
                bond1 = pos[bondedTo[N][0]] - r0
                bond2 = pos[bondedTo[N][1]] - r0
                rH = -(bond1+bond2)
                rH = 1.09 * rH / np.linalg.norm(rH)
                atoms.append(Atom('H', r0+rH))
                nH = nH+1

        # saturate C's on armchair edges
        daves_super_saturate(atoms)

        # create an image of the sheet
        atoms1 = atoms.copy()
        atoms1.rotate("x", np.pi/2.0)
        write(combo[0]+"_"+combo[1]+".png", atoms1)

        # # run Orca
        # mult = "1" if nH%2==0 else "2"
        # data = make_orca(atoms, filename=combo[0]+"_"+combo[1]+".inp",multiplicity=mult, method=method, geometry_opt=optimize_geometry, output=ORCA_filepath+"/NH_combinations/"+combo[0]+"_"+combo[1]+".out")
        #
        # try:
        #     moenergies_array = data.moenergies[0]
        #     scf.append(data.scfenergies)
        #     homo.append(moenergies_array[data.homos])
        #     lumo.append(moenergies_array[data.homos+1])
        #
        #     # write results
        #     results.write("##########################\n")
        #     results.write("Hydrogenation pattern:\t"+combo[0]+"_"+combo[1]+"\n")
        #     results.write("##########################\n")
        #     results.write("Total SCF energy in eV:\t"+str(data.scfenergies)+"\n")
        #     results.write("Molecular orbital energy of HOMO in eV:\t"+str(moenergies_array[data.homos])+"\n")
        #     results.write("Molecular orbital energy of LUMO in eV:\t"+str(moenergies_array[data.homos+1])+"\n")
        # except AttributeError:
        #     results.write("##########################\n")
        #     results.write("Hydrogenation pattern:\t"+combo[0]+"_"+combo[1]+"\n")
        #     results.write("##########################\n")
        #     results.write("SCF NOT CONVERGED\n")

    cleanup(ORCA_filepath+"/NH_combinations")

def NH_combinations_plot():
    os.chdir(ORCA_filepath + "/NH_combinations")
    results = open("results.txt", "r")
    loc = 0
    binaries = []
    scf = []
    homo = []
    lumo = []
    binary = ''
    for line in results:
        if loc == 0:
            if re.search("Hydrogenation pattern", line):
                words = line.rstrip().split()
                binary = words[len(words)-1]
                loc = 1
        elif loc == 1:
            if re.search("Total SCF energy in eV", line):
                words = line.rstrip().split()
                scf.append(float(words[len(words)-1].strip('[').strip(']')))
                binaries.append(binary)
                #loc = 2
                loc = 0
            elif re.search("SCF NOT CONVERGED", line):
                loc = 0
        # elif loc == 2:
        #     if re.search("Molecular orbital energy of HOMO in eV", line):
        #         words = line.rstrip().split()
        #         homo.append(float(words[len(words)-1].strip('[').strip(']')))
        #         loc = 3
        # elif loc == 3:
        #     if re.search("Molecular orbital energy of LUMO in eV", line):
        #         words = line.rstrip().split()
        #         lumo.append(float(words[len(words)-1].strip('[').strip(']')))
        #         loc = 0

    nH = [binary.count('1') for binary in binaries]
    # scf = list(np.array(scf_eV)*EV_TO_KCAL)
    # sheet_nH = [scf[i]+(len(binaries[0])-nH[i])*H2_TOTAL_ENERGY/2 for i in range(0,len(scf))]
    # ref_to_0 = [sheet_nH[i]-sheet_nH[0] for i in range(0,len(sheet_nH))]

    fig1 = plt.figure()
    plt.plot(nH, scf, 'o')
    ax1 = fig1.add_subplot(1,1,1)
    ax1.set_xlabel("n H atoms adsorbed")
    ax1.set_ylabel("Total SCF energy in eV")
    # for i,binary in enumerate(binaries):
#     ax.text(binary.count('1'), scf[i], binary)
    plt.savefig("scf.jpg")

    # fig2 = plt.figure()
    # plt.plot(nH, homo, 'o')
    # ax2 = fig2.add_subplot(1,1,1)
    # ax2.set_xlabel("n H atoms adsorbed")
    # ax2.set_ylabel("Molecular orbital energy of HOMO in eV")
    # plt.savefig("homo.jpg")
    #
    # fig3 = plt.figure()
    # plt.plot(nH, lumo, 'o')
    # ax3 = fig3.add_subplot(1,1,1)
    # ax3.set_xlabel("n H atoms adsorbed")
    # ax3.set_ylabel("Molecular orbital energy of LUMO in eV")
    # plt.savefig("lumo.jpg")

    # sort configurations by number of hydrogens
    sorted_binaries = [[] for i in range(0,len(binaries[0])+1)]
    sorted_scf = [[] for i in range(0,len(binaries[0])+1)]
    for i,n in enumerate(nH):
        sorted_binaries[n].append(binaries[i])
        sorted_scf[n].append(scf[i])
    # find the configuration with the minimum energy for each number of hydrogens
    min_binaries = []
    min_scf = []
    for i in range(0,len(sorted_binaries)):
        min_scf.append(np.amin(sorted_scf[i]))
        min_binaries.append(sorted_binaries[i][np.argmin(sorted_scf[i])])
    # iterate over pH and potential
    pHs = [i for i in np.linspace(-2,16,200)]
    potentials = [i for i in np.linspace(-1.8,2.2,200)]
    min_nH = []
    for ipH,pH in enumerate(pHs):
        min_nH.append([])
        for ipotential,potential in enumerate(potentials):
            # find the number of hydrogens for the most stable sheet at these conditions
            min_nH[ipH].append(np.argmin([min_scf[n] - min_scf[0] - n*H2_TOTAL_ENERGY/(2.0*EV_TO_KCAL) + n*0.059*pH
                                                  + n*potential for n in range(0,len(min_scf))]))
    # create a Pourbaix-esque diagram showing the regions where different sheets are most stable
    fig4 = plt.figure()
    cmap = plt.get_cmap("rainbow", np.amax(min_nH)-np.amin(min_nH)+1)
    pourbaix = plt.pcolor(np.array(pHs), np.array(potentials), np.array(min_nH), cmap=cmap, vmin=np.amin(min_nH)-0.5,
                          vmax=np.amax(min_nH)+0.5)
    plt.xticks(np.arange(-2,16,1))
    plt.yticks(np.arange(-1.8,2.2,0.2))
    ax4 = fig4.add_subplot(1,1,1)
    ax4.set_xlabel("pH")
    ax4.set_ylabel("Potential (eV)")
    cb4 = plt.colorbar(pourbaix, ticks=np.arange(np.amin(min_nH),np.amax(min_nH)+1,1))
    cb4.set_label("Number of adsorbed hydrogens")
    plt.savefig("pourbaix.jpg")

def plot_orbitals(data_path):
    os.chdir(data_path)
    os.popen("mkdir orbitals")
    os.chdir("orbitals")
    path_contents = os.listdir(data_path)
    all_eorbs = []
    all_efermi = []
    all_coeffs = []
    all_homos = []
    all_pi_filled = []
    all_pi_empty = []
    all_sig_filled = []
    all_sig_empty = []
    E = np.arange(-50,10,0.01)
    for filename in path_contents:
        if filename[-4:] == ".out":
            data = parse(data_path+filename)
            try:
                eorbs = data.moenergies[0]
                efermi = (eorbs[data.homos[0]] + eorbs[data.homos[0]+1]) / 2
                print(filename+' '+str(efermi))
                # sets width of the gaussian we will put on each orbital energy
                width = 1
                # E = np.arange(np.amin(eorbs)-3.0*width, np.amax(eorbs)+3.0*width, 0.01)

                # rho_filled = np.zeros(E.shape)
                # for eorb in [e for e in eorbs if e <= efermi]:
                #     rho_filled += np.exp(-np.power(np.subtract(E, eorb), 2) / (width ** 2))
                #
                # rho_empty = np.zeros(E.shape)
                # for eorb in [e for e in eorbs if e > efermi]:
                #     rho_empty += np.exp(-np.power(np.subtract(E, eorb), 2) / (width ** 2))

                orbnames = data.aonames
                pyorb = None
                for iname,name in enumerate(orbnames):
                    if name[-2:] == "PY":
                        pyorb = iname
                        break

                coeffs = data.mocoeffs[0]

                plt.figure()
                pi_filled = np.zeros(E.shape)
                pi_empty = np.zeros(E.shape)
                sig_filled = np.zeros(E.shape)
                sig_empty = np.zeros(E.shape)
                for i,eorb in enumerate(eorbs):
                    if eorb <= efermi:
                        if coeffs[i][pyorb] == 0:
                            plt.plot([0, 2],[eorb,eorb], 'b-', alpha=0.5)
                            sig_filled += np.exp(-np.power(np.subtract(E, eorb), 2) / (width ** 2))
                        else:
                            plt.plot([0, 2], [eorb, eorb], 'b-')
                            pi_filled += np.exp(-np.power(np.subtract(E, eorb), 2) / (width ** 2))
                    else:
                        if coeffs[i][pyorb] == 0:
                            plt.plot([0, 2], [eorb,eorb], 'r-', alpha=0.5)
                            sig_empty += np.exp(-np.power(np.subtract(E, eorb), 2) / (width ** 2))
                        else:
                            plt.plot([0, 2], [eorb, eorb], 'r-')
                            pi_empty += np.exp(-np.power(np.subtract(E, eorb), 2) / (width ** 2))

                plt.plot(pi_filled, E, 'b-', label="Filled pi orbitals")
                plt.plot(pi_empty, E, 'r-', label="Empty pi orbitals")
                plt.plot(sig_filled, E, 'b-', alpha=0.5, label="Filled sigma orbitals")
                plt.plot(sig_empty, E, 'r-', alpha=0.5, label="Empty sigma orbitals")
                plt.legend(loc="best")
                plt.savefig(filename[:-4]+".jpg")

                all_eorbs.append(eorbs)
                all_efermi.append(efermi)
                all_coeffs.append(coeffs)
                all_homos.append(data.homos[0])
                all_pi_filled.append(pi_filled)
                all_pi_empty.append(pi_empty)
                all_sig_filled.append(sig_filled)
                all_sig_empty.append(sig_empty)
            except AttributeError:
                print("AttributeError")

    for i in range(0, len(all_coeffs)-1):

        plt.figure()
        plt.plot(np.subtract(all_pi_filled[i+1],all_pi_filled[i]), E, 'b-', label="Filled pi orbitals")
        plt.plot(np.subtract(all_pi_empty[i+1],all_pi_empty[i]), E, 'r-', label="Empty pi orbitals")
        plt.plot(np.subtract(all_sig_filled[i+1],all_sig_filled[i]), E, 'b-', alpha=0.5, label="Filled sigma orbitals")
        plt.plot(np.subtract(all_sig_empty[i+1],all_sig_empty[i]), E, 'r-', alpha=0.5, label="Empty sigma orbitals")
        for iorb, eorb in enumerate(all_eorbs[0]):
            if eorb <= all_efermi[0]:
                if all_coeffs[0][iorb][pyorb] == 0:
                    plt.plot([-10,-8], [eorb, eorb], 'b-', alpha=0.5)
                else:
                    plt.plot([-10,-8], [eorb, eorb], 'b-')
            else:
                if all_coeffs[0][iorb][pyorb] == 0:
                    plt.plot([-10,-8], [eorb, eorb], 'r-', alpha=0.5)
                else:
                    plt.plot([-10,-8], [eorb, eorb], 'r-')

        for iorb, eorb in enumerate(all_eorbs[i + 1]):
            if eorb <= all_efermi[i + 1]:
                if all_coeffs[i + 1][iorb][pyorb] == 0:
                    plt.plot([8,10], [eorb, eorb], 'b-', alpha=0.5)
                else:
                    plt.plot([8,10], [eorb, eorb], 'b-')
            else:
                if all_coeffs[i + 1][iorb][pyorb] == 0:
                    plt.plot([8,10], [eorb, eorb], 'r-', alpha=0.5)
                else:
                    plt.plot([8,10], [eorb, eorb], 'r-')
        plt.legend(loc="best")
        plt.savefig("delta"+str(i+1)+".jpg")

        if len(all_coeffs[0]) < len(all_coeffs[i+1]):
            diff = len(all_coeffs[i+1]) - len(all_coeffs[0])
            all_coeffs[0] = np.vstack((all_coeffs[0], np.zeros((diff,len(all_coeffs[0])))))
            all_coeffs[0] = np.column_stack((all_coeffs[0], np.zeros((len(all_coeffs[0]),diff))))
        elif len(all_coeffs[0]) > len(all_coeffs[i+1]):
            diff = len(all_coeffs[0]) - len(all_coeffs[i+1])
            all_coeffs[i+1] = np.vstack((all_coeffs[i+1], np.zeros((diff,len(all_coeffs[i+1])))))
            all_coeffs[i+1] = np.column_stack((all_coeffs[i+1], np.zeros((len(all_coeffs[i+1]),diff))))

        connections = np.dot(all_coeffs[0], all_coeffs[i+1].transpose()) ** 2

        for k in range(0,10):

            plt.figure()

            for iorb, eorb in enumerate(all_eorbs[0]):
                if eorb <= all_efermi[0]:
                    if all_coeffs[0][iorb][pyorb] == 0:
                        plt.plot([0, 2], [eorb, eorb], 'b-', alpha=0.3, label="filled")
                    else:
                        plt.plot([0, 2], [eorb, eorb], 'b-', label="filled")
                else:
                    if all_coeffs[0][iorb][pyorb] == 0:
                        plt.plot([0, 2], [eorb, eorb], 'r-', alpha=0.1, label="empty")
                    else:
                        plt.plot([0, 2], [eorb, eorb], 'r-', label="empty")

            for iorb, eorb in enumerate(all_eorbs[i + 1]):
                if eorb <= all_efermi[i + 1]:
                    if all_coeffs[i + 1][iorb][pyorb] == 0:
                        plt.plot([4, 6], [eorb, eorb], 'b-', alpha=0.3, label="filled")
                    else:
                        plt.plot([4, 6], [eorb, eorb], 'b-', label="filled")
                else:
                    if all_coeffs[i + 1][iorb][pyorb] == 0:
                        plt.plot([4, 6], [eorb, eorb], 'r-', alpha=0.1, label="empty")
                    else:
                        plt.plot([4, 6], [eorb, eorb], 'r-', label="empty")

            iorb = all_homos[i + 1] - 5 + k
            maxcon = np.amax([connections[m][iorb] for m in range(0,len(connections))])
            for m in range(0,len(connections)-2):
                ratio = connections[m][iorb]/maxcon
                if ratio > 0.1:
                    plt.plot([2,4], [all_eorbs[0][m], all_eorbs[i+1][iorb]], 'k-')
                elif ratio > 0.01:
                    plt.plot([2,4], [all_eorbs[0][m], all_eorbs[i+1][iorb]], 'k-', alpha=0.2)

            plt.title("Sheet -> SheetH_"+str(2*(i+1)))
            plt.savefig("step"+str(i+1)+"orb"+str(k+1)+".jpg")
            plt.close()

def nitrogenate_inner_zigzag(nx=5, nz=3, method="am1", optimize_geometry=0, make_symmetric=1):
    all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list, descriptions = ([] for dummy_var in range(4))

    order = []
    for i in range(0, nx/2):
        order.append(nx/2-i-1)
        order.append(nx-i-2)

    os.popen("mkdir " + ORCA_filepath + "/nitrogenate_inner_zigzag")
    os.chdir(ORCA_filepath + "/nitrogenate_inner_zigzag")

    atoms = build_sheet(nx, nz, symmetry=1)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for i, atom in enumerate(atoms):
        plt.plot(atom.position[0], atom.position[2], '.')
        plt.text(atom.position[0], atom.position[2], str(i))
    plt.savefig("map.png")

    #find the inner atoms of the zigzag edges
    inner_front = []
    inner_back = []
    for i in range(0, nx - 1):
        inner_front.append((nz + 2) * (i + 1) + nz * i)
        inner_back.append((nz + 2) * (i + 1) + nz * (i + 1) + nz - 1)

    for step,position in enumerate(order):
        symbols = atoms.get_chemical_symbols()
        symbols[inner_front[position]] = 'N'
        symbols[inner_back[position]] = 'N'
        atoms.set_chemical_symbols(symbols)

        daves_super_saturate(atoms)

        atoms1 = atoms.copy()
        atoms1.rotate("x", np.pi/2.0)
        write("sheet"+str(step)+".png", atoms1)

        data = make_orca(atoms, filename="sheet"+str(step)+".inp", output="sheet"+str(step)+".out", geometry_opt=optimize_geometry)

        try:
            moenergies_array = data.moenergies[0]

            with open("results.txt", "a+") as e:
                e.write("\n##########################\n")
                e.write("sheet"+str(step))
                e.write("\n##########################\n")
                e.write("Total SCF energy in eV:\t")
                e.write(str(data.scfenergies))
                e.write("\nMolecular orbital energy of HOMO in eV:\t")
                e.write(str(moenergies_array[data.homos]))
                e.write("\nMolecular orbital energy of LUMO in eV:\t")
                e.write(str(moenergies_array[data.homos + 1]))
        except AttributeError:
            with open("results.txt", "a+") as e:
                e.write("\n##########################\n")
                e.write("sheet"+str(step))
                e.write("\n##########################\n")
                e.write("SCF failed")

    plot_orbitals(ORCA_filepath+"/nitrogenate_inner_zigzag/")

# def form_nitroso(nx=5, nz=3, method="am1", optimize_geometry=0, make_symmetric=1):
#     all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list, descriptions = ([[]] for dummy_var in range(4))
#     cell = 0
#
#     order = []
#     for i in range(0, nx/2):
#         order.append(nx/2-i-1)
#         order.append(nx-i-1)
#
#     atoms = build_sheet(nx, nz, symmetry=1)
#     edge_carbon_index = get_outer_zigzag_atoms(nx,nz)
#
#     # replace outer C's on zigzag with N's
#     for edge_carbon in edge_carbon_index:
#         symbols = atoms.get_chemical_symbols()
#         symbols[edge_carbon] = 'N'
#         atoms.set_chemical_symbols(symbols)
#     # saturate C's on armchair edges
#     daves_super_saturate(atoms)
#     # calculate energy of initial sheet
#     atoms1 = atoms.copy()
#     atoms1.rotate("x", np.pi / 2.0)
#     description = "nz" + str(nz) + "cell" + str(cell + 1) + "step0"
#     os.popen("mkdir " + ORCA_filepath + "/nitroso")
#     os.chdir(ORCA_filepath + "/nitroso")
#     write(description + ".png", atoms1)
#     data = make_orca(atoms, filename=description + ".inp", multiplicity="1", method=method,
#                      geometry_opt=optimize_geometry,
#                      output=ORCA_filepath + "/nitroso/" + description + ".out")
#     moenergies_array = data.moenergies[0]
#     all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
#     all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos]))
#     all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos + 1]))
#     descriptions.append(description)
#     # get bond structure data
#     pos = atoms.get_positions()
#     tree = KDTree(atoms.get_positions())
#     list_tree = list(tree.query_pairs(1.430))
#     bondedTo = [[] for i in range(len(atoms))]
#     for bond in list_tree:
#         bondedTo[bond[0]].append(bond[1])
#         bondedTo[bond[1]].append(bond[0])
#     Zs = atoms.get_atomic_numbers()
#     # find N atoms and their x and z positions
#     Ns = edge_carbon_index
#     xN = [pos[i][0] for i in edge_carbon_index]
#     uxN = np.unique(xN)
#     # sort into pairs of atoms at the same x position
#     sorted_Ns = [[] for i in range(len(uxN))]
#     for N in Ns:
#         sorted_Ns[list(uxN).index(pos[N][0])].append(N)
#     # bond O to N
#     for step,position in enumerate(order):
#         for k in range(0,cell+1):
#             for m in range(0,2):
#                 iatom = sorted_Ns[position+k*cell*nx][m]
#                 r0 = pos[iatom]
#                 bond1 = pos[bondedTo[iatom][0]] - r0
#                 bond2 = pos[bondedTo[iatom][1]] - r0
#                 rH = -(bond1 + bond2)
#                 rH = 1.5 * rH / np.linalg.norm(rH)
#                 # http://www.sciencedirect.com/science/article/pii/0022328X9186479A
#                 atoms.append(Atom('O', r0 + rH))
#         atoms1 = atoms.copy()
#         atoms1.rotate("x", np.pi/2.0)
#         description = "nz" + str(nz) + "cell" + str(cell + 1) + "step" + str(step + 1)
#         write(description + ".png", atoms1)
#         data = make_orca(atoms, filename=description + ".inp", multiplicity="1", method=method,
#                          geometry_opt=optimize_geometry,
#                          output=ORCA_filepath + "/nitroso/" + description + ".out")
#         moenergies_array = data.moenergies[0]
#         all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
#         all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos]))
#         all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos + 1]))
#         descriptions.append(description)
#
#         with open("results.txt", "a+") as e:
#             e.write("\n##########################\n")
#             e.write(description)
#             e.write("\n##########################\n")
#             e.write("Total SCF energy in eV:\t")
#             e.write(str(data.scfenergies))
#             e.write("\nMolecular orbital energy of HOMO in eV:\t")
#             e.write(str(moenergies_array[data.homos]))
#             e.write("\nMolecular orbital energy of LUMO in eV:\t")
#             e.write(str(moenergies_array[data.homos + 1]))
#
#     cleanup(ORCA_filepath+"/nitroso/")
#
#     plot_orbitals(ORCA_filepath+"/nitroso/")

# def form_ether(nx=5, nz=3, method="am1", optimize_geometry=0, make_symmetric=1):
#
#     cell = 0
#     all_N_SCF_energy_list = [[]]
#     all_N_HOMO_energy_list = [[]]
#     all_N_LUMO_energy_list = [[]]
#     descriptions = [[]]
#
#     # determine the order of positions to replace C with N
#     order = []
#     for i in range(0, nx/2):
#         order.append(nx/2-i-1)
#         order.append(nx-i-1)
#
#     atoms = build_sheet(nx, nz, symmetry=1)
#     edge_carbon_index = get_outer_zigzag_atoms(nx,nz)
#
#     atoms_temp = atoms.copy()
#     daves_super_saturate(atoms_temp)
#     atoms1 = atoms_temp.copy()
#     atoms1.rotate("x", np.pi/2.0)
#     description = "nz"+str(nz)+"cell"+str(cell+1)+"step0"
#     os.popen("mkdir " + ORCA_filepath + "/ether")
#     os.chdir(ORCA_filepath + "/ether")
#     write(description+".png", atoms1)
#     data = make_orca(atoms_temp, filename=description + ".inp", multiplicity="1", method=method,
#                      geometry_opt=optimize_geometry,
#                      output=ORCA_filepath + "/ether/" + description + ".out")
#     moenergies_array = data.moenergies[0]
#     all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
#     all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos]))
#     all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos + 1]))
#     descriptions.append(description)
#
#     # get bond structure data
#     pos = atoms.get_positions()
#     tree = KDTree(atoms.get_positions())
#     # find pairs of atoms that are bonded; they are bonded if their distance apart is less than 1.430
#     list_tree = list(tree.query_pairs(1.430))
#     bondedTo = [[] for i in range(len(atoms))]
#     for bond in list_tree:
#         bondedTo[bond[0]].append(bond[1])
#         bondedTo[bond[1]].append(bond[0])
#     Cs = edge_carbon_index
#     xC = [pos[i][0] for i in edge_carbon_index]
#     uxC = np.unique(xC)
#     # sort into pairs of atoms at the same x position
#     sorted_Cs = [[] for i in range(len(uxC))]
#     for C in Cs:
#         sorted_Cs[list(uxC).index(pos[C][0])].append(C)
#     # replace with O
#     for step,position in enumerate(order):
#         for k in range(0,cell+1):
#             for m in range(0,2):
#                 iatom = sorted_Cs[position+k*cell*nx][m]
#                 symbols = atoms.get_chemical_symbols()
#                 symbols[iatom] = 'O'
#                 atoms.set_chemical_symbols(symbols)
#                 atoms_temp = atoms.copy()
#             daves_super_saturate(atoms_temp)
#             atoms1 = atoms_temp.copy()
#             atoms1.rotate("x", np.pi/2.0)
#             description = "nz"+str(nz)+"cell"+str(cell+1)+"step"+str(step+1)
#             write(description+".png", atoms1)
#             data = make_orca(atoms_temp, filename=description+".inp", multiplicity="1", method=method, geometry_opt=optimize_geometry, output=ORCA_filepath+"/ether/"+description+".out")
#             moenergies_array = data.moenergies[0]
#             all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
#             all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos]))
#             all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos + 1]))
#             descriptions.append(description)
#
#             with open("results.txt", "a+") as e:
#                 e.write("\n##########################\n")
#                 e.write(description)
#                 e.write("\n##########################\n")
#                 e.write("Total SCF energy in eV:\t")
#                 e.write(str(data.scfenergies))
#                 e.write("\nMolecular orbital energy of HOMO in eV:\t")
#                 e.write(str(moenergies_array[data.homos]))
#                 e.write("\nMolecular orbital energy of LUMO in eV:\t")
#                 e.write(str(moenergies_array[data.homos + 1]))
#
#     cleanup(ORCA_filepath+"/ether/")
#
#     plot_orbitals(ORCA_filepath+"/ether/")

# def form_ketone(nx=5, nz=3, method="am1", optimize_geometry=0, make_symmetric=1):
#     all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list, descriptions = ([[]] for dummy_var in range(4))
#     cell = 0
#
#     order = []
#     for i in range(0, nx/2):
#         order.append(nx/2-i-1)
#         order.append(nx-i-1)
#
#     atoms = build_sheet(nx, nz, symmetry=1)
#     edge_carbon_index = get_outer_zigzag_atoms(nx,nz)
#     # calculate energy of initial sheet
#     atoms_temp = atoms.copy()
#     daves_super_saturate(atoms_temp)
#     atoms1 = atoms_temp.copy()
#     atoms1.rotate("x", np.pi / 2.0)
#     description = "nz" + str(nz) + "cell" + str(cell + 1) + "step0"
#     os.popen("mkdir " + ORCA_filepath + "/ketone")
#     os.chdir(ORCA_filepath + "/ketone")
#     write(description + ".png", atoms1)
#     data = make_orca(atoms_temp, filename=description + ".inp", multiplicity="1", method=method,
#                      geometry_opt=optimize_geometry,
#                      output=ORCA_filepath + "/ketone/" + description + ".out")
#     moenergies_array = data.moenergies[0]
#     all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
#     all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos]))
#     all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos + 1]))
#     descriptions.append(description)
#     # get bond structure data
#     pos = atoms.get_positions()
#     tree = KDTree(atoms.get_positions())
#     list_tree = list(tree.query_pairs(1.430))
#     bondedTo = [[] for i in range(len(atoms))]
#     for bond in list_tree:
#         bondedTo[bond[0]].append(bond[1])
#         bondedTo[bond[1]].append(bond[0])
#     Zs = atoms.get_atomic_numbers()
#     # find x and z positions of zigzag carbons
#     Cs = edge_carbon_index
#     xC = [pos[i][0] for i in edge_carbon_index]
#     uxC = np.unique(xC)
#     # sort into pairs of atoms at the same x position
#     sorted_Cs = [[] for i in range(len(uxC))]
#     for C in Cs:
#         sorted_Cs[list(uxC).index(pos[C][0])].append(C)
#     # bond O to C
#     for step,position in enumerate(order):
#         for k in range(0,cell+1):
#             for m in range(0,2):
#                 iatom = sorted_Cs[position+k*cell*nx][m]
#                 r0 = pos[iatom]
#                 bond1 = pos[bondedTo[iatom][0]] - r0
#                 bond2 = pos[bondedTo[iatom][1]] - r0
#                 rH = -(bond1 + bond2)
#                 rH = 1.2 * rH / np.linalg.norm(rH)
#                 atoms.append(Atom('O', r0 + rH))
#
#         atoms_temp = atoms.copy()
#
#         daves_super_saturate(atoms_temp)
#
#         atoms1 = atoms_temp.copy()
#         atoms1.rotate("x", np.pi/2.0)
#         description = "nz" + str(nz) + "cell" + str(cell + 1) + "step" + str(step + 1)
#         write(description+".png", atoms1)
#         data = make_orca(atoms_temp, filename=description+".inp", output=ORCA_filepath+"/ketone/"+description+".out", geometry_opt=optimize_geometry)
#         moenergies_array = data.moenergies[0]
#         all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
#         all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos]))
#         all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos + 1]))
#         descriptions.append(description)
#
#         with open("results.txt", "a+") as e:
#             e.write("\n##########################\n")
#             e.write(description)
#             e.write("\n##########################\n")
#             e.write("Total SCF energy in eV:\t")
#             e.write(str(data.scfenergies))
#             e.write("\nMolecular orbital energy of HOMO in eV:\t")
#             e.write(str(moenergies_array[data.homos]))
#             e.write("\nMolecular orbital energy of LUMO in eV:\t")
#             e.write(str(moenergies_array[data.homos + 1]))
#
#
#     cleanup(ORCA_filepath + "/ketone/")
#
#     plot_orbitals(ORCA_filepath + "/ketone/")

def form_alcohol(nx=5, nz=3, method="am1", optimize_geometry=0, zigzag=0, atoms=None):

    all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list, descriptions = ([[]] for dummy_var in range(4))
    cell = 0

    atoms = atoms or build_sheet(nx, nz, symmetry=1)
    edge_carbon_index = get_outer_zigzag_atoms(nx,nz)
    # calculate energy of initial sheet
    # atoms_temp = atoms.copy()
    # daves_super_saturate(atoms_temp)
    # atoms1 = atoms_temp.copy()
    # atoms1.rotate("x", np.pi / 2.0)
    # description = "nz" + str(nz) + "cell" + str(cell + 1) + "step0"
    directory = "/alcohol/" if zigzag == 0 else "/alcohol_zigzag/"
    os.popen("mkdir " + ORCA_filepath + directory)
    os.chdir(ORCA_filepath + directory)
    # write(description + ".png", atoms1)
    # data = make_orca(atoms_temp, filename=description + ".inp", multiplicity="1", method=method,
    #                  geometry_opt=optimize_geometry,
    #                  output=ORCA_filepath + directory + description + ".out")
    # moenergies_array = data.moenergies[0]
    # all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
    # all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos[0]]))
    # all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos[0] + 1]))
    # descriptions.append(description)
    # get bond structure data
    pos = atoms.get_positions()
    tree = KDTree(atoms.get_positions())
    list_tree = list(tree.query_pairs(1.430))
    bondedTo = [[] for i in range(len(atoms))]
    for bond in list_tree:
        bondedTo[bond[0]].append(bond[1])
        bondedTo[bond[1]].append(bond[0])
    Zs = atoms.get_atomic_numbers()

    if zigzag == 1:

        order = []
        for i in range(0, nx / 2):
            order.append(nx / 2 - i - 1)
            order.append(nx - i - 1)

        # find x and z positions of zigzag carbons
        Cs = edge_carbon_index
        xC = [pos[i][0] for i in edge_carbon_index]
        uxC = np.unique(xC)
        # sort into pairs of atoms at the same x position
        sorted_Cs = [[] for i in range(len(uxC))]
        for C in Cs:
            sorted_Cs[list(uxC).index(pos[C][0])].append(C)
        # bond O to C
        for step,position in enumerate(order):
            for bond in list_tree:
                bondedTo[bond[0]].append(bond[1])
                bondedTo[bond[1]].append(bond[0])
            for k in range(0,cell+1):
                for m in range(0,2):
                    iatom = sorted_Cs[position+k*cell*nx][m]
                    r0 = pos[iatom]
                    bond1 = pos[bondedTo[iatom][0]] - r0
                    bond2 = pos[bondedTo[iatom][1]] - r0
                    CO_bond = -(bond1 + bond2)
                    CO_bond = 1.4 * CO_bond / np.linalg.norm(CO_bond)
                    atoms.append(Atom('O', r0 + CO_bond))
                    angle = 109 * 2 * np.pi / 360
                    rotation = [[np.cos(angle), 0, np.sin(angle)],
                                [0, 1, 0],
                                [-np.sin(angle), 0, np.cos(angle)]]
                    OH_bond = np.dot(-CO_bond, rotation)
                    OH_bond = 0.96 * OH_bond / np.linalg.norm(OH_bond)
                    atoms.append(Atom('H', r0 + CO_bond + OH_bond))

            atoms_temp = atoms.copy()
            daves_super_saturate(atoms_temp)

            atoms1 = atoms_temp.copy()
            atoms1.rotate("x", np.pi/2.0)
            description = "nz" + str(nz) + "cell" + str(cell + 1) + "step" + str(step + 1)
            write(description+".png", atoms1)
            data = make_orca(atoms_temp, filename=description+".inp", output=ORCA_filepath+"/alcohol/"+description+".out", geometry_opt=optimize_geometry)
            moenergies_array = data.moenergies[0]
            all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
            all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos[0]]))
            all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos[0] + 1]))
            descriptions.append(description)

            with open("results.txt", "a+") as e:
                e.write("\n##########################\n")
                e.write(description)
                e.write("\n##########################\n")
                e.write("Total SCF energy in eV:\t")
                e.write(str(data.scfenergies))
                e.write("\nMolecular orbital energy of HOMO in eV:\t")
                e.write(str(moenergies_array[data.homos]))
                e.write("\nMolecular orbital energy of LUMO in eV:\t")
                e.write(str(moenergies_array[data.homos + 1]))


        cleanup(ORCA_filepath + "/alcohol_zigzag/")

        # plot_orbitals(ORCA_filepath + "/alcohol_zigzag/")

    else:
        daves_super_saturate(atoms)

        # find all C atoms
        Cs = []
        for iZ,Z in enumerate(Zs):
            if Z == 6:
                Cs.append(iZ)

        for step in range(0,4):
            # select a random C atom
            C = -1
            while C < 0:
                C = np.random.choice(Cs)
                C = C if len(bondedTo[C]) < 4 else -1

            # add a hydroxyl group to the carbon
            CO_bond = np.array([0,1.4,0]) if np.random.random() > 0.5 else np.array([0,-1.4,0])
            atoms.append(Atom('O', pos[C] + CO_bond))
            angle = 109 * 2 * np.pi / 360
            rotation = [[1, 0, 0],
                        [0, np.cos(angle), -np.sin(angle)],
                        [0, np.sin(angle), np.cos(angle)]]
            OH_bond = np.dot(-CO_bond, rotation)
            OH_bond = 0.96 * OH_bond / np.linalg.norm(OH_bond)
            atoms.append(Atom('H', pos[C] + CO_bond + OH_bond))

            if optimize_geometry == 0:
                xyz = open("step"+str(step+1)+".xyz", "w")
                xyz.write(str(len(atoms))+"\n")
                xyz.write("comment\n")
                symbols_temp = atoms.get_chemical_symbols()
                pos_temp = atoms.get_positions()
                for iatom in range(0,len(atoms)):
                    xyz.write(symbols_temp[iatom]+"\t"+str(pos_temp[iatom][0])+"\t"+str(pos_temp[iatom][1])+"\t"+str(pos_temp[iatom][2])+"\n")
                xyz.close()

        #     atoms1 = atoms.copy()
        #     atoms1.rotate("x", np.pi/2.0)
        #     description = "nz" + str(nz) + "cell" + str(cell + 1) + "step" + str(step + 1)
        #     write(description + ".png", atoms1)
        #     mult = "1" if step%2 == 1 else "2"
        #     data = make_orca(atoms, filename=description + ".inp",
        #                      output=ORCA_filepath + "/alcohol/" + description + ".out", geometry_opt=optimize_geometry,
        #                      multiplicity=mult)
        #     moenergies_array = data.moenergies[0]
        #     all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
        #     all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos[0]]))
        #     all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos[0] + 1]))
        #     descriptions.append(description)
        #
        #     with open("results.txt", "a+") as e:
        #         e.write("\n##########################\n")
        #         e.write(description)
        #         e.write("\n##########################\n")
        #         e.write("Total SCF energy in eV:\t")
        #         e.write(str(data.scfenergies))
        #         e.write("\nMolecular orbital energy of HOMO in eV:\t")
        #         e.write(str(moenergies_array[data.homos]))
        #         e.write("\nMolecular orbital energy of LUMO in eV:\t")
        #         e.write(str(moenergies_array[data.homos + 1]))
        #
        # cleanup(ORCA_filepath + "/alcohol/")

        # plot_orbitals(ORCA_filepath + "/alcohol/")

    return atoms

def form_epoxide(nx=5, nz=3, method="am1", optimize_geometry=0, zigzag=0, atoms=None):
    all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list, descriptions = ([[]] for dummy_var in range(4))
    cell = 0

    atoms = atoms or build_sheet(nx, nz, symmetry=1)
    edge_carbon_index = get_outer_zigzag_atoms(nx, nz)
    daves_super_saturate(atoms)

    # calculate energy of initial sheet
    # atoms1 = atoms.copy()
    # atoms1.rotate("x", np.pi / 2.0)
    # description = "nz" + str(nz) + "cell" + str(cell + 1) + "step0"
    directory = "/epoxide/" if zigzag == 0 else "/epoxide_zigzag/"
    os.popen("mkdir " + ORCA_filepath + directory)
    os.chdir(ORCA_filepath + directory)
    # write(description + ".png", atoms1)
    # data = make_orca(atoms, filename=description + ".inp", multiplicity="1", method=method,
    #                  geometry_opt=optimize_geometry,
    #                  output=ORCA_filepath + directory + description + ".out")
    # moenergies_array = data.moenergies[0]
    # all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
    # all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos[0]]))
    # all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos[0] + 1]))
    # descriptions.append(description)

    # get bond structure data
    pos = atoms.get_positions()
    tree = KDTree(atoms.get_positions())
    list_tree = list(tree.query_pairs(1.430))
    bondedTo = [[] for i in range(len(atoms))]
    for bond in list_tree:
        bondedTo[bond[0]].append(bond[1])
        bondedTo[bond[1]].append(bond[0])
    Zs = atoms.get_atomic_numbers()

    if zigzag == 1:

        order = []
        for i in range(0, nx/2):
            order.append(nx/2-i-1)
            order.append(nx-i-1)

        # find x and z positions of zigzag carbons
        Cs = edge_carbon_index
        xC = [pos[i][0] for i in edge_carbon_index]
        uxC = np.unique(xC)
        # sort into pairs of atoms at the same x position
        sorted_Cs = [[] for i in range(len(uxC))]
        for C in Cs:
            sorted_Cs[list(uxC).index(pos[C][0])].append(C)
        # add epoxide bridges
        for step,position in enumerate(order):
            for k in range(0,cell+1):
                for m in range(0,2):
                    iatom = sorted_Cs[position+k*cell*nx][m]
                    r0 = pos[iatom]

                    if bondedTo[iatom][0] < iatom and bondedTo[iatom][1] < iatom:
                        inner_carbon = np.amax(bondedTo[iatom])
                    else:
                        inner_carbon = np.amin(bondedTo[iatom])
                    CC_bond = pos[inner_carbon] - r0
                    atoms.append(Atom('O', r0 + 0.5*CC_bond + [0,np.sin(np.pi/3.0),0]))

            atoms1 = atoms.copy()
            atoms1.rotate("x", np.pi/2.0)
            description = "nz" + str(nz) + "cell" + str(cell + 1) + "step" + str(step + 1)
            write(description+".png", atoms1)
            data = make_orca(atoms, filename=description+".inp", output=ORCA_filepath+"/epoxide/"+description+".out", geometry_opt=optimize_geometry)
            moenergies_array = data.moenergies[0]
            all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
            all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos[0]]))
            all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos[0] + 1]))
            descriptions.append(description)

            with open("results.txt", "a+") as e:
                e.write("\n##########################\n")
                e.write(description)
                e.write("\n##########################\n")
                e.write("Total SCF energy in eV:\t")
                e.write(str(data.scfenergies))
                e.write("\nMolecular orbital energy of HOMO in eV:\t")
                e.write(str(moenergies_array[data.homos]))
                e.write("\nMolecular orbital energy of LUMO in eV:\t")
                e.write(str(moenergies_array[data.homos + 1]))

            if optimize_geometry == 0:
                xyz = open(description+".xyz", "w")
                xyz.write(str(len(atoms))+"\n")
                xyz.write("comment\n")
                symbols = atoms.get_chemical_symbols()
                pos = atoms.get_positions()
                for iatom in range(0,len(atoms)):
                    xyz.write(symbols[iatom]+"\t"+str(pos[iatom][0])+"\t"+str(pos[iatom][1])+"\t"+str(pos[iatom][2])+"\n")
                xyz.close()

        cleanup(ORCA_filepath+"/epoxide_zigzag/")

        # plot_orbitals(ORCA_filepath+"/epoxide_zigzag/")

    else:
        # find all C atoms
        Cs = []
        for iZ,Z in enumerate(Zs):
            if Z == 6:
                Cs.append(iZ)

        for step in range(0,4):
            # select a random bonded pair
            C1 = -1
            while C1 < 0:
                C1 = np.random.choice(Cs)
                C1 = C1 if len(bondedTo[C1]) < 4 else -1
            C2 = -1
            while C2 < 0:
                C2 = np.random.choice(bondedTo[C1])
                C2 = C2 if Zs[C2] == 6 and len(bondedTo[C2]) < 4 else -1
            Cs.remove(C1)
            Cs.remove(C2)

            # add an epoxide bridge over the pair
            CC_bond = pos[C2] - pos[C1]
            y_disp = [0,np.sin(np.pi/3.0),0] if np.random.random() > 0.5 else [0,-np.sin(np.pi/3.0),0]
            atoms.append(Atom('O', pos[C1] + 0.5*CC_bond + y_disp))

            # atoms1 = atoms.copy()
            # atoms1.rotate("x", np.pi/2.0)
            # description = "nz" + str(nz) + "cell" + str(cell + 1) + "step" + str(step + 1)
            # write(description + ".png", atoms1)
            # data = make_orca(atoms, filename=description + ".inp",
            #                  output=ORCA_filepath + "/epoxide/" + description + ".out", geometry_opt=optimize_geometry)
            # moenergies_array = data.moenergies[0]
            # all_N_SCF_energy_list[cell].append(data.scfenergies[len(data.scfenergies) - 1])
            # all_N_HOMO_energy_list[cell].append(float(moenergies_array[data.homos[0]]))
            # all_N_LUMO_energy_list[cell].append(float(moenergies_array[data.homos[0] + 1]))
            # descriptions.append(description)
            #
            # with open("results.txt", "a+") as e:
            #     e.write("\n##########################\n")
            #     e.write(description)
            #     e.write("\n##########################\n")
            #     e.write("Total SCF energy in eV:\t")
            #     e.write(str(data.scfenergies))
            #     e.write("\nMolecular orbital energy of HOMO in eV:\t")
            #     e.write(str(moenergies_array[data.homos]))
            #     e.write("\nMolecular orbital energy of LUMO in eV:\t")
            #     e.write(str(moenergies_array[data.homos + 1]))

            if optimize_geometry == 0:
                xyz = open("step" + str(step+1) + ".xyz", "w")
                xyz.write(str(len(atoms)) + "\n")
                xyz.write("comment\n")
                symbols = atoms.get_chemical_symbols()
                pos = atoms.get_positions()
                for iatom in range(0, len(atoms)):
                    xyz.write(symbols[iatom] + "\t" + str(pos[iatom][0]) + "\t" + str(pos[iatom][1]) + "\t" + str(
                        pos[iatom][2]) + "\n")
                xyz.close()

        cleanup(ORCA_filepath + "/epoxide/")

        # plot_orbitals(ORCA_filepath + "/epoxide/")

    return atoms

def form_carboxylic_acid(nx=5, nz=3, method="am1", optimize_geometry=0, atoms=None):

    all_N_SCF_energy_list, all_N_HOMO_energy_list, all_N_LUMO_energy_list, descriptions = ([[]] for dummy_var in range(4))
    cell = 0

    atoms = atoms or build_sheet(nx, nz, symmetry=1)
    edge_carbon_index = get_outer_zigzag_atoms(nx,nz)
    # calculate energy of initial sheet
    atoms_temp = atoms.copy()
    daves_super_saturate(atoms_temp)
    description = "nz" + str(nz) + "cell" + str(cell + 1) + "step0"
    directory = "/carboxylic_acid/"
    os.popen("mkdir " + ORCA_filepath + directory)
    os.chdir(ORCA_filepath + directory)
    # get bond structure data
    pos = atoms.get_positions()
    tree = KDTree(atoms.get_positions())
    list_tree = list(tree.query_pairs(1.430))
    bondedTo = [[] for i in range(len(atoms))]
    for bond in list_tree:
        bondedTo[bond[0]].append(bond[1])
        bondedTo[bond[1]].append(bond[0])
    Zs = atoms.get_atomic_numbers()

    for step in range(0,4):
        # select a random C atom on the zigzag edges
        C = np.random.choice(edge_carbon_index)

        # remove adjacent carbons, where there would be no room for an additional group
        ordered_carbons = list(np.unique(edge_carbon_index))
        mod = C%2
        index = ordered_carbons.index(C)
        for i in range(index+1,len(ordered_carbons)):
            if ordered_carbons[i]%2 == mod:
                edge_carbon_index.remove(ordered_carbons[i])
                break
        for i in range(index-1,-1,-1):
            if ordered_carbons[i]%2 == mod:
                edge_carbon_index.remove(ordered_carbons[i])
                break

        edge_carbon_index.remove(C)

        # add a carboxylic acid group
        r0 = pos[C]
        bond1 = pos[bondedTo[C][0]] - r0
        bond2 = pos[bondedTo[C][1]] - r0
        CC_bond = -(bond1 + bond2)
        CC_bond = 1.4 * CC_bond / np.linalg.norm(CC_bond)
        atoms.append(Atom('C', r0 + CC_bond))
        angle1 = 120 * 2.0 * np.pi / 360
        rotation1 = [[np.cos(angle1), 0, np.sin(angle1)],
                     [0, 1, 0],
                     [-np.sin(angle1), 0, np.cos(angle1)]]
        CO_double_bond = np.dot(-CC_bond, rotation1)
        CO_double_bond = 1.2 * CO_double_bond / np.linalg.norm(CO_double_bond)
        atoms.append(Atom('O', r0 + CC_bond + CO_double_bond))
        angle2 = -120 * 2.0 * np.pi / 360
        rotation2 = [[np.cos(angle2), 0, np.sin(angle2)],
                     [0, 1, 0],
                     [-np.sin(angle2), 0, np.cos(angle2)]]
        CO_bond = np.dot(-CC_bond, rotation2)
        CO_bond = 1.4 * CO_bond / np.linalg.norm(CO_bond)
        atoms.append(Atom('O', r0 + CC_bond + CO_bond))
        angle3 = 109 * 2 * np.pi / 360
        rotation3 = [[np.cos(angle3), 0, np.sin(angle3)],
                     [0, 1, 0],
                     [-np.sin(angle3), 0, np.cos(angle3)]]
        OH_bond = np.dot(-CO_bond, rotation3)
        OH_bond = 0.96 * OH_bond / np.linalg.norm(OH_bond)
        atoms.append(Atom('H', r0 + CC_bond + CO_bond + OH_bond))

        if optimize_geometry == 0:
            xyz = open("step"+str(step+1)+".xyz", "w")
            xyz.write(str(len(atoms))+"\n")
            xyz.write("comment\n")
            symbols_temp = atoms.get_chemical_symbols()
            pos_temp = atoms.get_positions()
            for iatom in range(0,len(atoms)):
                xyz.write(symbols_temp[iatom]+"\t"+str(pos_temp[iatom][0])+"\t"+str(pos_temp[iatom][1])+"\t"
                          +str(pos_temp[iatom][2])+"\n")
            xyz.close()

        atoms1 = atoms.copy()
        atoms1.rotate("x", np.pi/2.0)
        write(str(step+1)+".png", atoms1)

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


random_structure(rings=20, pyrroles=0, nitrogens=5, alcohols=0, nCOOH=1)

def cleanup(directory):
    os.chdir(directory)
    files = os.listdir(directory)
    for filename in files:
        if filename[-4:] == ".gbw" or filename[-4:] == ".inp" or filename[-4:] == ".tmp" or filename[-5:] == ".prop"\
                or filename[-7:] == ".engrad" or filename[-4:] == ".opt" or filename[-4:] == ".trj":
            os.remove(directory+filename)

def parse_mulliken(filename):
    file = open(filename)
    loc = 0
    charges = []
    for line in file:
        if loc == 0:
            if re.search("MULLIKEN ATOMIC CHARGES", line):
                loc = 1
        elif loc == 1:
            if re.search("-----", line):
                loc = 2
        elif loc == 2:
            if re.search("Sum of atomic charges", line):
                break
            else:
                words = line.rstrip().split()
                charges.append(float(words[len(words)-1]))
    return charges

def plot_mulliken(data_path):
    os.chdir(data_path)
    os.popen("mkdir mulliken")
    os.chdir("mulliken")
    path_contents = os.listdir(data_path)
    all_charges = []
    all_names = []
    all_x_pos = []
    all_z_pos = []
    for filename in path_contents:
        if filename[-4:] == ".out":
            charges = parse_mulliken(data_path+filename)
            data = parse(data_path+filename)
            fig = plt.figure()
            plt.xlabel("Atom X Position on Sheet")
            plt.ylabel("Atom Z Position on Sheet")
            # cmap = plt.get_cmap("rainbow", np.amax(charges)-np.amin(charges))
            ax = fig.add_subplot(111)
            plt.title(filename[:-4])
            plt.axis("equal")
            pos = data.atomcoords[0]
            x_pos = [pos[i][0] for i in range(0,len(pos))]
            z_pos = [pos[i][2] for i in range(0,len(pos))]
            p = ax.scatter(x_pos, z_pos, s=100, c=charges, marker='o', edgecolors="none", label="something")
            cbar = plt.colorbar(p)
            cbar.set_label("Mulliken Atomic Charge")
            plt.savefig(filename[:-4]+".png")
            plt.clf()
            all_charges.append(charges)
            all_names.append(filename[:-4])
            all_x_pos.append(x_pos)
            all_z_pos.append(z_pos)
    for i in range(0, len(all_charges)-1):
        fig = plt.figure()
        deltas = np.subtract(all_charges[i+1][:-2], all_charges[i])
        plt.xlabel("Atom X Position on Sheet")
        plt.ylabel("Atom Y Position on Sheet")
        ax = fig.add_subplot(111)
        plt.title(all_names[i] + " -> " + all_names[i+1])
        plt.axis("equal")
        p = ax.scatter(all_x_pos[i], all_z_pos[i], s=100, c=deltas, marker='o', edgecolors="none", label="something")
        cbar = plt.colorbar(p)
        cbar.set_label("Change in Mulliken Atomic Charge")
        plt.savefig("step"+str(i+1)+".png")
        plt.clf()

def condense_output(directory):
    os.chdir(directory)
    path_contents = os.listdir(directory)

    results = []

    for filename in path_contents:
        if filename[-4:] == ".out":
            file = ccopen(directory + filename)
            data = file.parse()
            try:
                natom = data.natom
                scf = data.scfenergies[len(data.scfenergies) - 1]
                pos = data.atomcoords[len(data.atomcoords) - 1]
                tree = KDTree(pos)
                list_tree = list(tree.query_pairs(1.430))
                bondedTo = [[] for i in range(data.natom)]
                for bond in list_tree:
                    bondedTo[bond[0]].append(bond[1])
                    bondedTo[bond[1]].append(bond[0])
                Zs = data.atomnos
                CH_index = []
                NH_index = []
                H_attachment = []
                for iatom in range(0, natom):
                    if (Zs[iatom] == 6 and len(bondedTo[iatom]) == 4) or (
                            Zs[iatom] == 7 and len(bondedTo[iatom]) == 3):
                        H_attachment.append(1)
                    elif Zs[iatom] == 6 or Zs[iatom] == 7:
                        H_attachment.append(0)
                moenergies_array = data.moenergies[0]
                homo = float(moenergies_array[data.homos[0]])
                lumo = float(moenergies_array[data.homos[0] + 1])
                print(filename)

                results.append({
                    "name": filename[:-4],
                    "natom": natom,
                    "scf": scf,
                    "H_attachment": H_attachment,
                    "Zs": Zs,
                    "homo": homo,
                    "lumo": lumo
                })

            except AttributeError:
                print(filename + "\tAttributeError")

        os.remove(filename)

    os.rmdir(directory)
    pickle.dump(results, open(directory[:-1], "w"))

def pickle_complete_output(directory):
    os.chdir(directory)
    path_contents = os.listdir(directory)

    results = []

    for filename in path_contents:
        if filename[-4:] == ".out":
            file = ccopen(directory + filename)
            data = file.parse()
            results.append(data)
            os.remove(filename)

    pickle.dump(results, open("out", "w"))

def draw_sheet_map(nx=5, nz=3):
    atoms = build_sheet(nx, nz, symmetry=1)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for i, atom in enumerate(atoms):
        plt.plot(atom.position[0], atom.position[2], '.')
        plt.text(atom.position[0], atom.position[2], str(i))
    plt.savefig(str(nx)+"x"+str(nz)+"map.png")

###Tkinter Section Below
################################################################################################################
################################################################################################################
#Requires 'hovering cursor' tooltip for all major components

# root = Tk()
# root.title("Graphene Nitrogenator")
# Label(root, text="GUI for Nitrogenated Graphene Energy Calculations").pack()


class sheet_param:
    def __init__(self, master):
        global horizon_sheet_variable
        horizon_sheet_variable = StringVar()
        global vertical_sheet_variable
        vertical_sheet_variable = StringVar()
        param_frame = Frame(master)
        Label(param_frame, text="horizontal").grid(row=1, column=1)
        Label(param_frame, text="vertical").grid(row=1, column=3)
        Label(param_frame, text="Sheet Dimensions").grid(row=2, sticky=E)
        Label(param_frame, text="x").grid(row=2, column=2)
        Entry(param_frame, width=2, textvariable=horizon_sheet_variable).grid(row=2, column=1)
        Entry(param_frame, width=2, textvariable=vertical_sheet_variable).grid(row=2, column=3)
        param_frame.pack()


class checkbox:
    def __init__(self,master):
        self.cbox_frame = Frame(master)
        self.cbox_frame.pack()
        
    def symtrc_cbox(self):
        global symmetry_var
        symmetry_var = IntVar()
        Checkbutton(self.cbox_frame, text="Make molecule symmetric", variable=symmetry_var).pack(anchor=W)

    def unsaturate_cbox(self):
        global unsat_var
        unsat_var = IntVar()
        Checkbutton(self.cbox_frame, text="Leave molecule unsaturated", variable=unsat_var, onvalue=1, offvalue=0).pack(anchor=W)
        #unsat_var_string = str(unsat_var.get())

    def opt_geo_cbox(self):
        opt_geo_var = IntVar()
        Checkbutton(self.cbox_frame, text="Optimize Geometry", variable=opt_geo_var, onvalue=1, offvalue=0).pack(anchor=W)
        global selected_geometry_opt
        selected_geometry_opt = int(opt_geo_var.get())

    def nitrogenate_all_cbox(self):
        global N_all_cbox_var
        N_all_cbox_var = IntVar()
        Checkbutton(self.cbox_frame, text="Nitrogenate all zig-zag positions", variable=N_all_cbox_var, onvalue=1, offvalue=0).pack(anchor=W)
#can't seem to get checkboxes to align


class drop_list:
    def __init__(self, master):
        self.drop_list_frame = Frame(master)
        self.drop_list_frame.pack()
        self.method_var = StringVar()
        self.calculator_var = StringVar()

    def calc_combobox(self):
        self.calculator_var.set("ORCA")
        OptionMenu(self.drop_list_frame, self.calculator_var, "ORCA", "NWChem", "Gaussian").grid(row=0, column=1)
        Label(self.drop_list_frame, text="Calculator").grid(row=0, sticky=E)
        global selected_calculator
        selected_calculator = str(self.calculator_var.get())

    def calc_method(self):
        self.method_var.set("am1")
        OptionMenu(self.drop_list_frame, self.method_var, "am1", "DFT").grid(row=1, column=1)
        Label(self.drop_list_frame, text="Method").grid(row=1, sticky=E)
        global selected_calc_method
        selected_calc_method = str(self.method_var.get())



class button:
    def __init__(self, master):
        self.button_frame = Frame(master)
        self.button_frame.pack(expand="yes")

    def bottom_buttons(self):
        Button(self.button_frame, text="Calculate", command=gui_calculate).grid(row=0, column=3)
        Button(self.button_frame, text="View in Avogadro", command=gui_view).grid(row=0, column=0)

    def test_button(self):
        Button(self.button_frame, text="Test Button", command=gui_nitrogenate_all).grid(row=1)

    def print_var(self):
        string = str(horizon_sheet_variable.get())
        print(string)


def build_param_frame(master):
    lbl_bld_frame = LabelFrame(master, text="Graphene Builder Parameters", padx=5, pady =5)
    lbl_bld_frame.pack(expand="yes")
    sheet_param(lbl_bld_frame)
    checkbox(lbl_bld_frame).symtrc_cbox()
    checkbox(lbl_bld_frame).unsaturate_cbox()
    checkbox(lbl_bld_frame).nitrogenate_all_cbox()

def calc_param_frame(master):
    lbl_calc_frame = LabelFrame(master, text="Calculation Parameters", padx=5, pady=5)
    lbl_calc_frame.pack(expand="yes")
    #drop_list(lbl_calc_frame).calc_combobox()
    drop_list(lbl_calc_frame).calc_method()
    checkbox(lbl_calc_frame).opt_geo_cbox()

def gui_view():
        horizontal_dimension = int(horizon_sheet_variable.get())
        vertical_dimension = int(vertical_sheet_variable.get())
        symmetry_int = int(symmetry_var.get())
        atoms = build_sheet(horizontal_dimension, vertical_dimension, symmetry=symmetry_int)
        unsat_int = int(unsat_var.get())
        if unsat_int==0:
            daves_super_saturate(atoms)
        elif unsat_int==1:
            pass
        view(atoms, viewer="avogadro")

def gui_calculate():
    horizontal_dimension = int(horizon_sheet_variable.get())
    vertical_dimension = int(vertical_sheet_variable.get())
    symmetry_int = int(symmetry_var.get())
    global atoms
    atoms = build_sheet(horizontal_dimension, vertical_dimension, symmetry=symmetry_int)
    # global unsat_int
    # unsat_int = int(unsat_var.get())

    if N_all_cbox_var==0:
        calc_edge_nitrogens(horizontal_dimension, vertical_dimension, method=selected_calc_method, optimize_geometry=selected_geometry_opt, make_symmetric=symmetry_int)
    elif N_all_cbox_var==1:
        nitrogenate_all_zig_zag(horizontal_dimension, vertical_dimension, method=selected_calc_method, optimize_geometry=selected_geometry_opt, make_symmetric=symmetry_int)

def gui_nitrogenate_all():
    horizontal_dimension = int(horizon_sheet_variable.get())
    vertical_dimension = int(vertical_sheet_variable.get())
    symmetry_int = int(symmetry_var.get())
    global atoms
    atoms = build_sheet(horizontal_dimension, vertical_dimension, symmetry=symmetry_int)
    global unsat_int
    unsat_int = int(unsat_var.get())
    nitrogenate_all_zig_zag(horizontal_dimension, vertical_dimension, method=selected_calc_method, optimize_geometry=selected_geometry_opt, make_symmetric=symmetry_int)


#saturated_nitrogenate_all_zig_zag(nx_min=6, nx_max=6, nz_min=6, nz_max=6, method="am1", optimize_geometry=0, make_symmetric=1)
        #build_param_frame(root)
        #calc_param_frame(root)
        #button(root).bottom_buttons()
                #button(root).test_button()
        #root.mainloop()

#atoms = build_sheet(3, 3, symmetry=1)
#nitrogenate(atoms, 34)
#daves_super_saturate(atoms)
#view(atoms, viewer="avogadro")
#calc_edge_nitrogens(3, 3, method="am1", optimize_geometry=False, make_symmetric=1)
#print data.atomcharges
#############################################Concurrent Bug List##########################################################
#asymmetric molecules are viewed as symmetric in .png energy map files for 3x3 sheet
