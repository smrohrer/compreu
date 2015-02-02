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
from Tkinter import *

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
        print "Making graphene sheet symmetric"
        global to_be_removed
        to_be_removed = []
        sym_fix_int = nz-1
        for i in xrange(2, nz+sym_fix_int, 2):
            to_be_removed.append(2*i-2)
            to_be_removed.append(2*i-1)
        
        to_be_removed = to_be_removed[::-1]

        for entry in to_be_removed:
            atoms.pop(entry)
    else:
        print "Molecule was not made symmetric"
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

def make_orca(atoms, filename="filename.inp", charge="0", multiplicity="1", method="am1", geometry_opt=0, output="temp.out"):
    if method == "am1" and geometry_opt == 0:
        out=''
        parameters0= '{0}\t{1}\t{2}\n{3}'.format("%method", "method", method, "end")
        parameters1='{0}\t{1}\t{2}\t{3}\n'.format("*", "xyz", charge, multiplicity)
        out=out+parameters0+parameters1
        end_of_atom_coordinates="*"

    elif method == "am1" and geometry_opt == 1:
        out=''
        parameters0= '{0}\t{1}\t{2}\t{3}\t{4}\n{5}'.format("%method", "method", method, "method", "OPT", "end")
        parameters1='{0}\t{1}\t{2}\t{3}\n'.format("*", "xyz", charge, multiplicity)
        out=out+parameters0+parameters1
        end_of_atom_coordinates="*"

    elif method == "DFT" and geometry_opt== 0:
        out=''
        parameters0= '{0}\t{1}\n'.format("!", "DFT-Energy")
        parameters1='{0}\t{1}\t{2}\t{3}\n'.format("*", "xyz", charge, multiplicity)
        out=out+parameters0+parameters1
        end_of_atom_coordinates="*"

    elif method == "DFT" and geometry_opt == 1:
        out=''
        parameters0= '{0}\t{1}\n'.format("!", "Good-Opt")
        parameters1='{0}\t{1}\t{2}\t{3}\n'.format("*", "xyz", charge, multiplicity)
        out=out+parameters0+parameters1
        end_of_atom_coordinates="*"

    with open(filename, 'w') as f:
        f.write(out)
        f.write(print_atoms(atoms))
        f.write(end_of_atom_coordinates)
    subprocess.call(ORCA_filepath + "/orca/orca "+ filename + " > " + output, shell=True)
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


def calc_edge_nitrogens(nx="1", nz="1", method="am1", optimize_geometry=0, make_symmetric=0, sub_all_zigzag=0):
    if make_symmetric == 1:
        print "atoms sheet at start of calc is symmetric"
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

    for number in xrange(0, 2*nx):
        addition.append(number)
    for number in xrange(2, 2*(nx+2), 2):
        multiplication.append(number)
        multiplication.append(number)
    for value in xrange(0, len(addition)):
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
            print "nitrogenate make atoms sheet is symmetric"
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
    for i in xrange(3):
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

def nitrogenate_all_zig_zag(nx, nz, method="am1", optimize_geometry=0, make_symmetric=0):
    if make_symmetric == 1:
        atoms = build_sheet(nx, nz, symmetry=1)
        symmetry_folder_string = "_symmetric"
    else:
        atoms = build_sheet(nx, nz, symmetry=0)
        symmetry_folder_string = "_asymmetric"

    addition = [0]
    multiplication = []
    edge_carbon_index =[]

    for number in xrange(0, 2*nx):
        addition.append(number)
    for number in xrange(2, 2*(nx+2), 2):
        multiplication.append(number)
        multiplication.append(number)
    for value in xrange(0, len(addition)):
        edge_carbon_index.append(nz*multiplication[value]+addition[value])
    edge_carbon_index.pop(1)

    if make_symmetric == 1:
        edge_carbon_index[:] = [x-len(to_be_removed) for x in edge_carbon_index]
    else:
        pass
    for edge_carbon in edge_carbon_index:
        #nitrogenate(atoms, edge_carbon)
        symbols = atoms.get_chemical_symbols()
        symbols[edge_carbon] = 'N'
        atoms.set_chemical_symbols(symbols)
    view(atoms, viewer="avogadro")




###Tkinter Section Below
################################################################################################################
################################################################################################################
#Requires 'hovering cursor' tooltip for all major components

root = Tk()
root.title("Graphene Nitrogenator")
Label(root, text="GUI for Nitrogenated Graphene Energy Calculations").pack()


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
        Checkbutton(self.cbox_frame, text="Optimize Geometry", variable=opt_geo_var).pack(anchor=W)
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
        Button(self.button_frame, text="Nitrogenate All", command=gui_nitrogenate_all).grid(row=1)

    def print_var(self):
        string = str(horizon_sheet_variable.get())
        print string


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
    global unsat_int
    unsat_int = int(unsat_var.get())

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
    



build_param_frame(root)
calc_param_frame(root)
button(root).bottom_buttons()

button(root).test_button()
root.mainloop()

#atoms = build_sheet(3, 3, symmetry=1)
#nitrogenate(atoms, 34)
#daves_super_saturate(atoms)
#view(atoms, viewer="avogadro")
#calc_edge_nitrogens(3, 3, method="am1", optimize_geometry=False, make_symmetric=1)
#print data.atomcharges
#############################################Concurrent Bug List##########################################################
#asymmetric molecules are viewed as symmetric in .png energy map files for 3x3 sheet