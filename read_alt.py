#get eigenvalues in ev, fillings (for homo/lumo), total energy, orbital energy

"""https://gitlab.com/ase/ase/blob/master/ase/calculators/dftb.py
http://www.dftb-plus.info/fileadmin/DFTB-Plus/public/recipes/html/basics/firstcalc.html"""


#import os
#import numpy as np
#from ase.calculators.dftb import Dftb
#from ase.calculators.calculator import FileIOCalculator, kpts2mp

def read_alt():
    """these results are pulled from the "detailed.out" file. 
    So far I've pulled eigenvalues and fillings, which can be used for homo/lumo 
    Note that these are spin unpolarized calculations"""

    myfile2 = open("detailed.out", "r")
    lines2 = myfile2.readlines() #changed to lines2
    myfile2.close()
    
    eigenvalues = 0
    fillings = 0
    eigenvalues = read_eigenvalues(lines2)
    fillings = read_fillings(lines2)

    return [eigenvalues,fillings]#why won't it return RIP

def read_eigenvalues(lines2):
    
    # eigenvalues line index
    for iline, line in enumerate(lines2):
        evstring = 'Eigenvalues /eV'
        if line.find(evstring) >= 0:
            index_eigenvalues_start = iline + 1
            break
    try:
        eigenvals = []
        i = index_eigenvalues_start
        try:
            while lines2[i].split()[0] != "\n":
                eigenvals.append(lines2[i].split()[0])
                i += 1
        except: 
            return eigenvals

    except:
        raise RuntimeError('Problem in reading eigenvalues')

def read_fillings(lines2):
    # eigenvalues line index
    for iline, line in enumerate(lines2):
        fstring = 'Fillings'
        if line.find(fstring) >= 0:
            index_fillings_start = iline + 1
            break
    try:
        fillings = []
        i = index_fillings_start
        try:
            while lines2[i].split()[0] != "\n":
                fillings.append(lines2[i].split()[0])
                i += 1
        except: 
            return fillings
    except:
        raise RuntimeError('Problem in reading fillings')
    

    
  
