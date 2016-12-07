from __future__ import print_function
from dftb_wrapper import dftb_calc
from structure_building import random_structure, name_structure
import timeit
import random

# ref: struct = random_structure(rings=10, pyrroles=0, nitrogens=2, alcohols=0, COOH=0)


for i in range(100):
    nRings = random.randint(5,20)
    kwords = {'rings' : nRings,
            'pyrroles' : random.randint(0,nRings/5),
            'nitrogens' : random.randint(0,nRings/2),
            'alcohols' : random.randint(0,nRings/5),
            'COOH' : random.randint(0,nRings/5),
            'epoxide' : 0}

    test = random_structure(**kwords)
    name = name_structure(test,**kwords)

    path = '/home/egottlie/dftb_calcs/'
    calc = dftb_calc(path,name,test)
    
