from dftb_wrapper import dftb_calc
from structure_building import random_structure, name_structure


# ref: struct = random_structure(rings=10, pyrroles=0, nitrogens=2, alcohols=0, COOH=0)

kwords = {'rings' : 8,
          'pyrroles' : 0,
          'nitrogens' : 2,
          'alcohols' : 1,
          'COOH' : 1,
          'epoxide' : 0}

test = random_structure(**kwords)
name = name_structure(test,**kwords)

path = '/home/egottlie/dftb_calcs/'
calc = dftb_calc(path,name,test)
