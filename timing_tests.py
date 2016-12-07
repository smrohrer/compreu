from __future__ import print_function
from dftb_wrapper import dftb_calc
from structure_building import random_structure, name_structure
import timeit
import random

# ref: struct = random_structure(rings=10, pyrroles=0, nitrogens=2, alcohols=0, COOH=0)


for i in range(20):
    kwords = {'rings' : random.randint(8,12),
            'pyrroles' : random.randint(0,2),
            'nitrogens' : random.randint(0,2),
            'alcohols' : random.randint(0,2),
            'COOH' : random.randint(0,1),
            'epoxide' : 0}

    test = random_structure(**kwords)
    name = name_structure(test,**kwords)

    path = '/home/egottlie/dftb_calcs/'
    start_time = timeit.default_timer()
    calc = dftb_calc(path,name,test)
    elapsed = timeit.default_timer() - start_time
    
    f = open(name+',time','w')
    print(elapsed,file=f)
    f.close()

