from dftb_wrapper import dftb_calc
from structure_building import random_structure
import os

test = random_structure(rings=5,nitrogens=0,epoxide=1)
path = '/home/egottlie/dftb_calcs/'
os.chdir(path)
calc = dftb_calc('/home/egottlie/dftb_calcs/test_calc',test)
test.set_calculator(calc)
calc.calculate(test)
