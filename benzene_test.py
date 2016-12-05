from dftb_wrapper import dftb_calc
from structure_building import random_structure

test = random_structure(nitrogens=1)
calc = dftb_calc('test_calc',test)
test.set_calculator(calc)
calc.calculate(test)
