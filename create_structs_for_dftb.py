from ase.calculators.dftb import Dftb
from structure_building import random_structure


struct = random_structure(rings=10, pyrroles=0, nitrogens=2, alcohols=0, nCOOH=0)

calc = Dftb(label='')

struct.set_calculator(calc)
calc.calculate(system)
