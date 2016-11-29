from orca import Orca
from structure_building import random_structure

struct = random_structure(rings=10)
orc=Orca()
orc.set_atoms(struct)
orc.write_inp(struct)
