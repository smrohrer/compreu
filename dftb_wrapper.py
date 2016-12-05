from ase.calculators.dftb import Dftb

def dftb_calc(input_label,sys):
    # Assuming C and H are always present
    args = {}
    if 'N' in sys.get_chemical_symbols():
        args['Hamiltonian_MaxAngularMomentum_N'] = '"p"'
    if 'O' in sys.get_chemical_symbols():
        args['Hamiltonian_MaxAngularMomentum_O'] = '"p"'
    calc = Dftb(label=input_label,
                atoms=sys,
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='"p"',
                Hamiltonian_MaxAngularMomentum_H='"s"',
                **args)
    return calc
