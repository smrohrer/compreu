from ase.calculators.dftb import Dftb
import os

def dftb_calc(path, calc_folder, sys):
    os.chdir(path)
    if not os.path.exists(calc_folder):
        os.makedirs(calc_folder)
    os.chdir(calc_folder)
    # Assuming C and H are always present
    args = {}
    if 'N' in sys.get_chemical_symbols():
        args['Hamiltonian_MaxAngularMomentum_N'] = '"p"'
    if 'O' in sys.get_chemical_symbols():
        args['Hamiltonian_MaxAngularMomentum_O'] = '"p"'

    calc = Dftb(label=calc_folder,
                atoms=sys,
                run_manyDftb_steps=True,
                Driver_='ConjugateGradient',
                Driver_MaxForceComponent='1E-4',
                Driver_MaxSteps=1000,
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='"p"',
                Hamiltonian_MaxAngularMomentum_H='"s"',
                **args)

    sys.set_calculator(calc)
    calc.calculate(sys)
    os.chdir(os.pardir)
    return calc
