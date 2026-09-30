# Drop-in replacement for calculate_eda.py using gpu4pyscf ALMO-EDA2 instead of Q-Chem.
# Output JSON has the same keys and units (kcal/mol) as calculate_eda.py.

import glob
import json

import pyscf
from gpu4pyscf.properties import eda

from bytemol.core import Molecule

KJ_TO_KCAL = 1 / 4.184


def run_eda(atoms, net_charges):
    mols = [
        pyscf.M(atom=[(s, tuple(p)) for s, p in zip(a.symbols, a.positions)],
                basis='def2-TZVPD',
                charge=c,
                unit='Angstrom') for a, c in zip(atoms, net_charges)
    ]
    res, _ = eda.eval_ALMO_EDA_2_energies(mols, xc='wB97M-V')
    assert res['unit'] == 'kJ/mol'
    values = {
        'ELEC': res['electrostatic'],
        'CLS ELEC': res['classical electrostatic'],
        'ELEC_PAULI': res['electrostatic'] + res['pauli'],
        'DISP': res['dispersion'],
        'FROZEN': res['frozen'],
        'POLARIZATION': res['polarization'],
        'CHARGE TRANSFER': res['charge transfer'],
        'PREPARATION': 0.0,
        'TOTAL': res['total'],
    }
    assert abs(values['ELEC_PAULI'] + values['DISP'] - values['FROZEN']) < 1e-3
    assert abs(values['FROZEN'] + values['POLARIZATION'] + values['CHARGE TRANSFER'] - values['TOTAL']) < 1e-3
    return {k: v * KJ_TO_KCAL for k, v in values.items()}


def main(input_dir):
    mols = [Molecule.from_xyz(xyz) for xyz in glob.glob(f'{input_dir}/*.xyz')]
    assert len(mols) == 2
    atoms, net_charges, names = [], [], []
    for mol in mols:
        atoms.append(mol.conformers[0].to_ase_atoms())
        net_charges.append(int(sum(mol.formal_charges)))
        names.append(mol.name)
    eda_results = run_eda(atoms, net_charges)
    with open(f'{names[0]}_{names[1]}.json', 'w') as f:
        json.dump(eda_results, f, indent=4)  # kcal/mol


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, default='ACT_EC_dimer/conf_0')
    args = parser.parse_args()
    main(args.input_dir)
