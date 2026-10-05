# Drop-in replacement for generate_conf.py using gpu4pyscf instead of Q-Chem.
# Same level of theory (B3LYP-D3BJ/def2-SVPD, geomeTRIC, GAU convergence) and same output xyz.

import ase
import ase.io
import pyscf
from gpu4pyscf.dft import rks
from pyscf.geomopt.geometric_solver import optimize

from bytemol.core import Molecule


def main(mapped_smiles, mol_name):
    mol = Molecule.from_mapped_smiles(mapped_smiles, nconfs=1, name=mol_name)
    net_charge = int(sum(mol.formal_charges))
    atoms = mol.conformers[0].to_ase_atoms()
    pmol = pyscf.M(atom=[(s, tuple(p)) for s, p in zip(atoms.symbols, atoms.positions)],
                   basis='def2-SVPD',
                   charge=net_charge,
                   spin=0,
                   unit='Angstrom',
                   cart=False)
    if pmol.natm > 1:
        mf = rks.RKS(pmol, xc='b3lyp')
        mf.disp = 'd3bj'
        mf.conv_tol = 1e-10
        mf.max_cycle = 200
        pmol = optimize(mf, conv_params={'convergence_set': 'GAU'})
    atoms = ase.Atoms(symbols=[pmol.atom_symbol(i) for i in range(pmol.natm)],
                      positions=pmol.atom_coords(unit='Angstrom'))
    atoms.info['mapped_smiles'] = mapped_smiles
    atoms.info['name'] = mol_name
    ase.io.write(f"./{mol_name}.xyz", atoms)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--mol_name", type=str, default="ACT")
    parser.add_argument("--mapped_smiles",
                        type=str,
                        default="[C:1]([C:2](=[O:3])[C:4]([H:8])([H:9])[H:10])([H:5])([H:6])[H:7]")
    args = parser.parse_args()
    main(args.mapped_smiles, args.mol_name)
