import cProfile
import time

import numpy as np
from Bio.PDB import PDBIO
from tabulate import tabulate

from paramagpy import protein


class Tracker:

    def __init__(self, rs):
        self.rs = rs
        self.res_atoms = {}

    def save_atom_coords(self):
        for x in self.rs.get_atoms():
            if self.res_atoms.get(x) is None:
                self.res_atoms[x] = []
            self.res_atoms[x] += [x.get_coord()]

    def print_atom_coords(self):
        for atom in self.res_atoms:
            print(f"{atom} : {self.res_atoms[atom]}")

    def print_atom_linkage(self):
        atoms = list(self.rs.get_atoms())
        for atom in atoms:
            print(f"{atom} : {[x for x in atom.bonded_to()]}")


def rad(deg): return (deg / 180) * np.pi


def save_structure(structure, filename):
    # Save structure in PDB format

    io = PDBIO()
    io.set_structure(structure)
    io.save(filename)


def main():
    timer = {'start_program': time.perf_counter() * 1E3, 'start_pdb_parser': time.perf_counter() * 1E3}
    prot = protein.load_pdb('../data_files/lys.pdb')
    timer['end_pdb_parser'] = time.perf_counter() * 1E3

    res = prot[0]['A'][55]
    trk = Tracker(res)

    timer['start_atom_linkage'] = time.perf_counter() * 1E3
    trk.print_atom_linkage()
    timer['end_atom_linkage'] = time.perf_counter() * 1E3

    trk.save_atom_coords()

    # Default: [rad(47.63), rad(162.49), rad(-176.18), rad(141.70)
    timer['start_set_dihedral'] = time.perf_counter() * 1E3
    res.set_dihedral(np.array([rad(170), rad(170), rad(-70), rad(50)]))
    # res.set_dihedral(np.array([-2.0943951023931957, 1.3962634015954627, -2.0943951023932383, 0]))
    timer['end_set_dihedral'] = time.perf_counter() * 1E3

    # Will give back original if the below line is executed
    timer['start_set_delta_dihedral'] = time.perf_counter() * 1E3
    # res.set_delta_dihedral(np.array([rad(10), 0, 0, rad(-5)]))
    timer['end_set_delta_dihedral'] = time.perf_counter() * 1E3

    # Full sweep

    pr = cProfile.Profile()

    pr.enable()
    timer['start_set_delta_dihedral_full'] = time.perf_counter() * 1E3
    rot_param = np.array([[-7 / 9 * np.pi, np.pi, 9] * 4]).reshape(4, 3)
    res.grid_search_rotamer(rot_param)
    timer['end_set_delta_dihedral_full'] = time.perf_counter() * 1E3
    pr.disable()
    pr.print_stats(sort="cumtime")

    trk.save_atom_coords()
    trk.print_atom_coords()

    timer['end_program'] = time.perf_counter() * 1E3

    # Print benchmark results
    print()
    print()
    print('=' * 50)
    print('Benchmark Results')
    print('=' * 50)
    print(tabulate([['Program', timer['end_program'] - timer['start_program']],
                    ['PDB Parsing', timer['end_pdb_parser'] - timer['start_pdb_parser']],
                    ['Atom Linkages', timer['end_atom_linkage'] - timer['start_atom_linkage']],
                    ['Set Dihedral', timer['end_set_dihedral'] - timer['start_set_dihedral']],
                    ['Set delta-dihedral', timer['end_set_delta_dihedral'] - timer['start_set_delta_dihedral']],
                    ['Set delta-dihedral Full',
                     timer['end_set_delta_dihedral_full'] - timer['start_set_delta_dihedral_full']]],
                   headers=['Name', 'Duration (ms)'], tablefmt="fancy_grid"))

    save_structure(prot, "lys_conformer.pdb")


def _test():
    prot = protein.load_pdb('../data_files/1ubq_k6.pdb')
    res = prot[0]['A'][6]

    res['CG'].position = np.array([2.43427963e-09, 3.93229787e-09, 1.43778842e-09])
    res['CD'].position = np.array([2.37847472e-09, 3.91072686e-09, 1.29911211e-09])
    res['CE'].position = np.array([2.58744302e-09, 3.91582533e-09, 1.44857787e-09])

    res['CG'].position = np.array([2.53592950e-09, 3.82112766e-09, 1.20278980e-09])
    res['CD'].position = np.array([2.61213267e-09, 3.88580137e-09, 1.08955993e-09])
    res['CE'].position = np.array([2.56654200e-09, 3.87964746e-09, 1.34237880e-09])

    save_structure(prot, "1ubq_k6_conformer.pdb")


if __name__ == '__main__':
    main()
    _test()
