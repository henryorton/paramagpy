import cProfile

import numpy as np
from bonded_to import Tracker
from pcs_reader import PCSReader

from paramagpy import protein, metal, fit


def rad(deg): return (deg / 180) * np.pi


def fit_metal():
    # Pre-computed list containing the actual names of amide-H in each res (Use <CustomAtom>.bonded_to() to get this)
    amide_h_list = [[None], ['H07'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'],
                    ['H02'],
                    ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], [], ['H02'], ['H02'], ['H02'],
                    ['H02'],
                    ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'],
                    ['H02'],
                    ['H02'], [], [], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'],
                    ['H02'],
                    ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'],
                    ['H02'],
                    ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'],
                    ['H02'],
                    ['H02'], ['H02'], ['H02'], ['H02']]

    # Load the PDB file
    prot = protein.load_pdb('../data_files/1ubq.pdb')

    # Rename amide-H's to 'H'
    for chain in prot[0]:
        for res in chain:
            if res.id[1] > 75:
                continue
            new_h = amide_h_list[res.id[1]]
            if len(new_h) > 0:
                new_h = res[new_h[0]]
                new_h.id = new_h.fullname = new_h.name = 'H'
                res.child_list.append(new_h)
                res.child_dict['H'] = new_h

    # Reading PCS Values
    pcs = PCSReader('../data_files/s57c_pcsexp_pcscalc.list', lambda x: 'H' in x['Atom'])
    pcs_data_exp = pcs.get_result()

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    # Set the starting position to an atom close to the metal
    m_start = metal.Metal(position=prot[0]['A'][57]['CB'].position)

    # Calculate an initial tensor from an SVD grid search
    m_guess, calc, q_fac = fit.svd_gridsearch_fit_metal_from_pcs([m_start], [parsed_data], radius=10, points=10)

    # Refine the tensor using non-linear regression
    m_fit, calc, q_fac = fit.nlr_fit_metal_from_pcs(m_guess, [parsed_data])
    m_fit[0].save('s57c_1ubq_PCS_tensor.txt')
    print(m_guess, calc, q_fac)


def fit_rotamer():
    # Load the PDB file
    prot = protein.load_pdb('../data_files/1ubq.pdb')

    # Reading PCS Values
    pcs = PCSReader('../data_files/s57c_pcsexp_pcscalc.list')
    pcs_data_exp = pcs.get_result()

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    # Define an initial tensor
    chi_tensor = {"axrh": (-9.178E-32, -1.452E-32),
                  "position": (21.165E-10, 11.566E-10, 10.674E-10),
                  "eulers": (0.95855228, 1.05032169, 2.02433759)}
    mtl = metal.Metal(**chi_tensor)

    pos_array = np.empty((len(parsed_data), 3))
    for idx, data in enumerate(parsed_data):
        pos_array[idx] = data[0].position

    pcs_data_calc = mtl.fast_pcs(pos_array)
    for idx, data in enumerate(parsed_data):
        print(data[0], data[1], pcs_data_calc[idx])

    pcs_fitter = fit.PCSToRotamer(prot[0], mtl, parsed_data)

    trk = Tracker(pcs_fitter.model['A'][6])
    trk.print_atom_linkage()
    trk.save_atom_coords()

    pr = cProfile.Profile()

    def staggered_positions_search():
        pr.clear()
        pr.enable()
        # Grid search 5 degrees about each staggered position
        min_pcs = pcs_fitter.run_staggered_positions_search('A', 6, 0.174533, 5)
        pr.disable()
        pr.print_stats(sort="cumtime")
        print(min_pcs)

    def grid_search():
        pr.clear()
        pr.enable()
        rot_param = np.array([[0, -2 / 9 * np.pi, 9]] * 4)
        pcs_fitter.set_rotation_parameter('A', 6, rot_param)
        result1 = pcs_fitter.run_grid_search()
        pr.disable()
        pr.print_stats(sort="cumtime")
        print(result1)

    def pairwise_grid_search():
        pr.clear()
        pr.enable()
        rot_param = np.array([[0, -2 / 6 * np.pi, 6]] * 4)
        pcs_fitter.set_rotation_parameter('A', 6, rot_param)
        result2 = pcs_fitter.run_pairwise_grid_search()
        pr.disable()
        pr.print_stats(sort="cumtime")
        print(result2)

    staggered_positions_search()
    grid_search()
    pairwise_grid_search()

    trk.save_atom_coords()
    trk.print_atom_coords()


if __name__ == '__main__':
    # fit_metal()
    fit_rotamer()
