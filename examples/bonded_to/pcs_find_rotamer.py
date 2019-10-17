import cProfile

import numpy as np
from Bio.PDB import PDBIO
from bonded_to import Tracker

from paramagpy import protein, metal, fit, dataparse


def rad(deg): return (deg / 180) * np.pi


def save_structure(structure, filename):
    # Save structure in PDB format

    io = PDBIO()
    io.set_structure(structure)
    io.save(filename)


def fit_metal():
    # Pre-computed list containing the actual names of amide-H in each res (Use <CustomAtom>.bonded_to() to get this)
    # amide_h_list = [[None], ['H07'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'],
    #                 ['H02'],
    #                 ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], [], ['H02'], ['H02'], ['H02'],
    #                 ['H02'],
    #                 ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'],
    #                 ['H02'],
    #                 ['H02'], [], [], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'],
    #                 ['H02'],
    #                 ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'],
    #                 ['H02'],
    #                 ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'], ['H02'],
    #                 ['H02'],
    #                 ['H02'], ['H02'], ['H02'], ['H02']]

    # Load the PDB file
    prot = protein.load_pdb('../data_files/1UBQ_H.pdb')

    # Rename amide-H's to 'H'
    # for chain in prot[0]:
    #     for res in chain:
    #         if res.id[1] > 75:
    #             continue
    #         new_h = amide_h_list[res.id[1]]
    #         if len(new_h) > 0:
    #             new_h = res[new_h[0]]
    #             new_h.id = new_h.fullname = new_h.name = 'H'
    #             res.child_list.append(new_h)
    #             res.child_dict['H'] = new_h

    # Reading PCS Values
    pcs_data_exp = dataparse.read_pcs("../data_files/s57c_pcsexp.npc", lambda x: x[1] == 'H')

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    # Set the starting position to an atom close to the metal
    m_start = metal.Metal(position=prot[0][' '][57]['CB'].position)

    # Calculate an initial tensor from an SVD grid search
    m_guess, calc, q_fac = fit.svd_gridsearch_fit_metal_from_pcs([m_start], [parsed_data], radius=10, points=10)

    # Refine the tensor using non-linear regression
    m_fit, calc, q_fac = fit.nlr_fit_metal_from_pcs(m_guess, [parsed_data])
    m_fit[0].save('s57c_1ubq_H_PCS_tensor.txt')
    print(m_guess, calc, q_fac)


def fit_rotamer():
    # Load the PDB file
    prot = protein.load_pdb('../data_files/1UBQ_H.pdb')

    # Reading PCS Values
    pcs_data_exp = dataparse.read_pcs("../data_files/s57c_pcsexp.npc")

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    # Define an initial tensor
    chi_tensor = {"axrh": (-9.182E-32, -1.495E-32),
                  "position": (21.142E-10, 11.562E-10, 10.689E-10),
                  "eulers": (0.95579466, 1.05276515, 2.01690248)}
    mtl = metal.Metal(**chi_tensor)

    pcs_fitter = fit.PCSToRotamer(prot[0], mtl, parsed_data)

    trk = Tracker(pcs_fitter.model[' '][6])
    # trk.print_atom_linkage()
    trk.save_atom_coords()

    bins = np.arange(0, 16)
    bins[13] = 12
    bins[15] = 13

    pr = cProfile.Profile()

    def staggered_positions_search():
        pr.clear()
        pr.enable()
        # Grid search 5 degrees about each staggered position
        min_pcs = pcs_fitter.run_staggered_positions_search(' ', 6, 3*0.174533, 7, bins=bins)
        pr.disable()
        pr.print_stats(sort="cumtime")
        print(min_pcs)

    def grid_search():
        pr.clear()
        pr.enable()
        rot_param = np.array([[0, -2 / 9 * np.pi, 9]] * 4)
        pcs_fitter.set_rotation_parameter(' ', 6, rot_param, bins)
        prot, result1 = pcs_fitter.run_grid_search(top_n=2, structure='grid_search_struct')
        pr.disable()
        pr.print_stats(sort="cumtime")
        print(result1)
        save_structure(prot, 'grid_search_result.pdb')

    def pairwise_grid_search():
        pr.clear()
        pr.enable()
        rot_param = np.array([[0, -2 / 6 * np.pi, 6]] * 4)
        pcs_fitter.set_rotation_parameter(' ', 6, rot_param, bins)
        prot, result2 = pcs_fitter.run_pairwise_grid_search(top_n=2, structure='pairwise_grid_search_struct')
        pr.disable()
        pr.print_stats(sort="cumtime")
        print(result2)
        save_structure(prot, 'pairwise_grid_search_result.pdb')

    # staggered_positions_search()
    # grid_search()
    pairwise_grid_search()

    trk.save_atom_coords()
    trk.print_atom_coords()


if __name__ == '__main__':
    # fit_metal()
    fit_rotamer()
