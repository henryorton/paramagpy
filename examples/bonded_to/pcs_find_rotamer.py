import bisect

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBIO
from matplotlib import gridspec

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
    prot = protein.load_pdb('../data_files/1ig5/4icb_H_aligned.pdb')

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
    pcs_data_exp = dataparse.read_pcs('../data_files/1ig5/obtained_pcs_Nov26_2019.npc',
                                      lambda x: (2 <= int(x[0]) <= 41 or 58 <= int(x[0]) <= 70) and x[1][0] == 'H')

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    # Set the starting position to an atom close to the metal
    m_start = metal.Metal(position=np.array((4.3785, 18.0035, 17.1465)) * 1E-10)
    # m_start = metal.Metal(position=prot[0]['A'][('H_ CA', 77, ' ')]['CA'].position)

    # Calculate an initial tensor from an SVD grid search
    m_guess, calc, q_fac = fit.svd_gridsearch_fit_metal_from_pcs([m_start], [parsed_data], radius=0.8, points=10)

    # Refine the tensor using non-linear regression
    m_fit, calc, q_fac = fit.nlr_fit_metal_from_pcs(m_guess, [parsed_data])
    m_guess[0].save("../data_files/4icb/Global_Chi_Tensor_Metal.txt")
    print(m_guess, calc, q_fac)


def fit_rotamer(pdb, res_no, res_name, dev, bins, steps):
    # Load the PDB file
    res = str(res_no) + res_name

    # prot = protein.load_pdb('../data_files/' + pdb + '/' + res + '_H.pdb')
    prot = protein.load_pdb('../data_files/' + pdb + '/mut_H.pdb')

    # Reading PCS Values
    pcs_data_exp = dataparse.read_pcs('../data_files/' + pdb + '/' + res + '_pcsexp.npc', err=False)
    # bins_h = bins[:int(len(bins) / 2)]

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    mtl = metal.load_tensor('../data_files/' + pdb + '/' + res + '_Local_Chi_Tensor.txt')
    # mtl = metal.load_tensor('../data_files/' + pdb + '/Global_Chi_Tensor_Metal.txt')

    pcs_fitter = fit.PCSToRotamer(prot[0], mtl, parsed_data)

    def staggered_positions_search():
        n = 10
        # Grid search 5 degrees about each staggered position
        result1 = pcs_fitter.run_staggered_positions_search('A', res_no, 15 / 180 * np.pi, n, bins=bins, top_n=-1)

        # pr.print_stats(sort="cumtime")

        # print(min_pcs)
        x1_i, x2_i = prot[0]['A'][res_no].get_dihedral_full() * 180 / np.pi
        x1_o, x2_o = result1[0][1] * 180 / np.pi

        rmsd, chi_angles, pcs_values = list(zip(*result1))
        chi_angles = np.array(chi_angles).transpose() * 180 / np.pi  # Convert to degrees
        x, y = chi_angles[0], chi_angles[1]
        z = -1 * np.array(rmsd)
        xyz = list(zip(x, y, z, pcs_values))  # Get in image format
        list.sort(xyz)
        x, y, z, p = list(zip(*xyz))
        x = np.array(x).reshape(3 * n, 3 * n)
        y = np.array(y).reshape(3 * n, 3 * n)
        z = np.array(z).reshape(3 * n, 3 * n)

        fig, ax = plt.subplots()
        ax.contourf(x, y, z)
        ax.scatter(x1_i, x2_i, label='Crystal Structure', zorder=1)  # Mark the angles from the crystal structure
        ax.scatter(x1_o, x2_o, label='Optimum', zorder=1)  # Mark the optimized value

        plt.show()

    def grid_search():
        # Get dihedral angles in PDB file
        x1_i, x2_i = prot[0]['A'][res_no].get_dihedral_full() * 180 / np.pi

        # Set precision for chi angles (steps / full rotation)
        x1, x2 = (steps,) * 2
        rot_param = np.array(
            [[0, -(2 / x1) * np.pi, x1], [0, - (2 / x2) * np.pi, x2]])
        pcs_fitter.set_rotation_parameter('A', res_no, rot_param, bins)
        # pcs_fitter_h.set_rotation_parameter('A', res_no, rot_param, bins_h)

        # Run the grid search tool
        result1 = pcs_fitter.run_grid_search(top_n=-1)
        # result1_r = pcs_fitter.run_grid_search(top_n=-1, racs=False)

        # save_structure(prot1, 'grid_search_result_1.pdb')

        # Process the results
        def find_angle(chi_list, chi_in):
            return bisect.bisect_left(chi_list, chi_in)

        chi_1 = np.linspace(-180, 180, x1 + 1)[1:]
        chi_2 = np.linspace(-180, 180, x2 + 1)[1:]
        x1_o, x2_o = result1['A'][res_no][0][1] * 180 / np.pi
        # x1_oh, x2_oh = result2['A'][res_no][0][1] * 180 / np.pi
        chi_60, chi_180, chi_300 = find_angle(chi_1, 60), find_angle(chi_1, 180), find_angle(chi_1, -60)
        chi_1_i, chi_1_o = find_angle(chi_1, x1_i), find_angle(chi_1, x1_o)
        chi_2_i, chi_2_o = find_angle(chi_2, x2_i), find_angle(chi_2, x2_o)

        rmsd, chi_angles, pcs_values = list(zip(*result1['A'][res_no]))
        chi_angles = np.array(chi_angles).transpose() * 180 / np.pi  # Convert to degrees
        x, y = chi_angles[0], chi_angles[1]
        z = -1 * np.array(rmsd)
        xyz = list(zip(x, y, z, pcs_values))  # Get in image format
        list.sort(xyz)
        x, y, z, p = list(zip(*xyz))
        x = np.array(x).reshape(x1, x2)
        y = np.array(y).reshape(x1, x2)
        z = np.array(z).reshape(x1, x2)
        # p = np.array(p).reshape(x1, x2, -1)

        # fig, ax = plt.subplots()
        # ax.plot(chi_1, p[:, 0, 4])
        # ax.plot(chi_1, np.transpose(np.ones(len(chi_1)) * np.array([[0.01], [-0.01]])))

        # PCS distance contour
        fig, ax = plt.subplots(1, 2, figsize=(8.4, 6.5))
        gs = gridspec.GridSpec(1, 2, width_ratios=[4, 0.2])
        ax = plt.subplot(gs[0], aspect=1)
        lvls = int((np.amax(z) - np.amin(z)) / dev)
        ax2 = ax.contourf(x, y, z, lvls, cmap='RdGy')

        ax.autoscale(False)  # Scatter plot can scale the contour plot, we don't want that
        # zorder=1 ensures scatter plot is overlaid on top of contour, and not the other way
        ax.scatter(x1_i, x2_i, label='Crystal Structure', zorder=1)  # Mark the angles from the crystal structure
        ax.annotate(f'{x1_i:.2f}, {x2_i:.2f} [{z[chi_1_i, chi_2_i]:.4f}]', (x1_i, x2_i))
        # ax.scatter(x1_oh, x2_oh, label='Optimum (H only)', zorder=1)  # Mark the optimized value
        ax.scatter(x1_o, x2_o, label='Optimum', zorder=1)  # Mark the optimized value
        ax.annotate(f'{x1_o:.2f}, {x2_o:.2f} [{z[chi_1_o, chi_2_o]:.4f}]', (x1_o, x2_o))

        ax.legend(bbox_to_anchor=(1, 0), loc="lower right", bbox_transform=fig.transFigure, ncol=2)
        ax.set_xlabel(r'$\chi_1 (\degree)$')
        ax.set_ylabel(r'$\chi_2 (\degree)$')
        # ax.set_title(r'$\chi_1, \chi_2$ vs RMSD of PCS values for ' + res)

        # Colorbar
        ax = plt.subplot(gs[1])
        plt.colorbar(ax2, cax=ax)
        ax.set_ylabel('RMSD (ppm)')

        # plt.savefig('../data_files/' + pdb + '/plots/racs/' + res + '_Contour.png')
        # plt.close(fig)

        # from mpl_toolkits.mplot3d import Axes3D
        # fig, ax = plt.subplots()
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_wireframe(x, y, z)

        # Slices at staggered X1
        # ax = plt.subplot(gs[2], aspect=360 / (np.amax(z) - np.amin(z)))
        #
        # ax.set_prop_cycle(cycler(color=['#2ca02c', '#d62728', '#9467bd', '#1f77b4', '#ff7f0e']))
        # ax.plot(chi_2, z[chi_60], label=r'$\chi_1 = 60\degree$')
        # ax.plot(chi_2, z[chi_180], label=r'$\chi_1 = 180\degree$')
        # ax.plot(chi_2, z[chi_300], label=r'$\chi_1 = -60\degree$')
        # ax.plot(chi_2, z[chi_1_i], label=r'$\chi_1 = ' + str(int(x1_i)) + r'\degree$')
        # ax.plot(chi_2, z[chi_1_o], label=r'$\chi_1 = ' + str(int(x1_o)) + r'\degree$')
        #
        # ax.set_ylim(np.amin(z), np.amax(z))
        # ax.set_xlim(-180, 180)
        # # ax.legend()
        # ax.set_xlabel(r'$\chi_2 (\degree)$')
        # ax.set_ylabel('RMSD (ppm)')
        # ax.set_title(r'$\chi_2$ vs RMSD of PCS values for ' + res)

        # X2 effects averaged
        # fig, ax = plt.subplots()
        # ax.plot(chi_1, np.mean(z, axis=1))
        # ax.axvline(x=x1_i)
        # ax.set_xlabel(r'$\chi_1 (\degree)$')
        # ax.set_ylabel('RMSD (ppm)')

        # Slices of individual deviations
        # H @ optimum
        # fig, ax = plt.subplots()
        # bin_count = np.bincount(bins)
        # bin_count[bin_count == 0] = 1
        # p_binned = np.empty((len(chi_2), len(bin_count)))
        # for idx in range(len(chi_2)):
        #     p_binned[idx] = np.bincount(bins, weights=p[find_angle(chi_1, x1_o), idx]) / bin_count

        # for idx, _pcs in enumerate(pcs_data_exp):
        #     if _pcs[1][0] == 'H':
        #         ax.plot(chi_2, p[find_angle(chi_1, x1_o), :, idx], label=_pcs[1])
        #         if bin_count[bins[idx]] > 1:
        #             ax.plot(chi_2, p_binned[:, bins[idx]], label=_pcs[1][:-1] + '*')
        #             bin_count[bins[idx]] = 1
        # ax.legend()
        # ax.set_xlabel(r'$\chi_2 (\degree)$')
        # ax.set_ylabel('Deviation (Calc-Exp)')
        # ax.set_title(r'$\chi_2$ vs Deviation at ' + r'$\chi_1 = $' + str(int(x1_o)) + r'$\degree$ of H PCS values for '
        #              + res)

        # C @ optimum
        # fig, ax = plt.subplots()
        # for idx, _pcs in enumerate(pcs_data_exp):
        #     if _pcs[1][0] == 'C':
        #         ax.plot(chi_2, p[find_angle(chi_1, x1_o), :, idx], label=_pcs[1])
        #         if bin_count[bins[idx]] > 1:
        #             ax.plot(chi_2, p_binned[:, bins[idx]], label=_pcs[1][:-1] + '*')
        #             bin_count[bins[idx]] = 1
        #
        # ax.legend()
        # ax.set_xlabel(r'$\chi_2 (\degree)$')
        # ax.set_ylabel('RACS (ppm)')
        # ax.set_title(r'$\chi_2$ vs RACS at ' + r'$\chi_1 = $' + str(int(x1_o)) + r'$\degree$ for ' + res)
        # if cnt < 1:
        #     ax.clear()
        # plt.savefig('../data_files/' + pdb + '/plots/racs/' + res + '_RACS_Optimum.png')
        # plt.close(fig)

        # H @ Crystal
        # fig, ax = plt.subplots()
        # bin_count = np.bincount(bins)
        # bin_count[bin_count == 0] = 1
        # p_binned = np.empty((len(chi_2), len(bin_count)))
        # for idx in range(len(chi_2)):
        #     p_binned[idx] = np.bincount(bins, weights=p[find_angle(chi_1, x1_i), idx]) / bin_count
        #
        # # for idx, _pcs in enumerate(pcs_data_exp):
        #     if _pcs[1][0] == 'H':
        #         ax.plot(chi_2, p[find_angle(chi_1, x1_i), :, idx], label=_pcs[1])
        #         if bin_count[bins[idx]] > 1:
        #             ax.plot(chi_2, p_binned[:, bins[idx]], label=_pcs[1][:-1] + '*')
        #             bin_count[bins[idx]] = 1
        # ax.legend()
        # ax.set_xlabel(r'$\chi_2 (\degree)$')
        # ax.set_ylabel('Deviation (Calc-Exp)')
        # ax.set_title(r'$\chi_2$ vs Deviation at ' + r'$\chi_1 = $' + str(int(x1_i)) + r'$\degree$ of H PCS values for '
        #              + res)

        # C @ Crystal
        # fig, ax = plt.subplots()
        # for idx, _pcs in enumerate(pcs_data_exp):
        #     if _pcs[1][0] == 'C':
        #         ax.plot(chi_2, p[find_angle(chi_1, x1_i), :, idx], label=_pcs[1])
        #         if bin_count[bins[idx]] > 1:
        #             ax.plot(chi_2, p_binned[:, bins[idx]], label=_pcs[1][:-1] + '*')
        #             bin_count[bins[idx]] = 1
        #
        # ax.legend()
        # ax.set_xlabel(r'$\chi_2 (\degree)$')
        # ax.set_ylabel('RACS (ppm)')
        # ax.set_title(r'$\chi_2$ vs RACS at ' + r'$\chi_1 = $' + str(int(x1_i)) + r'$\degree$ for ' + res)

        plt.tight_layout()
        plt.show()
        # plt.savefig('../data_files/' + pdb + '/plots/racs/' + res + '_RACS_Crystal.png')
        # plt.close(fig)

    def pairwise_grid_search():
        rot_param = np.array([[-1.2922218986820855, -1.2922218986820855, 2], [np.pi / 3, -2 / 3 * np.pi, 2]])
        pcs_fitter.set_rotation_parameter('A', 13, rot_param, bins)
        prot, result2 = pcs_fitter.run_pairwise_grid_search(top_n=2, structure='pairwise_grid_search_struct')
        print(result2)
        save_structure(prot, 'pairwise_grid_search_result.pdb')

    def test():
        # prot = protein.load_pdb('../data_files/4icb/36Phe_fixed.pdb')
        pcs_calc = [None, None]
        res = prot[0]['A'][res_no]
        res.set_dihedral(np.array([-80, 0]) * np.pi / 180)
        save_structure(prot, 'test_1.pdb')
        coord_matrix = np.empty((len(res.child_list), 3))

        for _a in enumerate(res.child_list):
            print(_a)

        for idx, _a in enumerate(res.child_list):
            coord_matrix[idx] = _a.position
        pcs_calc[0] = mtl.fast_pcs(coord_matrix)

        res.set_delta_dihedral(np.array([0, np.pi]))
        save_structure(prot, 'test_2.pdb')
        coord_matrix = np.empty((len(res.child_list), 3))

        for idx, _a in enumerate(res.child_list):
            coord_matrix[idx] = _a.position

        pcs_calc[1] = mtl.fast_pcs(coord_matrix)

        print(pcs_calc)

    # test()
    # staggered_positions_search()
    grid_search()

    # pairwise_grid_search()


if __name__ == '__main__':
    # fit_metal()
    aroms = [('1ig5', 13, 'Tyr', 0.019, np.array([0, 0, 1, 1, 2, 2, 3, 3])),  # 0
             ('1ig5', 36, 'Phe', 0.006, np.array([0, 0, 1, 1, 2, 3, 3, 4, 4, 5])),  # 1
             ('1ig5', 63, 'Phe', 0.008, np.array([0, 0, 1, 2, 2, 3])),  # 2
             ('1ig5', 66, 'Phe', 0.035, np.array([0, 0, 1, 1, 2, 3, 3, 4, 4, 5])),  # 3
             ('4icb', 13, 'Tyr', 0.018, np.array([0, 0, 1, 1, 2, 2, 3, 3])),  # 4
             ('4icb', 36, 'Phe', 0.011, np.array([0, 0, 1, 1, 2, 3, 3, 4, 4, 5])),  # 5
             ('4icb', 63, 'Phe', 0.027, np.array([0, 0, 1, 2, 2, 3])),  # 6
             ('4icb', 66, 'Phe', 0.055, np.array([0, 0, 1, 1, 2, 3, 3, 4, 4, 5]))]  # 7
    run_for = np.arange(1, 2)
    for _r in run_for:
        print(f"Fitting {aroms[_r]}")
        fit_rotamer(*aroms[_r], 72)

    # fit_rotamer('1ig5', 9, 'Ile', 0.019, np.array([0, 0, 0, 1, 2, 2, 2, 3]), 180)
