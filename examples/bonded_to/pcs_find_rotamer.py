import bisect
import cProfile

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBIO
from bonded_to import Tracker
from matplotlib import gridspec, cycler

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
    prot = protein.load_pdb('../data_files/4icb/mut_H.pdb')

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
    pcs_data_exp = dataparse.read_pcs('../data_files/4icb/36Phe_Local_PCS.npc')

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    # Set the starting position to an atom close to the metal
    m_start = metal.Metal(position=prot[0]['A'][('H_ CA', 77, ' ')]['CA'].position)

    # Calculate an initial tensor from an SVD grid search
    m_guess, calc, q_fac = fit.svd_gridsearch_fit_metal_from_pcs([m_start], [parsed_data], radius=0.5, points=10)

    # Refine the tensor using non-linear regression
    # m_fit, calc, q_fac = fit.nlr_fit_metal_from_pcs(m_guess, [parsed_data])
    m_guess[0].save("../data_files/4icb/13Tyr_Local_Chi_Tensor.txt")
    print(m_guess, calc, q_fac)


def fit_rotamer(pdb, res_no, res_name, dev, bins, steps):
    # Load the PDB file
    res = str(res_no) + res_name

    # prot = protein.load_pdb('../data_files/' + pdb + '/' + res + '_H.pdb')
    prot = protein.load_pdb('../data_files/' + pdb + '/mut_H.pdb')

    # Reading PCS Values
    pcs_data_exp = dataparse.read_pcs('../data_files/' + pdb + '/' + res + '_pcsexp.npc')

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    # Define an initial tensor
    # chi_tensor = {"axrh": (9.777E-32, 5.653E-32),
    #               "position": (4.124E-10, 18.763E-10, 17.169E-10),
    #               "eulers": (2.53730731, 1.87498976, 0.91317372)}
    # mtl = metal.Metal(**chi_tensor)
    # mtl.save('../data_files/1ig5/' + res + '_Local_Chi_Tensor.txt')

    mtl = metal.load_tensor('../data_files/' + pdb + '/' + res + '_Local_Chi_Tensor.txt')

    pcs_fitter = fit.PCSToRotamer(prot[0], mtl, parsed_data)

    trk = Tracker(pcs_fitter.model['A'][res_no])
    # trk.print_atom_linkage()
    trk.save_atom_coords()

    pr = cProfile.Profile()

    def staggered_positions_search():
        pr.clear()
        pr.enable()
        # Grid search 5 degrees about each staggered position
        min_pcs = pcs_fitter.run_staggered_positions_search('A', 13, 3 * 0.174533, 7, bins=bins)
        pr.disable()
        pr.print_stats(sort="cumtime")
        print(min_pcs)

    def grid_search():
        pr.clear()
        pr.enable()

        # Get dihedral angles in PDB file
        x1_i, x2_i = prot[0]['A'][res_no].get_dihedral_full() * 180 / np.pi

        # Set precision for chi angles (steps / full rotation)
        x1, x2 = (steps,) * 2
        rot_param = np.array(
            [[0, -(2 / x1) * np.pi, x1], [0, - (2 / x2) * np.pi, x2]])
        pcs_fitter.set_rotation_parameter('A', res_no, rot_param, bins)

        # Run the grid search tool
        result1 = pcs_fitter.run_grid_search(top_n=-1)
        # prot1, result1 = pcs_fitter.run_grid_search(top_n=1, structure='grid_search_result.pdb')

        pr.disable()
        pr.print_stats(sort="cumtime")

        # print(result1)
        # save_structure(prot1, 'grid_search_result_1.pdb')

        # Process the results
        def find_angle(chi_list, chi_in):
            return bisect.bisect_left(chi_list, chi_in)

        chi_1 = np.linspace(-180, 180, x1 + 1)[1:]
        chi_2 = np.linspace(-180, 180, x2 + 1)[1:]
        x1_o, x2_o = result1['A'][res_no][0][1] * 180 / np.pi
        chi_60, chi_180, chi_300 = find_angle(chi_1, 60), find_angle(chi_1, 180), find_angle(chi_1, -60)
        chi_1_i, chi_1_o = find_angle(chi_1, x1_i), find_angle(chi_1, x1_o)
        chi_2_i, chi_2_o = find_angle(chi_2, x2_i), find_angle(chi_2, x2_o)

        rmsd, chi_angles, pcs_values, rdc_values = list(zip(*result1['A'][res_no]))
        chi_angles = np.array(chi_angles).transpose() * 180 / np.pi  # Convert to degrees
        x, y = chi_angles[0], chi_angles[1]
        z = -1 * np.array(rmsd)
        xyz = list(zip(x, y, z, pcs_values, rdc_values))  # Get in image format
        list.sort(xyz)
        x, y, z, p, r = list(zip(*xyz))
        x = np.array(x).reshape(x1, x2)
        y = np.array(y).reshape(x1, x2)
        z = np.array(z).reshape(x1, x2)
        p = np.array(p).reshape(x1, x2, -1)
        r = np.array(r).reshape(x1, x2, -1)

        # fig, ax = plt.subplots()
        # ax.plot(chi_1, p[:, 0, 4])
        # ax.plot(chi_1, np.transpose(np.ones(len(chi_1)) * np.array([[0.01], [-0.01]])))

        # RDC values contour
        # fig, ax = plt.subplots()
        # r1 = (r[:, :, 0] + r[:, :, 1]) / 2
        # ax.contour(x, y, r1, 20, colors='black')
        # plt.imshow(np.transpose(r1), extent=[x[0][0], x[-1][-1], y[0][0], y[-1][-1]], origin='lower', cmap='RdGy',
        #            alpha=0.5)
        # plt.colorbar()
        # ax.set_xlabel(r'$\chi_1 (\degree)$')
        # ax.set_ylabel(r'$\chi_2 (\degree)$')
        # ax.set_title(r'$\chi_1, \chi_2$ vs CE* - HE* RDC values for ' + res)
        #
        # fig, ax = plt.subplots()
        # r1 = (r[:, :, 2] + r[:, :, 3]) / 2
        # ax.contour(x, y, r1, 20, colors='black')
        # plt.imshow(np.transpose(r1), extent=[x[0][0], x[-1][-1], y[0][0], y[-1][-1]], origin='lower', cmap='RdGy',
        #            alpha=0.5)
        # plt.colorbar()
        # ax.set_xlabel(r'$\chi_1 (\degree)$')
        # ax.set_ylabel(r'$\chi_2 (\degree)$')
        # ax.set_title(r'$\chi_1, \chi_2$ vs CD* - HD* RDC values for ' + res)

        result2 = pcs_fitter.run_pairwise_grid_search(result=result1, top_n=5)

        # PCS distance contour
        fig, ax = plt.subplots(1, 3, figsize=(16.4, 6.5))
        gs = gridspec.GridSpec(1, 3, width_ratios=[4, 0.2, 4])
        ax = plt.subplot(gs[0], aspect=1)
        lvls = int((np.amax(z) - np.amin(z)) / dev)
        ax2 = ax.contourf(x, y, z, lvls, cmap='RdGy')

        ax.autoscale(False)  # Scatter plot can scale the contour plot, we don't want that
        # zorder=1 ensures scatter plot is overlaid on top of contour, and not the other way
        ax.scatter(x1_i, x2_i, label='Crystal Structure', zorder=1)  # Mark the angles from the crystal structure
        # ax.annotate(f'{x1_i:.2f}, {x2_i:.2f} [{z[chi_1_i, chi_2_i]:.4f}]', (x1_i, x2_i))
        ax.scatter(x1_o, x2_o, label='Optimum', zorder=1)  # Mark the optimized value
        # ax.annotate(f'{x1_o:.2f}, {x2_o:.2f} [{z[chi_1_o, chi_2_o]:.4f}]', (x1_o, x2_o))

        # ax.legend(bbox_to_anchor=(1, 0), loc="lower right", bbox_transform=fig.transFigure, ncol=2)
        ax.set_xlabel(r'$\chi_1 (\degree)$')
        ax.set_ylabel(r'$\chi_2 (\degree)$')
        # ax.set_title(r'$\chi_1, \chi_2$ vs RMSD of PCS values for ' + res)

        # Colorbar
        ax = plt.subplot(gs[1])
        plt.colorbar(ax2, cax=ax)
        ax.set_ylabel('RMSD (ppm)')

        # from mpl_toolkits.mplot3d import Axes3D
        # fig, ax = plt.subplots()
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_wireframe(x, y, z)

        # Slices at staggered X1
        ax = plt.subplot(gs[2], aspect=360 / (np.amax(z) - np.amin(z)))

        ax.set_prop_cycle(cycler(color=['#2ca02c', '#d62728', '#9467bd', '#1f77b4', '#ff7f0e']))
        ax.plot(chi_2, z[chi_60], label=r'$\chi_1 = 60\degree$')
        ax.plot(chi_2, z[chi_180], label=r'$\chi_1 = 180\degree$')
        ax.plot(chi_2, z[chi_300], label=r'$\chi_1 = -60\degree$')
        ax.plot(chi_2, z[chi_1_i], label=r'$\chi_1 = ' + str(int(x1_i)) + r'\degree$')
        ax.plot(chi_2, z[chi_1_o], label=r'$\chi_1 = ' + str(int(x1_o)) + r'\degree$')

        ax.set_ylim(np.amin(z), np.amax(z))
        ax.set_xlim(-180, 180)
        # ax.legend()
        ax.set_xlabel(r'$\chi_2 (\degree)$')
        ax.set_ylabel('RMSD (ppm)')
        # ax.set_title(r'$\chi_2$ vs RMSD of PCS values for ' + res)

        # Slices of individual deviations
        # H @ optimum
        # fig, ax = plt.subplots()
        # bin_count = np.bincount(bins)
        # bin_count[bin_count == 0] = 1
        # p_binned = np.empty((len(chi_2), len(bin_count)))
        # for idx in range(len(chi_2)):
        #     p_binned[idx] = np.bincount(bins, weights=p[find_angle(chi_1, x1_o), idx]) / bin_count
        #
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
        # cnt = 0
        # for idx, _pcs in enumerate(pcs_data_exp):
        #     if _pcs[1][0] == 'C':
        #         cnt += 1
        #         ax.plot(chi_2, p[find_angle(chi_1, x1_o), :, idx], label=_pcs[1])
        #         if bin_count[bins[idx]] > 1:
        #             ax.plot(chi_2, p_binned[:, bins[idx]], label=_pcs[1][:-1] + '*')
        #             bin_count[bins[idx]] = 1
        #
        # ax.legend()
        # ax.set_xlabel(r'$\chi_2 (\degree)$')
        # ax.set_ylabel('Deviation (Calc-Exp)')
        # ax.set_title(r'$\chi_2$ vs Deviation at ' + r'$\chi_1 = $' + str(int(x1_o)) + r'$\degree$ of C PCS values for '
        #              + res)
        # if cnt < 1:
        #     ax.clear()

        # H @ Crystal
        # fig, ax = plt.subplots()
        # bin_count = np.bincount(bins)
        # bin_count[bin_count == 0] = 1
        # p_binned = np.empty((len(chi_2), len(bin_count)))
        # for idx in range(len(chi_2)):
        #     p_binned[idx] = np.bincount(bins, weights=p[find_angle(chi_1, x1_i), idx]) / bin_count
        #
        # for idx, _pcs in enumerate(pcs_data_exp):
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
        # cnt = 0
        # for idx, _pcs in enumerate(pcs_data_exp):
        #     if _pcs[1][0] == 'C':
        #         cnt += 1
        #         ax.plot(chi_2, p[find_angle(chi_1, x1_i), :, idx], label=_pcs[1])
        #         if bin_count[bins[idx]] > 1:
        #             ax.plot(chi_2, p_binned[:, bins[idx]], label=_pcs[1][:-1] + '*')
        #             bin_count[bins[idx]] = 1
        #
        # ax.legend()
        # ax.set_xlabel(r'$\chi_2 (\degree)$')
        # ax.set_ylabel('Deviation (Calc-Exp)')
        # ax.set_title(r'$\chi_2$ vs Deviation at ' + r'$\chi_1 = $' + str(int(x1_i)) + r'$\degree$ of C PCS values for '
        #              + res)
        # if cnt < 1:
        #     ax.clear()

        plt.tight_layout()
        plt.show()
        # plt.savefig('../data_files/' + pdb + '/plots/' + res + '.png')

    def pairwise_grid_search():
        pr.clear()
        pr.enable()
        rot_param = np.array([[-1.2922218986820855, -1.2922218986820855, 2], [np.pi / 3, -2 / 3 * np.pi, 2]])
        pcs_fitter.set_rotation_parameter('A', 13, rot_param, bins)
        prot, result2 = pcs_fitter.run_pairwise_grid_search(top_n=2, structure='pairwise_grid_search_struct')
        pr.disable()
        pr.print_stats(sort="cumtime")
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

    trk.save_atom_coords()
    trk.print_atom_coords()


if __name__ == '__main__':
    # fit_metal()
    aroms = [('1ig5', 13, 'Tyr', 0.019, np.array([0, 0, 1, 1])),  # 0
             ('1ig5', 36, 'Phe', 0.006, np.array([0, 0, 1, 1, 2])),  # 1
             ('1ig5', 63, 'Phe', 0.008, np.array([0, 0, 1])),  # 2
             ('1ig5', 66, 'Phe', 0.035, np.array([0, 0, 1, 1, 2])),  # 3
             ('4icb', 13, 'Tyr', 0.018, np.array([0, 0, 1, 1])),  # 4
             ('4icb', 36, 'Phe', 0.011, np.array([0, 0, 1, 1, 2])),  # 5
             ('4icb', 63, 'Phe', 0.027, np.array([0, 0, 1])),  # 6
             ('4icb', 66, 'Phe', 0.055, np.array([0, 0, 1, 1, 2]))]  # 7
    run_for = np.arange(1, 2)
    for _r in run_for:
        print(f"Fitting {aroms[_r]}")
        fit_rotamer(*aroms[_r], 72)
