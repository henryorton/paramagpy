import bisect
import cProfile

import matplotlib.pyplot as plt
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
    prot = protein.load_pdb('/home/go/Workspace/paramagpy_extra/4icbH_mut_H.pdb')

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
    pcs_data_exp = dataparse.read_pcs("/home/go/Workspace/Calbindin_Assigments/obtained_pcs.npc")

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    # Set the starting position to an atom close to the metal
    m_start = metal.Metal(position=prot[0]['A'][56]['CB'].position)

    # Calculate an initial tensor from an SVD grid search
    m_guess, calc, q_fac = fit.svd_gridsearch_fit_metal_from_pcs([m_start], [parsed_data], radius=10, points=10)

    # Refine the tensor using non-linear regression
    m_fit, calc, q_fac = fit.nlr_fit_metal_from_pcs(m_guess, [parsed_data])
    m_fit[0].save('/home/go/Workspace/Calbindin_Assigments/calb_PCS_tensor_NH_56.txt')
    print(m_guess, calc, q_fac)


def fit_rotamer():
    # Load the PDB file
    res_no, res_name = 36, 'Phe'
    res = str(res_no) + res_name

    prot = protein.load_pdb('../data_files/4icb/' + res + '_H.pdb')

    # Reading PCS Values
    pcs_data_exp = dataparse.read_pcs('../data_files/4icb/' + res + '_pcsexp.npc')

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    # Define an initial tensor
    chi_tensor = {"axrh": (-7.433E-32, -4.221E-32),
                  "position": (25.908E-10, 10.392E-10, 7.387E-10),
                  "eulers": (2.2780561, 2.54488204, 0.70586351)}
    mtl = metal.Metal(**chi_tensor)

    pcs_fitter = fit.PCSToRotamer(prot[0], mtl, parsed_data)

    trk = Tracker(pcs_fitter.model['A'][res_no])
    # trk.print_atom_linkage()
    trk.save_atom_coords()

    bins = np.array([3, 3, 4, 4, 5])

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
        x1, x2 = 36, 36
        rot_param = np.array(
            [[0, -(2 / x1) * np.pi, x1], [0, - (2 / x2) * np.pi, x2]])
        pcs_fitter.set_rotation_parameter('A', res_no, rot_param, bins)

        # Run the grid search tool
        result1 = pcs_fitter.run_grid_search(top_n=-1)
        # prot1, result1 = pcs_fitter.run_grid_search(top_n=1, structure='grid_search_result.pdb')

        pr.disable()
        pr.print_stats(sort="cumtime")

        # print(result1)
        # save_structure(prot1, 'grid_search_result.pdb')

        # Process the results
        rmsd, chi_angles, pcs_values = list(zip(*result1['A'][res_no]))
        chi_angles = np.array(chi_angles).transpose() * 180 / np.pi  # Convert to degrees
        x, y = chi_angles[0], chi_angles[1]
        z = -1 * np.array(rmsd)
        xyz = list(zip(x, y, z))  # Get in image format
        list.sort(xyz)
        x, y, z = list(zip(*xyz))
        x = np.array(x).reshape(x1, x2)
        y = np.array(y).reshape(x1, x2)
        z = np.array(z).reshape(x1, x2)

        fig, ax = plt.subplots()
        ax.contour(x, y, z, 3, colors='black')
        plt.imshow(np.transpose(z), extent=[x[0][0], x[-1][-1], y[0][0], y[-1][-1]], origin='lower', cmap='RdGy',
                   alpha=0.5)
        plt.colorbar()
        # ax.clabel(CS, fontsize=8, inline=True)

        ax.autoscale(False)  # Scatter plot can scale the contour plot, we don't want that
        # zorder=1 ensures scatter plot is overlaid on top of contour, and not the other way
        ax.scatter(x1_i, x2_i, label='Crystal Structure', zorder=1)  # Mark the angles from the crystal structure
        ax.scatter(*(result1['A'][res_no][0][1] * 180 / np.pi), label='Optimum', zorder=1)  # Mark the optimized value

        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_wireframe(x, y, z)

        ax.legend(bbox_to_anchor=(1, 0), loc="lower right", bbox_transform=fig.transFigure, ncol=2)
        ax.set_xlabel(r'$\chi_1 (\degree)$')
        ax.set_ylabel(r'$\chi_2 (\degree)$')
        ax.set_title(r'$\chi_1, \chi_2$ vs RMSD of PCS values for ' + res)

        # Slices at staggered X1
        fig, ax = plt.subplots()

        def find_angle(chi_list, chi_in):
            return bisect.bisect_left(chi_list, chi_in)

        chi_1 = np.linspace(-180, 180, x1 + 1)[1:]
        chi_2 = np.linspace(-180, 180, x2 + 1)[1:]
        chi_60, chi_180, chi_300 = find_angle(chi_1, 60), find_angle(chi_1, 180), find_angle(chi_1, -60)
        ax.plot(chi_2, z[chi_60], label=r'$\chi_1 = 60\degree$')
        ax.plot(chi_2, z[chi_180], label=r'$\chi_1 = 180\degree$')
        ax.plot(chi_2, z[chi_300], label=r'$\chi_1 = -60\degree$')
        ax.legend()
        ax.set_xlabel(r'$\chi_2 (\degree)$')
        ax.set_ylabel('RMSD')
        ax.set_title(r'$\chi_2$ vs RMSD of PCS values for ' + res)

        plt.show()

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
        pcs_calc = [None, None]
        res = prot[0]['A'][res_no]
        res.set_dihedral(np.array([-60, 95]) * np.pi/180)
        save_structure(prot, 'test_1.pdb')
        coord_matrix = np.empty((len(res.child_list), 3))

        for _a in enumerate(res.child_list):
            print(_a)

        for idx, _a in enumerate(res.child_list):
            coord_matrix[idx] = _a.position
        pcs_calc[0] = mtl.fast_pcs(coord_matrix)

        res.set_delta_dihedral(np.array([0, np.pi/2]))
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
    fit_rotamer()
