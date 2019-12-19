import bisect
import warnings

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

from paramagpy import protein, metal, fit, dataparse


def rad(deg): return (deg / 180) * np.pi


def fit_rotamer(pdb, res_no, res_name, lvls, bins, bins2, steps):
    # Load the PDB file
    res = str(res_no) + res_name

    # prot = protein.load_pdb('../data_files/' + pdb + '/' + res + '_H.pdb')
    prot = protein.load_pdb('../data_files/' + pdb + '/mut_H.pdb')

    # Reading PCS Values
    pcs_data_exp = dataparse.read_pcs('../data_files/' + pdb + '/' + res + '_pcsexp.npc', err=False)

    # Associate PCS data with atoms of the PDB
    parsed_data = prot.parse(pcs_data_exp)

    # RDC Data
    rdc_data = dataparse.read_rdc('../data_files/' + pdb + '/' + res + '_rdcexp.rdc')

    # CCR Data
    ccr_data = dataparse.read_ccr('../data_files/' + pdb + '/' + res + '_ccrexp.ccr')

    mtl = metal.load_tensor('../data_files/' + pdb + '/' + res + '_Local_Chi_Tensor.txt')
    # mtl = metal.load_tensor('../data_files/' + pdb + '/Global_Chi_Tensor_Metal.txt')

    pcs_fitter = fit.PCSToRotamer(prot[0], mtl, parsed_data)
    x1, x2 = (steps,) * 2
    # Get dihedral angles in PDB file
    x1_i, x2_i = prot[0]['A'][res_no].get_dihedral_full() * 180 / np.pi

    def find_angle(chi_list, chi_in):
        return bisect.bisect_left(chi_list, chi_in)

    def rdc_fit_func(residue, rdcs, _bins=None):
        if _bins is not None:
            _bins = np.zeros(len(rdcs), dtype=np.int64)

        _calc, _exp = np.empty(len(rdcs)), np.empty(len(rdcs))

        for i, _k in enumerate(rdcs):
            k = [x for x in _k]
            if k[0][0] == residue.id[1] and k[1][0] == residue.id[1]:
                _a1, _a2 = residue[k[0][1]], residue[k[1][1]]
                _calc[i] = residue._metal.atom_rdc(_a1, _a2)
                _exp[i] = rdcs[_k][0]
            else:
                warnings.warn('RDC for ' + k + ' doesn\'t belong to ' + residue)
                _calc[i], _exp[i] = 0, 0

        rdc_calc_binned = np.bincount(_bins, weights=_calc)
        rdc_exp_binned = np.bincount(_bins, weights=_exp)
        _bin_c = np.bincount(_bins)
        _bin_c[_bin_c == 0] = 1
        _x = (rdc_exp_binned - rdc_calc_binned) / _bin_c

        rdc_dist = np.linalg.norm(_x) / len(_x)

        return -rdc_dist, np.array(residue._dihedral_full), _calc

    def ccr_fit_func(residue, ccrs, _bins=None):
        if _bins is not None:
            _bins = np.zeros(len(ccrs), dtype=np.int64)

        _calc, _exp = np.empty(len(ccrs)), np.empty(len(ccrs))

        for i, _k in enumerate(ccrs):
            k = [x for x in _k]
            if k[0][0] == residue.id[1] and k[1][0] == residue.id[1]:
                _a1, _a2 = residue[k[0][1]], residue[k[1][1]]
                _calc[i] = residue._metal.atom_ccr(_a1, _a2)
                _exp[i] = ccrs[_k][0]
            else:
                warnings.warn('CCR for ' + k + ' doesn\'t belong to ' + residue)
                _calc[i], _exp[i] = 0, 0

        ccr_calc_binned = np.bincount(_bins, weights=_calc)
        ccr_exp_binned = np.bincount(_bins, weights=_exp)
        _bin_c = np.bincount(_bins)
        _bin_c[_bin_c == 0] = 1
        _x = (ccr_exp_binned - ccr_calc_binned) / _bin_c

        ccr_dist = np.linalg.norm(_x) / len(_x)

        return -ccr_dist, np.array(residue._dihedral_full), _calc

    def grid_search(**kwargs):
        # Set precision for chi angles (steps / full rotation)
        rot_param = np.array(
            [[0, -(2 / x1) * np.pi, x1], [0, - (2 / x2) * np.pi, x2]])
        pcs_fitter.set_rotation_parameter('A', res_no, rot_param, bins)
        # pcs_fitter_h.set_rotation_parameter('A', res_no, rot_param, bins_h)

        # Run the grid search tool
        result1 = pcs_fitter.run_grid_search(top_n=-1, **kwargs)

        # result1_r = pcs_fitter.run_grid_search(top_n=-1, racs=False)

        # save_structure(prot1, 'grid_search_result_1.pdb')

        # Process the results
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

        return x, y, z

    # PCSs
    _x, _y, z1 = grid_search()

    # RDCs
    _x, _y, z2 = grid_search(fit_func=rdc_fit_func, rdcs=rdc_data, _bins=bins2)

    # CCRs
    _x, _y, z3 = grid_search(fit_func=ccr_fit_func, ccrs=ccr_data, _bins=bins2)

    z1 = (z1 - np.amin(z1)) / (np.amax(z1) - np.amin(z1))
    z2 = (z2 - np.amin(z2)) / (np.amax(z2) - np.amin(z2))
    z3 = (z3 - np.amin(z3)) / (np.amax(z3) - np.amin(z3))
    # z1 = np.zeros(z1.shape)

    _z = np.sqrt(np.sum(np.array([z1, z2, z3])**2, axis=0))

    # Distance contour
    fig, ax = plt.subplots(1, 5)
    gs = gridspec.GridSpec(1, 5, width_ratios=[4, 4, 4, 4, 0.2])

    def plot_c(_ax, __z, name):
        ax2 = _ax.contourf(_x, _y, __z, lvls, cmap='RdGy')

        chi_1 = np.linspace(-180, 180, x1 + 1)[1:]
        chi_2 = np.linspace(-180, 180, x2 + 1)[1:]
        # x1_oh, x2_oh = result2['A'][res_no][0][1] * 180 / np.pi
        chi_1_i, chi_2_i = find_angle(chi_1, x1_i), find_angle(chi_2, x2_i)
        chi_1_o, chi_2_o = np.unravel_index(np.argmin(__z), __z.shape)
        x1_o, x2_o = chi_1[chi_1_o], chi_2[chi_2_o]

        _ax.autoscale(False)  # Scatter plot can scale the contour plot, we don't want that
        # zorder=1 ensures scatter plot is overlaid on top of contour, and not the other way
        _ax.scatter(x1_i, x2_i, label='Crystal Structure', zorder=1)  # Mark the angles from the crystal structure
        _ax.annotate(f'{x1_i:.2f}, {x2_i:.2f} [{__z[chi_1_i, chi_2_i]:.4f}]', (x1_i, x2_i))
        # ax.scatter(x1_oh, x2_oh, label='Optimum (H only)', zorder=1)  # Mark the optimized value
        _ax.scatter(x1_o, x2_o, label='Optimum', zorder=1)  # Mark the optimized value
        _ax.annotate(f'{x1_o:.2f}, {x2_o:.2f} [{__z[chi_1_o, chi_2_o]:.4f}]', (x1_o, x2_o))

        _ax.legend(bbox_to_anchor=(1, 0), loc="lower right", bbox_transform=fig.transFigure, ncol=2)
        _ax.set_xlabel(r'$\chi_1 (\degree)$')
        _ax.set_ylabel(r'$\chi_2 (\degree)$')
        _ax.set_title(r'$\chi_1, \chi_2$ vs RMSD of ' + name + ' values for ' + res)

        return ax2

    plot_c(plt.subplot(gs[0], aspect=1), z1, 'PCS')
    plot_c(plt.subplot(gs[1], aspect=1), z2, 'RDC')
    plot_c(plt.subplot(gs[2], aspect=1), z3, 'CCR')
    ax2 = plot_c(plt.subplot(gs[3], aspect=1), (_z - np.amin(_z)) / (np.amax(_z) - np.amin(_z)), 'PCS & RDC')

    # Colorbar
    ax = plt.subplot(gs[4])
    plt.colorbar(ax2, cax=ax)
    ax.set_ylabel('RMSD (Norm)')

    plt.tight_layout()
    plt.show()
    # plt.savefig('../data_files/' + pdb + '/plots/racs/' + res + '_RACS_Crystal.png')
    # plt.close(fig)


if __name__ == '__main__':
    aroms = [('1ig5', 13, 'Tyr', 40, np.array([0, 0, 1, 1, 2, 2, 3, 3]), np.array([0, 0])),  # 0
             ('1ig5', 36, 'Phe', 40, np.array([0, 0, 1, 1, 2, 3, 3, 4, 4, 5]), np.array([0, 0, 1, 1, 2])),  # 1
             ('1ig5', 63, 'Phe', 40, np.array([0, 0, 1, 2, 2, 3]), np.array([0, 0, 1])),  # 2
             ('1ig5', 66, 'Phe', 40, np.array([0, 0, 1, 1, 2, 3, 3, 4, 4, 5]), np.array([0, 0, 1])),  # 3
             ('4icb', 13, 'Tyr', 40, np.array([0, 0, 1, 1, 2, 2, 3, 3]), np.array([0, 0])),  # 4
             ('4icb', 36, 'Phe', 40, np.array([0, 0, 1, 1, 2, 3, 3, 4, 4, 5]), np.array([0, 0, 1, 1, 2])),  # 5
             ('4icb', 63, 'Phe', 40, np.array([0, 0, 1, 2, 2, 3]), np.array([0, 0, 1])),  # 6
             ('4icb', 66, 'Phe', 40, np.array([0, 0, 1, 1, 2, 3, 3, 4, 4, 5]), np.array([0, 0, 1]))]  # 7
    run_for = np.arange(5, 6)
    for _r in run_for:
        print(f"Fitting {aroms[_r]}")
        fit_rotamer(*aroms[_r], 72)
