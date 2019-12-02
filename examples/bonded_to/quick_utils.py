import numpy as np
from Bio.PDB import PDBIO

from paramagpy import protein


def save_structure(structure, filename):
    # Save structure in PDB format

    io = PDBIO()
    io.set_structure(structure)
    io.save(filename)


def fix_proton_names():
    prot = protein.load_pdb('/home/go/Downloads/mut_H.pdb')

    for res in prot[0]['A']:
        for _a in res:
            if not _a.disordered_flag and _a.element == 'C' and _a.name not in res.back_bone:
                _protons = list(filter(lambda x: x.element == 'H', _a.bonded_to()))
                if len(_protons) == 0:
                    continue
                if len(_protons) == 1:
                    _proton = _protons[0]
                    _name = 'H' + _a.id[1:]
                    _proton.id = _name
                    _proton.name = _name
                    _fid = list(_proton.full_id)
                    _fid[4] = (_name, _fid[4][1])
                    _proton.full_id = tuple(_fid)
                    _proton.fullname = _name.rjust(4)
                    continue
                for idx, _proton in enumerate(_protons):
                    _name = 'H' + _a.id[1:] + str(3 - idx)
                    _proton.id = _name
                    _proton.name = _name
                    _fid = list(_proton.full_id)
                    _fid[4] = (_name, _fid[4][1])
                    _proton.full_id = tuple(_fid)
                    _proton.fullname = _name.rjust(4)

            if not _a.disordered_flag and _a.element == 'N':
                _protons = list(filter(lambda x: x.element == 'H', _a.bonded_to()))
                if len(_protons) != 1:
                    continue
                _proton = _protons[0]
                _name = 'H'
                _proton.id = _name
                _proton.name = _name
                _fid = list(_proton.full_id)
                _fid[4] = (_name, _fid[4][1])
                _proton.full_id = tuple(_fid)
                _proton.fullname = _name.rjust(4)

    for res in prot[0]['A']:
        for _a in res:
            if not _a.disordered_flag and _a.element == 'C':
                print(f"{_a} : {_a.bonded_to()}")

    save_structure(prot, '/home/go/Workspace/paramagpy_extra/mut_H.pdb')


def nearest_backbone_protons():
    prot = protein.load_pdb('../data_files/4icb/mut_H.pdb')
    aroms = [13, 36, 63, 66]
    h_dict = {}
    ring = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
    centers = {}
    for _a in aroms:
        pos = np.mean(np.array([prot[0]['A'][_a][x].coord for x in ring]), axis=0)
        centers[_a] = pos

    for res in prot[0]['A']:
        if res.has_id('H'):
            h_dict[res.id[1]] = res['H']

    result = {}
    for src_res in aroms:
        result[src_res] = []
        for dst_res, dst_atom in h_dict.items():
            if src_res != dst_res:
                result[src_res].append((np.linalg.norm(centers[src_res] - dst_atom.coord), dst_res))
        result[src_res].sort()

    print(result[66])


def fix_arom_ring():
    def uv(v):
        return v / np.linalg.norm(v)

    prot = protein.load_pdb('../data_files/4icb/mut_H.pdb')
    arom = prot[0]['A'][36]

    prot = protein.load_pdb('../data_files/4icb/phe.pdb')
    phe = prot[0][' '][0]

    delta = arom['CB'].coord - phe['CB'].coord
    for _a in phe:
        _a.coord = _a.coord + delta

    v1 = phe['CA'].coord - phe['CB'].coord
    v2 = arom['CA'].coord - arom['CB'].coord
    uv1, uv2 = uv(v1), uv(v2)
    rot_axis = np.cross(uv1, uv2)
    theta = np.arccos(np.clip(uv1.dot(uv2), -1, 1))
    q = phe._CustomResidue__rot_quat(theta, rot_axis)
    atoms_position = {x: 999 for x in phe}
    phe._CustomResidue__rotate_atoms(atoms_position, 0, q, phe['CB'].coord)

    for _bb_a in ['HA', 'N', 'C']:
        phe[_bb_a].coord = arom[_bb_a].coord

    phe.set_dihedral(arom.get_dihedral_full())

    for _a in phe:
        print(_a.id, _a.coord, arom[_a.id].coord, np.linalg.norm(_a.coord - arom[_a.id].coord))

    save_structure(prot, '../data_files/4icb/36Phe_fixed.pdb')


if __name__ == '__main__':
    # fix_proton_names()
    nearest_backbone_protons()
    # fix_arom_ring()
