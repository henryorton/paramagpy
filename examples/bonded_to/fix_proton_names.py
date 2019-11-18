from Bio.PDB import PDBIO

from paramagpy import protein


def main():
    def save_structure(structure, filename):
        # Save structure in PDB format

        io = PDBIO()
        io.set_structure(structure)
        io.save(filename)

    prot = protein.load_pdb('/home/go/Workspace/paramagpy_extra/4icbH_mut_1.pdb')

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

    for res in prot[0]['A']:
        for _a in res:
            if not _a.disordered_flag and _a.element == 'C':
                print(f"{_a} : {_a.bonded_to()}")

    save_structure(prot, '/home/go/Workspace/paramagpy_extra/4icbH_mut_H.pdb')


def test():
    prot = protein.load_pdb('/home/go/Workspace/paramagpy_extra/4icbH_mut_H.pdb')
    res = prot[0]['A'][36]
    print(res['HE2'] - res['H'])
    # print(res['HB3'] - res['H01'])
    # print(res['HB2'] - res['H03'])
    # print(res['HB2'] - res['H01'])


if __name__ == '__main__':
    # main()
    test()
