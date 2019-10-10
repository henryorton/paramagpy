import numpy as np

from paramagpy import fit, protein, metal, dataparse

prot = protein.load_pdb('/home/go/Workspace/paramagpy/examples/data_files/1UBQ_H.pdb')
pcs_data_exp = dataparse.read_pcs('/home/go/Workspace/paramagpy/examples/data_files/s57c_pcsexp.npc')
parsed_data = prot.parse(pcs_data_exp)

chi_tensor = {"axrh": (-9.182E-32, -1.495E-32),
              "position": (21.142E-10, 11.562E-10, 10.689E-10),
              "eulers": (0.95579466, 1.05276515, 2.01690248)}
mtl = metal.Metal(**chi_tensor)

pcs_fitter = fit.PCSToRotamer(prot[0], mtl, parsed_data)

chain = ' '
res = 6
delta = 0.174533
steps = 5
bins = np.arange(0, 16)
bins[13] = 12
bins[15] = 13


def rad(deg): return (deg / 180) * np.pi


full_rot = np.array([rad(0), -rad(40), 9])

stag_pos = np.array(
    [[rad(60) - delta, rad(60) + delta, steps], [rad(180) - delta, rad(-180) + delta, steps],
     [rad(-60) - delta, rad(-60) + delta, steps]])
min_pcs = (float('inf'), 0, None)
for i in range(3 ** 3):
    _i = [int(c) for c in np.base_repr(i, 3, 3)[-3:]]
    rot_param = np.vstack([full_rot, np.array(stag_pos[_i])])
    pcs_fitter.set_rotation_parameter(chain, res, rot_param, bins)
    result = pcs_fitter.run_grid_search()
    min_pcs = min((-result[chain][res][0][0], min_pcs[1] + 1, result), min_pcs)

print(min_pcs)
