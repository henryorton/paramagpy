# PyMOL macro for loading tensor isosurface from paramagpy
# ax    | 1E-32 m^3 :    -8.152
# rh    | 1E-32 m^3 :    -4.911
# x     |   1E-10 m :    25.786
# y     |   1E-10 m :     9.515
# z     |   1E-10 m :     6.558
# a     |       deg :   125.841
# b     |       deg :   142.287
# g     |       deg :    41.758
# mueff |        Bm :     9.581
# shift |       ppm :     0.000
# B0    |         T :    18.790
# temp  |         K :   298.150
# t1e   |        ps :     0.189
# taur  |        ns :     4.250

import os, pymol
curdir = os.path.dirname(pymol.__script__)
set normalize_ccp4_maps, off
meshfile = os.path.join(curdir, './isomap.pml.ccp4')
cmd.load(meshfile, 'isomap', 1, 'ccp4')
isosurface pos_isomap, isomap, 1.0
isosurface neg_isomap, isomap, -1.0
set transparency, 0.5, pos_isomap
set transparency, 0.5, neg_isomap
set surface_color, blue, pos_isomap
set surface_color, red, neg_isomap
pseudoatom ori_isomap, pos=[25.786, 9.515, 6.558]
show spheres, ori_isomap
color pink, ori_isomap
cmd.load(os.path.join(curdir, '../data_files/4icbH_mut.pdb'),'4icbH_mut')
show_as cartoon, 4icbH_mut
