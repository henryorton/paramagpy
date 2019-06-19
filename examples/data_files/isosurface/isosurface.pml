# PyMOL macro for loading tensor isosurface from paramagpy
# ax    | 1E-32 m^3 :    31.882
# rh    | 1E-32 m^3 :    12.882
# x     |   1E-10 m :    26.133
# y     |   1E-10 m :     9.305
# z     |   1E-10 m :     6.222
# a     |       deg :   150.083
# b     |       deg :   151.389
# g     |       deg :    70.055
# mueff |        Bm :     0.000
# shift |       ppm :     0.000
# B0    |         T :    18.790
# temp  |         K :   298.150
# t1e   |        ps :     0.000
# taur  |        ns :     0.000

import os, pymol
curdir = os.path.dirname(pymol.__script__)
set normalize_ccp4_maps, off
meshfile = os.path.join(curdir, 'isosurface.pml.ccp4')
cmd.load(meshfile, 'isomap', 1, 'ccp4')
isosurface pos_isosurface, isomap, 1.0
isosurface neg_isosurface, isomap, -1.0
set transparency, 0.5, pos_isosurface
set transparency, 0.5, neg_isosurface
set surface_color, blue, pos_isosurface
set surface_color, red, neg_isosurface
pseudoatom ori_isosurface, pos=[26.133, 9.305, 6.222]
show spheres, ori_isosurface
color pink, ori_isosurface
cmd.load(os.path.join(curdir, '4icbH_mut.pdb'),'4icbH_mut')
show_as cartoon, 4icbH_mut
