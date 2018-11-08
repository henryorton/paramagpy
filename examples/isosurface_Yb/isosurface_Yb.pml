# PyMOL macro for loading tensor isosurface from paramagpy
# ax    | 1E-32 m^3 :    -5.184
# rh    | 1E-32 m^3 :    -2.427
# x     |   1E-10 m :    25.720
# y     |   1E-10 m :     9.498
# z     |   1E-10 m :     6.633
# a     |       deg :   128.884
# b     |       deg :   138.194
# g     |       deg :    86.397
# mueff |        Bm :     0.000
# shift |       ppm :     0.000
# B0    |         T :    18.790
# temp  |         K :   298.150
# t1e   |        ps :     0.000

import os, pymol
curdir = os.path.dirname(pymol.__script__)
set normalize_ccp4_maps, off
meshfile = os.path.join(curdir, 'isosurface_Yb.pml.ccp4')
cmd.load(meshfile, 'isomap', 1, 'ccp4')
isosurface pos_isosurface_Yb, isomap, 1.0
isosurface neg_isosurface_Yb, isomap, -1.0
set transparency, 0.5, pos_isosurface_Yb
set transparency, 0.5, neg_isosurface_Yb
set surface_color, blue, pos_isosurface_Yb
set surface_color, red, neg_isosurface_Yb
pseudoatom ori_isosurface_Yb, pos=[25.72, 9.498, 6.633]
show spheres, ori_isosurface_Yb
color pink, ori_isosurface_Yb
cmd.load(os.path.join(curdir, 'ho4icb43G.pdb'))
show_as cartoon, ho4icb43G
