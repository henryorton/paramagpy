# PyMOL macro for loading tensor isosurface from paramagpy
# ax    | 1E-32 m^3 :    -4.920
# rh    | 1E-32 m^3 :    -2.422
# x     |   1E-10 m :    25.600
# y     |   1E-10 m :     9.629
# z     |   1E-10 m :     6.542
# a     |       deg :   130.055
# b     |       deg :   137.996
# g     |       deg :    84.398
# mueff |        Bm :     0.000
# shift |       ppm :     0.000
# B0    |         T :    18.790
# temp  |         K :   298.150
# t1e   |        ps :     0.000

set normalize_ccp4_maps, off
load ./isosurface2.pml.ccp4, isomap, 1, ccp4
isosurface pos_isosurface2, isomap, 1.0
isosurface neg_isosurface2, isomap, -1.0
set transparency, 0.5, pos_isosurface2
set transparency, 0.5, neg_isosurface2
set surface_color, blue, pos_isosurface2
set surface_color, red, neg_isosurface2
pseudoatom ori_isosurface2, pos=[25.6, 9.629000000000001, 6.542]
show spheres, ori_isosurface2
color pink, ori_isosurface2
load /home/u5376227/Dropbox/PhD/git/paramagpy/examples/ho4icb43G.pdb
show_as cartoon, ho4icb43G