# PyMOL macro for loading tensor isosurface from paramagpy
# ax    | 1E-32 m^3 :   -17.200
# rh    | 1E-32 m^3 :    -7.600
# x     |   1E-10 m :     7.072
# y     |   1E-10 m :    12.678
# z     |   1E-10 m :    33.081
# a     |       deg :   154.694
# b     |       deg :    42.478
# g     |       deg :    17.144
# mueff |        Bm :     0.000
# shift |       ppm :     0.000
# B0    |         T :    18.790
# temp  |         K :   298.150
# t1e   |        ps :     0.000
# taur  |        ns :     0.000

import os, pymol
curdir = os.path.dirname(pymol.__script__)
set normalize_ccp4_maps, off
meshfile = os.path.join(curdir, './isomap.pml.ccp4')
cmd.load(meshfile, 'isomap', 1, 'ccp4')
isosurface pos_isomap, isomap, 150
isosurface neg_isomap, isomap, -150
set transparency, 0.5, pos_isomap
set transparency, 0.5, neg_isomap
set surface_color, blue, pos_isomap
set surface_color, red, neg_isomap
pseudoatom ori_isomap, pos=[7.07200003, 12.6780005, 33.081001300000004]
show spheres, ori_isomap
color pink, ori_isomap
cmd.load(os.path.join(curdir, '../data_files/1ig5/mut_H.pdb'),'mut_H')
show_as cartoon, mut_H
