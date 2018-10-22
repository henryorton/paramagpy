# PyMOL macro for loading tensor isosurface from paramagpy
# ax    | 1E-32 m^3 :     0.000
# rh    | 1E-32 m^3 :     0.000
# x     |   1E-10 m :     0.000
# y     |   1E-10 m :     0.000
# z     |   1E-10 m :     0.000
# a     |       deg :     0.000
# b     |       deg :     0.000
# g     |       deg :     0.000
# mueff |        Bm :     0.000
# shift |       ppm :     0.000
# B0    |         T :    18.790
# temp  |         K :   298.150
# t1e   |        ps :     0.000

set normalize_ccp4_maps, off
load ./buenosdias.pml.ccp4, isomap, 1, ccp4
isosurface pos_buenosdias, isomap, 1.0
isosurface neg_buenosdias, isomap, -1.0
set transparency, 0.5, pos_buenosdias
set transparency, 0.5, neg_buenosdias
set surface_color, blue, pos_buenosdias
 #<<<<Change colour hereset surface_color, red, neg_buenosdias
  #<<<<Change colour herepseudoatom ori_buenosdias, pos=[0.0, 0.0, 0.0]
show spheres, ori_buenosdias
color pink, ori_buenosdias
