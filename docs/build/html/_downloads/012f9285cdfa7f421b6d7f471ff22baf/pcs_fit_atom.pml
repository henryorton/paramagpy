# PyMOL script generated by paramagpy
import os, pymol
curdir = os.path.dirname(pymol.__script__)
set normalize_ccp4_maps, off
cmd.load(os.path.join(curdir, '../data_files/5ev6AH.pdb'),'5ev6')
show_as cartoon, 5ev6
cmd.load(os.path.join(curdir, 'H07Tmmap.ccp4'), 'H07Tmmap', 1, 'ccp4')
isodot H07Tmmap+0.04, H07Tmmap, 0.04
set dot_color, teal, H07Tmmap+0.04
set transparency, 0.5, H07Tmmap+0.04
cmd.load(os.path.join(curdir, 'H08Tmmap.ccp4'), 'H08Tmmap', 1, 'ccp4')
isodot H08Tmmap+0.016, H08Tmmap, 0.016
set dot_color, magenta, H08Tmmap+0.016
set transparency, 0.5, H08Tmmap+0.016
cmd.load(os.path.join(curdir, 'H07Tbmap.ccp4'), 'H07Tbmap', 1, 'ccp4')
isodot H07Tbmap+0.04, H07Tbmap, 0.04
set dot_color, blue, H07Tbmap+0.04
set transparency, 0.5, H07Tbmap+0.04
cmd.load(os.path.join(curdir, 'H08Tbmap.ccp4'), 'H08Tbmap', 1, 'ccp4')
isodot H08Tbmap+0.02, H08Tbmap, 0.02
set dot_color, red, H08Tbmap+0.02
set transparency, 0.5, H08Tbmap+0.02
set dot_radius, 0.05
show sticks, ////28 and sc.
show sticks, ////28/CA
set bg_rgb=[1,1,1]
set mesh_width, 0.5
zoom ////28/H07


set_view (     0.505656540,   -0.827194929,   -0.245069817,    -0.741597414,   -0.561904311,    0.366465807,    -0.440846384,   -0.003562994,   -0.897575319,     0.000152570,    0.000080852,  -36.169487000,    48.539413452,   83.819839478,   42.674442291,    26.907037735,   45.422363281,  -20.000000000 )

ray 1600
png pcs_fit_atom.png