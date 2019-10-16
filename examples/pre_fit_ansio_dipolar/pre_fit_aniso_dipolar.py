from paramagpy import protein, metal, fit, dataparse
import numpy as np

prot = protein.load_pdb('Tb.pdb')
# prot = protein.load_pdb('Dy.pdb')

pre_exp = dataparse.read_pre('Tb_R1_exp.pre')
pre_dip = dataparse.read_pre('Tb_R1_dip.pre')

# pre_exp = dataparse.read_pre('Dy_R1_exp.pre')
# pre_dip = dataparse.read_pre('Dy_R1_dip.pre')


exp = prot.parse(pre_exp)
dip = prot.parse(pre_dip)

m = metal.Metal(taur=0.42E-9, B0=1.0, temperature=300.0)

m.set_lanthanide('Tb')
# m.set_lanthanide('Dy')

gtensor = np.array([
	[1754.0, -859.0, -207.0],
	[-859.0, 2285.0, -351.0],
	[-207.0, -351.0, -196.0]]) * 1E-60


# gtensor = np.array([
# 	[  494.0,-1180.0, 451.0],
# 	[-1180.0, 2696.0, -30.0],
# 	[  451.0,  -30.0, 518.0]]) * 1E-60

m.g_tensor = gtensor


for atom, value, err in dip:
	r1 = m.g_sbm_r1(atom.position, atom.gamma)
	print(r1, value)


