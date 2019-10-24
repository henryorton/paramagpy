import sys
sys.path.append('../..')
from paramagpy import protein, metal, fit, dataparse
import numpy as np
import time

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

posarr = []
gamarr = []
values = []
for atom, value, err in dip:
	posarr.append(atom.position)
	gamarr.append(atom.gamma)
	values.append(value)

posarr = np.array(posarr)
gamarr = np.array(gamarr)
values = np.array(values)

l = len(dip)

t1 = time.time()
for atom, value, err in dip[:l]:
	r1 = m.g_sbm_r1(atom.position, atom.gamma)

t2 = time.time()
r1 = m.fast_g_sbm_r1(posarr[:l], gamarr[:l])

t3 = time.time()

print(t2-t1)
print(t3-t2)
