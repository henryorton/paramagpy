import sys
sys.path.append('../..')
from paramagpy import protein, metal, fit, dataparse
from matplotlib import pyplot as plt
import numpy as np

prot = protein.load_pdb('Tb.pdb')

pre_exp = dataparse.read_pre('Tb_R1_exp.pre')
pre_dip = dataparse.read_pre('Tb_R1_dip.pre')

exp = prot.parse(pre_exp)
dip = prot.parse(pre_dip)

m = metal.Metal(taur=0.42E-9, B0=1.0, temperature=300.0)
m.set_lanthanide('Tb')

m.g_tensor = np.array([
	[1754.0, -859.0, -207.0],
	[-859.0, 2285.0, -351.0],
	[-207.0, -351.0, -196.0]]) * 1E-60

m0 = metal.Metal(taur=0.42E-9, B0=1.0, temperature=300.0)
m0.set_lanthanide('Tb')

[mfit], [calc], [q] = fit.nlr_fit_metal_from_pre([m0], [exp], params=('t1e', 'gax', 'grh', 'a','b','g'), 
	sumIndices=None, rtypes=['r1'], usesbm=False, usegsbm=True, usedsa=True, 
	usecsa=False, progress=None)

atoms, r1, err = zip(*exp)
pos = np.array([a.position for a in atoms])
gam = np.array([a.gamma for a in atoms])

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.plot([0,3200],[0,3200], '-k')
ax.plot(r1, mfit.fast_sbm_r1(pos, gam), marker='o', lw=0, label='iso')
ax.plot(r1, mfit.fast_g_sbm_r1(pos, gam), marker='o', lw=0, label='aniso')
ax.plot(r1, m.fast_g_sbm_r1(pos, gam), marker='o', lw=0, label='literature fit')
ax.set_xlim(0,3200)
ax.set_ylim(0,3200)
ax.set_xlabel("R1 experimental /Hz")
ax.set_ylabel("R1 calcualted /Hz")
ax.set_title("Tb parashift R1")
ax.legend()
fig.savefig("pre_fit_aniso_diplar.png", dpi=200)