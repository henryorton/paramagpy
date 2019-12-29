from paramagpy import protein, metal, fit, dataparse
from matplotlib import pyplot as plt
import numpy as np

prot = protein.load_pdb('../data_files/parashift_Tb.pdb')
pre_exp = dataparse.read_pre('../data_files/parashift_Tb_R1_exp.pre')
exp = prot.parse(pre_exp)

m = metal.Metal(taur=0.42E-9, B0=1.0, temperature=300.0)
m.set_lanthanide('Tb')

m.g_tensor = np.array([
	[1754.0, -859.0, -207.0],
	[-859.0, 2285.0, -351.0],
	[-207.0, -351.0, -196.0]]) * 1E-60

m0 = metal.Metal(taur=0.42E-9, B0=1.0, temperature=300.0)
m0.set_lanthanide('Tb')

[mfit], [data] = fit.nlr_fit_metal_from_pre([m0], [exp], params=('t1e', 'gax', 'grh', 'a','b','g'), 
	rtypes=['r1'], usesbm=False, usegsbm=True, usedsa=True, 
	usecsa=False, progress=None)

pos = np.array([a.position for a in exp['atm']])
gam = np.array([a.gamma for a in exp['atm']])

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.plot([0,3200],[0,3200], '-k')
ax.plot(exp['exp'], mfit.fast_sbm_r1(pos, gam), marker='o', lw=0, label='iso')
ax.plot(exp['exp'], mfit.fast_g_sbm_r1(pos, gam), marker='o', lw=0, label='aniso')
ax.plot(exp['exp'], m.fast_g_sbm_r1(pos, gam), marker='o', lw=0, label='literature fit')
ax.set_xlim(0,3200)
ax.set_ylim(0,3200)
ax.set_xlabel("R1 experimental /Hz")
ax.set_ylabel("R1 calcualted /Hz")
ax.set_title("Tb parashift R1")
ax.legend()
fig.tight_layout()
fig.savefig("pre_fit_aniso_dipolar.png")