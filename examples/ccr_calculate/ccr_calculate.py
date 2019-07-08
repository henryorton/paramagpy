from paramagpy import protein, fit, dataparse, metal
import numpy as np

# Load the PDB file and get iron centre
prot = protein.load_pdb('../data_files/1bzrH.pdb')
ironAtom = prot[0]['A'][("H_HEM",154," ")]['FE']

# Make high and low spin Fe paramagnetic centers
met_cn = metal.Metal(position=ironAtom.position, 
					 B0=18.79, 
					 temperature=303.0,
					 taur=5.7E-9)
met_f = met_cn.copy()
met_cn.iso = 4.4E-32
met_f.iso = 30.1E-32

# Load experimental data
data_cn = prot.parse(dataparse.read_ccr("../data_files/myoglobin_cn.ccr"))
data_f = prot.parse(dataparse.read_ccr("../data_files/myoglobin_f.ccr"))

# Calculate the cross-correlated realxation
compare_cn = []
for H, N, value, error in data_cn:
	delta = met_cn.atom_ccr(H, N)
	compare_cn.append((value, delta))

compare_f = []
for H, N, value, error in data_f:
	delta = met_f.atom_ccr(H, N)
	compare_f.append((value, delta))

#### Plot the correlation ####
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

# Plot the data correlations
ax.scatter(*zip(*compare_cn), s=7, label="myo_cn")
ax.scatter(*zip(*compare_f), s=7, label="myo_f")

# Plot a diagonal
l, h = ax.get_xlim()
ax.plot([l,h],[l,h],'-k',zorder=0)
ax.set_xlim(l,h)
ax.set_ylim(l,h)

# Make axis labels and save figure
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
fig.savefig("ccr_calculate.png")
plt.show()


H = prot[0]['A'][23]['H']
N = prot[0]['A'][23]['N']

print(met_f.dipole_shift_tensor(H.position)*1E6)
print(N.dipole_shift_tensor(H.position)*1E6)

print(met_f.atom_ccr(H, N))


# print(met_f.dsa_r2(H.position, H.gamma))

# {{-2.44981292,  2.46216645, -4.32282316},
#  { 2.46216645, -1.73676364, -4.99385716},
#  {-4.32282316, -4.99385716,  4.18657656}}

# {{-231.53542458,    7.00415538,  130.74435749},
#  {   7.00415538,  138.77405716,   -2.47205773},
#  { 130.74435749,   -2.47205773,   92.76136742}}


# {{-214.72723312,  119.79115598,   89.36823381},
#  { 119.79115598,   98.36758395,  -30.26987692},
#  {  89.36823381,  -30.26987692,  116.35964917}}