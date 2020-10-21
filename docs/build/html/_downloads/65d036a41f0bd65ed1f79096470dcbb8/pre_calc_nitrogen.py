from paramagpy import protein, metal, dataparse

# Load the PDB file
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# Load PRE data
rawData = dataparse.read_pre('../data_files/calbindin_Tb_N_R1_600.pre')

# Parse PRE data
data = prot.parse(rawData)

# Load the fitted tensor and set relevant parameters
met = metal.load_tensor('../data_files/calbindin_Tb_HN_PCS_tensor.txt')
met.B0 = 14.1
met.T = 298.0
met.taur = 4.25E-9

# Loop over nitrogen atoms and calculate PRE with and without CSA
exp = []
cal = []
cal_csa = []
for atom, pre, err in data[['atm','exp','err']]:
	exp.append(pre)
	cal.append(met.atom_pre(atom, rtype='r1'))
	cal_csa.append(met.atom_pre(atom, rtype='r1', csa=atom.csa))

#### Plot the correlation ####
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

# Plot the data
ax.scatter(exp, cal, label="Standard Theory")
ax.scatter(exp, cal_csa, label="CSA x Curie spin")

# Plot a diagonal
l, h = ax.get_xlim()
ax.plot([l,h],[l,h],'grey',zorder=0)
ax.set_xlim(l,h)
ax.set_ylim(l,h)

# Make axis labels and save figure
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
fig.savefig("pre_calc_nitrogen.png")