from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# Load the PCS data
rawData = dataparse.read_pcs('../data_files/calbindin_Er_HN_PCS.npc')

# Associate PCS data with atoms of the PDB
parsedData = prot.parse(rawData)

# Load the fitted tensor
met = metal.load_tensor('../data_files/calbindin_Er_HN_PCS_tensor.txt')

# Make two lists of experimental and calculated values
exp = []
cal = []
for atom, exp_pcs, error in parsedData:
	exp.append(exp_pcs)
	cal.append(met.atom_pcs(atom))

# Calculate the Q-factor for the correlation
qfactor = fit.qfactor(exp, cal)

#### Plot the correlation ####
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

# Plot the data
ax.plot(exp, cal, marker='o', lw=0, ms=3, c='r', 
	label="Q-factor = {:5.4f}".format(qfactor))

# Plot a diagonal
l, h = ax.get_xlim()
ax.plot([l,h],[l,h],'-k',zorder=0)
ax.set_xlim(l,h)
ax.set_ylim(l,h)

# Axis labels
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
fig.savefig("pcs_plot_correlation.png")
