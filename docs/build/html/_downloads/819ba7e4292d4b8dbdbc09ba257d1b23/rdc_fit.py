from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/2kox.pdb')

# Load the RDC data
rawData1 = dataparse.read_rdc('../data_files/ubiquitin_a28c_c1_Tb_HN.rdc')
rawData2 = dataparse.read_rdc('../data_files/ubiquitin_s57c_c1_Tb_HN.rdc')

# Associate RDC data with atoms of the PDB
parsedData1 = prot.parse(rawData1)
parsedData2 = prot.parse(rawData2)

# Define an initial tensor
mStart1 = metal.Metal(B0=18.8, temperature=308.0)
mStart2 = metal.Metal(B0=18.8, temperature=308.0)

# Calculate the tensor using SVD
sol1 = fit.svd_fit_metal_from_rdc(mStart1, parsedData1)
sol2 = fit.svd_fit_metal_from_rdc(mStart2, parsedData2)

# Save the fitted tensor to file
sol1[0].save('ubiquitin_a28c_c1_Tb_tensor.txt')
sol2[0].save('ubiquitin_s57c_c1_Tb_tensor.txt')

#### Plot the correlation ####
from matplotlib import pyplot as plt
fig = plt.figure(figsize=(5,10))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.set_title('A28C-C1-Tb')
ax2.set_title('S57C-C1-Tb')

for sol, ax, pdata in zip([sol1,sol2], [ax1,ax2], 
	[parsedData1,parsedData2]):

	# Unpack the experimental values
	atoms1, atoms2, exp, err = zip(*pdata)
	atoms = list(zip(atoms1, atoms2))
	metal, calc, qfac = sol
	expEnsemble, calcEnsemble = fit.ensemble_average(atoms, exp, calc)

	# Plot all models
	ax.plot(exp, calc, marker='o', lw=0, ms=2, c='b', 
		alpha=0.5, label="All models: Q = {:5.4f}".format(qfac))

	# Plot the ensemble average
	ax.plot(expEnsemble, calcEnsemble, marker='o', lw=0, ms=2, c='r', 
		label="Ensemble Average: Q = {:5.4f}".format(qfac))

	# Plot a diagonal
	l, h = ax.get_xlim()
	ax.plot([l,h],[l,h],'-k',zorder=0)
	ax.set_xlim(l,h)
	ax.set_ylim(l,h)

	# Make axis labels and save figure
	ax.set_xlabel("Experiment")
	ax.set_ylabel("Calculated")
	ax.legend()

fig.tight_layout()
fig.savefig("rdc_fit.png")
