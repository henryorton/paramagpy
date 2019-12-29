from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# Load the PCS data
rawData = dataparse.read_pcs('../data_files/calbindin_Er_HN_PCS_errors.npc')

# Associate PCS data with atoms of the PDB
parsedData = prot.parse(rawData)

# Define an initial tensor
mStart = metal.Metal()

# Set the starting position to an atom close to the metal
mStart.position = prot[0]['A'][56]['CA'].position

# Calculate an initial tensor from an SVD gridsearch
[mGuess], [data] = fit.svd_gridsearch_fit_metal_from_pcs(
	[mStart],[parsedData], radius=10, points=10)

# Refine the tensor using non-linear regression
[mFit], [data] = fit.nlr_fit_metal_from_pcs([mGuess], [parsedData])

qfac = fit.qfactor(data)

# Save the fitted tensor to file
mFit.save('calbindin_Er_HN_PCS_tensor_errors.txt')

#### Plot the correlation ####
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

# Plot the data
ax.errorbar(data['exp'], data['cal'], xerr=data['err'], fmt='o', c='r', ms=2, 
	ecolor='k', capsize=3, label="Q-factor = {:5.4f}".format(qfac))

# Plot a diagonal
l, h = ax.get_xlim()
ax.plot([l,h],[l,h],'grey',zorder=0)
ax.set_xlim(l,h)
ax.set_ylim(l,h)

# Make axis labels and save figure
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
fig.savefig("pcs_fit_error.png")