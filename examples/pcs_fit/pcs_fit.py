from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# Load the PCS data
rawData = dataparse.read_pcs('../data_files/calbindin_Er_HN_PCS.npc')

# Associate PCS data with atoms of the PDB
parsedData = prot.parse(rawData)

# Define an initial tensor
mStart = metal.Metal()

# Set the starting position to an atom close to the metal
mStart.position = prot[0]['A'][56]['CA'].position

# Calculate an initial tensor from an SVD gridsearch
mGuess, calc, qfac = fit.svd_gridsearch_fit_metal_from_pcs(
    [mStart], [parsedData], radius=10, points=10)

# Refine the tensor using non-linear regression
mFit, calc, qfac = fit.nlr_fit_metal_from_pcs(mGuess, [parsedData])

# Save the fitted tensor to file
mFit[0].save('calbindin_Er_HN_PCS_tensor.txt')

#### Plot the correlation ####
from matplotlib import pyplot as plt

fig, ax = plt.subplots(figsize=(5, 5))

# Unpack the experimental values
atoms, experiment, errors = zip(*parsedData)

# Plot the data
ax.plot(experiment, calc[0], marker='o', lw=0, ms=3, c='r',
        label="Q-factor = {:5.4f}".format(qfac[0]))

# Plot a diagonal
l, h = ax.get_xlim()
ax.plot([l, h], [l, h], '-k', zorder=0)
ax.set_xlim(l, h)
ax.set_ylim(l, h)

# Make axis labels and save figure
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
fig.savefig("pcs_fit.png")
