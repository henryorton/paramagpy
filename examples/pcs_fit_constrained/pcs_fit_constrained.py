from paramagpy import protein, fit, dataparse, metal

# Load data
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')
rawData = dataparse.read_pcs('../data_files/calbindin_Er_HN_PCS.npc')
parsedData = prot.parse(rawData)
mStart = metal.Metal()

# Set the starting position to Calcium ion heteroatom in PDB
mStart.position = prot[0]['A'][('H_ CA', 77, ' ')]['CA'].position

# Calculate tensor by SVD
mFit, calc, qfac = fit.svd_gridsearch_fit_metal_from_pcs(
    [mStart], [parsedData], radius=0, points=1)

mFit[0].save('calbindin_Er_HN_PCS_tensor_position_constrained.txt')

# Calculate axially symmetric tensor by NRL
mFitAx, calcAx, qfacAx = fit.nlr_fit_metal_from_pcs(
    [mStart], [parsedData], params=('ax', 'b', 'g', 'x', 'y', 'z'))

mFitAx[0].save('calbindin_Er_HN_PCS_tensor_axially_symmetric.txt')

#### Plot the correlation ####
from matplotlib import pyplot as plt

fig, ax = plt.subplots(figsize=(5, 5))

# Unpack the experimental values
atoms, experiment, errors = zip(*parsedData)

# Plot the data
ax.plot(experiment, calc[0], marker='o', lw=0, ms=2, c='r',
        label="Position constrained: Q = {:5.4f}".format(qfac[0]))

ax.plot(experiment, calcAx[0], marker='o', lw=0, ms=2, c='b',
        label="Axially symmetric: Q = {:5.4f}".format(qfacAx[0]))

# Plot a diagonal
l, h = ax.get_xlim()
ax.plot([l, h], [l, h], '-k', zorder=0)
ax.set_xlim(l, h)
ax.set_ylim(l, h)

# Make axis labels and save figure
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
fig.savefig("pcs_fit_constrained.png")
