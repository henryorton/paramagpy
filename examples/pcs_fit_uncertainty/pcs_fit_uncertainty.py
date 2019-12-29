from paramagpy import protein, fit, dataparse, metal
import numpy as np

# Load the PDB file
prot = protein.load_pdb('../data_files/2bcb.pdb')

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

# Estimate uncertainty sourcing noise from the models of the PDB
[mod_all], [mod_std] = fit.fit_error_models(fit.nlr_fit_metal_from_pcs, 
	initMetals=[mFit], dataArrays=[parsedData])

mod_std.save('error_tensor_models.txt')

# Estimate uncertainty sourcing noise from experimental uncertainties
[mc_all], [mc_std] = fit.fit_error_monte_carlo(fit.nlr_fit_metal_from_pcs, 
	50, initMetals=[mFit], dataArrays=[parsedData])

mod_std.save('error_tensor_monte_carlo.txt')

# Estimate uncertainty sourcing noise from sample fractions
[bs_all], [bs_std] = fit.fit_error_bootstrap(fit.nlr_fit_metal_from_pcs, 
	50, 0.8, initMetals=[mFit], dataArrays=[parsedData])

mod_std.save('error_tensor_bootstrap.txt')

#### Plot Sanson-Flamsteed ####
from matplotlib import pyplot as plt

def transform(vector):
	x, y, z = vector
	theta = np.arctan2(y, x)
	phi = -np.arccos(z) + np.pi/2.
	return theta, phi

for name, mset in [('models',mod_all), ('monte_carlo',mc_all), ('bootstrap',bs_all)]:
	spcoords = []
	for m in mset:
		x, y, z = m.rotationMatrix.T
		spcoords.append(tuple(map(transform, [x,y,z])))
	points = zip(*spcoords)
	fig = plt.figure(figsize=(5, 3), dpi=100)
	ax = fig.add_subplot(111, projection='hammer')
	ax.set_xlabel("theta")
	ax.set_ylabel("phi")
	ax.set_title(name)
	ax.grid()
	for data, col, label in zip(points, ['r','g','b'], ['x','y','z']):
		theta, phi = zip(*data)
		ax.scatter(theta, phi, s=0.4, c=col, label=label, zorder=10)
	ax.legend()
	fig.savefig("{}.png".format(name))
