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
mGuess, calc, qfac = fit.svd_gridsearch_fit_metal_from_pcs(
	[mStart],[parsedData], radius=10, points=10)

# Refine the tensor using non-linear regression
mFit, calc, qfac = fit.nlr_fit_metal_from_pcs(mGuess, [parsedData])

# mets, stdm = fit.pcs_fit_error_bootstrap(mFit, [parsedData], 10, 0.95)

mets, stdm = fit.pcs_fit_error_monte_carlo(mFit, [parsedData], 50)



# self.errorTensor.set_params(devs.items())
# def transform(vector):
# 	x, y, z = vector
# 	theta = np.arctan2(y, x)
# 	phi = -np.arccos(z) + np.pi/2.
# 	return theta, phi

# spcoords = []
# for eulers in zip(stds['a'], stds['b'], stds['g']):
# 	rotationMatrix = metal.euler_to_matrix(np.array(eulers))
# 	x, y, z = rotationMatrix.T
# 	spcoords.append(tuple(map(transform, [x,y,z])))

# self.type_tensors()
# self.plot(zip(*spcoords))

# def plot(self, points=None):
# self.axes.clear()
# self.axes.set_xlabel("theta")
# self.axes.set_ylabel("phi")
# self.axes.grid()

# if points is not None:
# 	for data, col, label in zip(points, ['r','g','b'], ['x','y','z']):
# 		theta, phi = zip(*data)
# 		self.axes.scatter(theta, phi, s=0.4, c=col, label=label)
# 		# self.axes.plot(theta, phi, marker='x',
# 			# lw=0, ms=3, label=label, color=col)
# 	self.axes.legend()
# self.canvas.draw()


















# # Save the fitted tensor to file
# mFit[0].save('calbindin_Er_HN_PCS_tensor_errors.txt')

# #### Plot the correlation ####
# from matplotlib import pyplot as plt
# fig, ax = plt.subplots(figsize=(5,5))

# # Unpack the experimental values
# atoms, experiment, errors = zip(*parsedData)

# # Plot the data
# ax.errorbar(experiment, calc[0], xerr=errors, fmt='o', c='r', ms=2, 
# 	ecolor='k', capsize=3, label="Q-factor = {:5.4f}".format(qfac[0]))

# # Plot a diagonal
# l, h = ax.get_xlim()
# ax.plot([l,h],[l,h],'grey',zorder=0)
# ax.set_xlim(l,h)
# ax.set_ylim(l,h)

# # Make axis labels and save figure
# ax.set_xlabel("Experiment")
# ax.set_ylabel("Calculated")
# ax.legend()
# fig.savefig("pcs_fit_error.png")