from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
prot = protein.load_pdb('2bcb.pdb')

# Load the PCS data
rawData = dataparse.read_pcs('calbindin_Er_HN_PCS.npc')

# Associate PCS data with atoms of the PDB
parsedData = prot.parse(rawData)

# Define an initial tensor, and set the lanthanide to Er
mStart = metal.Metal()
# mStart.set_lanthanide('Er')

# Set the starting position to an atom close to the metal
mStart.position = prot[0]['A'][56]['CA'].position

# Calculate an initial tensor from an SVD gridsearch
mGuess = fit.svd_gridsearch_fit_metal_from_pcs(
	[mStart],[parsedData], radius=10, points=10)

print(mGuess[0].info())

# Refine the tensor using non-linear regression
fitParameters = ['x','y','z','ax','rh','a','b','g']
mFit = fit.nlr_fit_metal_from_pcs(mGuess, [parsedData], fitParameters)

# Print the info of the fitted tensor
print(mFit[0].info())




idxs = []
exps = []
calcs = []

qs = []
for m in prot:
	exp = []
	calc = []

	for a, e, r in parsedData:
		if a.parent.parent.parent.id==m.id:
			idxs.append(a.serial_number)
			exps.append(e)
			calcs.append(mFit[0].atom_pcs(a))
			exp.append(e)
			calc.append(mFit[0].atom_pcs(a))
	# sumIndices = fit.clean_indices(idxs)
	q = fit.qfactor(exp, calc)
	print(q)
	qs.append(q)

import numpy as np

sumIndices = fit.clean_indices(idxs)

print(np.mean(qs))
print(fit.qfactor(exps, calcs, sumIndices))


# from matplotlib import pyplot as plt

# xs = []
# ys = []
# for atom, exp_pcs, err in parsedData:
# 	xs.append(exp_pcs)
# 	ys.append(mFit[0].atom_pcs(atom))

# plt.scatter(xs, ys)
# plt.show()