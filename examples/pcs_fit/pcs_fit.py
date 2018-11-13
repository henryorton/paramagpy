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
mGuess = fit.svd_gridsearch_fit_metal_from_pcs(
	[mStart],[parsedData], radius=10, points=10)

# Refine the tensor using non-linear regression
fitParameters = ['x','y','z','ax','rh','a','b','g']
mFit = fit.nlr_fit_metal_from_pcs(mGuess, [parsedData], fitParameters)

# Print the info of the fitted tensor
print(mFit[0].info())

