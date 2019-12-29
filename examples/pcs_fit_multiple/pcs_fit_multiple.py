from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# Load the PCS data
rawData1 = dataparse.read_pcs('../data_files/calbindin_Tb_HN_PCS.npc')
rawData2 = dataparse.read_pcs('../data_files/calbindin_Er_HN_PCS.npc')
rawData3 = dataparse.read_pcs('../data_files/calbindin_Yb_HN_PCS.npc')

# Associate PCS data with atoms of the PDB
parsedData = []
for rd in [rawData1, rawData2, rawData3]:
	parsedData.append(prot.parse(rd))

# Make a list of starting tensors
mStart = [metal.Metal(), metal.Metal(), metal.Metal()]

# Set the starting position to an atom close to the metal
mStart[0].position = prot[0]['A'][56]['CA'].position

# Calculate initial tensors from an SVD gridsearch
mGuess, datas = fit.svd_gridsearch_fit_metal_from_pcs(
	mStart, parsedData, radius=10, points=10)

# Refine the tensors using non-linear regression
fitParameters = ['x','y','z','ax','rh','a','b','g']
mFit, datas = fit.nlr_fit_metal_from_pcs(mGuess, parsedData, fitParameters)

# Save the fitted tensors to files
for name, metal in zip(['Tb','Er','Yb'], mFit):
	metal.save("tensor_{}.txt".format(name))

#### Plot the correlation ####
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

# Plot the data
for d, name, colour in zip(datas, ['Tb','Er','Yb'],['r','g','b']):
	qfactor = fit.qfactor(d)
	ax.plot(d['exp'], d['cal'], marker='o', lw=0, ms=1, c=colour,
		label="{0:} - {1:5.3f}".format(name, qfactor))

# Plot a diagonal
l, h = ax.get_xlim()
ax.plot([l,h],[l,h],'-k',zorder=0)
ax.set_xlim(l,h)
ax.set_ylim(l,h)

# Axis labels
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
fig.savefig("pcs_fit_multiple.png")