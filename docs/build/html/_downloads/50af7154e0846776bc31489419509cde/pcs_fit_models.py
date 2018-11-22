from paramagpy import protein, fit, dataparse, metal

# Load data
prot = protein.load_pdb('../data_files/2bcb.pdb')
rawData = dataparse.read_pcs('../data_files/calbindin_Er_HN_PCS.npc')
mStart = metal.Metal()
mStart.position = prot[0]['A'][56]['CA'].position


#### Ensemble average fitting ####
parsedData = prot.parse(rawData)
mGuess, _, _ = fit.svd_gridsearch_fit_metal_from_pcs(
	[mStart],[parsedData], radius=10, points=10)
mFit, calc, qfac = fit.nlr_fit_metal_from_pcs(mGuess, [parsedData])
mFit[0].save('calbindin_Er_HN_PCS_tensor_ensemble.txt')


#### Single model fitting ####
# Loop over models, fit tensor and keep one with best Q-factor
minQfacMod = 1E50
for model in prot:
	parsedDataMod = prot.parse(rawData, models=model.id)
	mFitMod, calcMod, qfacMod = fit.nlr_fit_metal_from_pcs(
		mGuess, [parsedDataMod])
	if qfacMod[0] < minQfacMod:
		minMod = model.id
		minParsedDataMod = parsedDataMod
		minmFitMod = mFitMod
		mincalcMod = calcMod
		minQfacMod = qfacMod


# #### Plot the correlation ####
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

# Unpack the experimental values
atoms, exp, err = zip(*parsedData)
expEnsemble, calcEnsemble = fit.ensemble_average(atoms, exp, calc[0])
atomsMod, expMod, errMod = zip(*minParsedDataMod)

# Plot all models
ax.plot(exp, calc[0], marker='o', lw=0, ms=2, c='g', 
	alpha=0.5, label="All models: Q = {:5.4f}".format(qfac[0]))

# Plot the ensemble average
ax.plot(expEnsemble, calcEnsemble, marker='o', lw=0, ms=2, c='r', 
	alpha=0.5, label="Ensemble Average: Q = {:5.4f}".format(qfac[0]))

# Plot the model with minimum Q-factor
ax.plot(expMod, mincalcMod[0], marker='o', lw=0, ms=2, c='b', 
	alpha=0.5, label="Best Model ({0:}): Q = {1:5.4f}".format(
		minMod, minQfacMod[0]))

# Plot a diagonal
l, h = ax.get_xlim()
ax.plot([l,h],[l,h],'-k',zorder=0)
ax.set_xlim(l,h)
ax.set_ylim(l,h)

# Make axis labels and save figure
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
fig.savefig("pcs_fit_models.png")
