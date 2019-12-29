from paramagpy import protein, fit, dataparse, metal

# Load data
prot = protein.load_pdb('../data_files/2bcb.pdb')
rawData = dataparse.read_pcs('../data_files/calbindin_Er_HN_PCS.npc')
mStart = metal.Metal()
mStart.position = prot[0]['A'][56]['CA'].position


#### Ensemble average fitting ####
parsedData = prot.parse(rawData)
[mGuess], [data] = fit.svd_gridsearch_fit_metal_from_pcs(
	[mStart],[parsedData], radius=10, points=10, ensembleAverage=True)
[mFit], [data] = fit.nlr_fit_metal_from_pcs([mGuess], [parsedData], ensembleAverage=True)

qfac = fit.qfactor(data, ensembleAverage=True)
dataEAv = fit.ensemble_average(data)

mFit.save('calbindin_Er_HN_PCS_tensor_ensemble.txt')


#### Single model fitting ####
[mFitMod], [dataMod] = fit.nlr_fit_metal_from_pcs(
	[mGuess], [parsedData])
qs = {}
for mdl in set(dataMod['mdl']):
	qs[mdl] = fit.qfactor(dataMod[dataMod['mdl']==mdl])

minModel, minQfac = sorted(qs.items(), key=lambda x: x[1])[0]
minData = dataMod[dataMod['mdl']==minModel]


# #### Plot the correlation ####
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

# Plot all models
ax.plot(data['exp'], data['cal'], marker='o', lw=0, ms=2, c='g', 
	alpha=0.5, label="All models: Q = {:5.4f}".format(qfac))

# Plot the ensemble average
ax.plot(dataEAv['exp'], dataEAv['cal'], marker='o', lw=0, ms=2, c='r', 
	alpha=0.5, label="Ensemble Average: Q = {:5.4f}".format(qfac))

# Plot the model with minimum Q-factor
ax.plot(minData['exp'], minData['cal'], marker='o', lw=0, ms=2, c='b', 
	alpha=0.5, label="Best Model ({0:}): Q = {1:5.4f}".format(
		minModel, minQfac))

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
