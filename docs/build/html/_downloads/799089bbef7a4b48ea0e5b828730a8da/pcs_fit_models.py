from paramagpy import protein, fit, dataparse, metal

# Load data
prot = protein.load_pdb('../data_files/2bcb.pdb')
rawData = dataparse.read_pcs('../data_files/calbindin_Er_HN_PCS.npc')
parsedData = prot.parse(rawData)

# Set metal starting position
mStart = metal.Metal()
mStart.position = prot[0]['A'][56]['CA'].position

#### Averaged fit to all models ####
[mGuess], [data] = fit.svd_gridsearch_fit_metal_from_pcs([mStart], [parsedData], radius=10, points=10, ensembleAverage=False)
[mFit], [data] = fit.nlr_fit_metal_from_pcs([mGuess], [parsedData], ensembleAverage=False)
qfac = fit.qfactor(data, ensembleAverage=False)
avg = qfac, data, mFit

#### Ensembled averaged fit to all models ####
[mGuess], [data] = fit.svd_gridsearch_fit_metal_from_pcs([mStart], [parsedData], radius=10, points=10, ensembleAverage=True)
[mFit], [data] = fit.nlr_fit_metal_from_pcs([mGuess], [parsedData], ensembleAverage=True)
qfac = fit.qfactor(data, ensembleAverage=True)
e_avg = qfac, data, mFit

#### Seperate fit for each model ####
sep = {}
for model in prot:
    singleModelData = parsedData[parsedData['mdl']==model.id]
    [mGuess], [data] = fit.svd_gridsearch_fit_metal_from_pcs([mStart], [singleModelData], radius=10, points=10)
    [mFit], [data] = fit.nlr_fit_metal_from_pcs([mGuess], [singleModelData])
    qfac = fit.qfactor(data)
    sep[model.id] = qfac, data, mFit


#### Plot the correlation ####
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

# Plot averaged fit correlation
qfac, data, mFit = avg
mFit.save('calbindin_Er_HN_PCS_tensor_average.txt')
ax.plot(data['exp'], data['cal'], marker='o', lw=0, ms=2, c='g', 
    alpha=0.5, label="Averaged Fit: Q = {:5.4f}".format(qfac))

# Plot ensemble averaged fit correlation
qfac, data, mFit = e_avg
mFit.save('calbindin_Er_HN_PCS_tensor_ensemble_average.txt')
# Ensemble average the data to get a single point for each model
data = fit.ensemble_average(data)
ax.plot(data['exp'], data['cal'], marker='o', lw=0, ms=2, c='r', 
    alpha=0.5, label="Ensemble Avg. Fit: Q = {:5.4f}".format(qfac))

# Plot best fit model correlation
# Sort fits by Qfactor and take smallest
model, (qfac, data, mFit) = sorted(sep.items(), key=lambda x: x[1][0])[0]
mFit.save('calbindin_Er_HN_PCS_tensor_best_model.txt')
ax.plot(data['exp'], data['cal'], marker='o', lw=0, ms=2, c='b', 
    alpha=0.5, label="Best Fit Model {}: Q = {:5.4f}".format(model, qfac))

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

