#!/usr/bin/python3


PDB_FILE = "/home/u5376227/4icbH_mut.pdb"
START_ATOM = 56, 'CA'
SVD_RADIUS = 10.0
SVD_DENSITY = 1.0
USE_RACS = False
USE_RADS = False
FITTED_TENSOR_FILE_NAME = "pcs_tensor.txt"
BACK_CALCULATED_PCS_FILE_NAME = 'pcscal.npc'
CORRELATION_PLOT = True
CORRELATION_PLOT_FILE_NAME = "pcs_fit.png"
CORRELATION_PLOT_ANNOTATE = False
PYMOL_ISOSURFACE = True



from paramagpy import protein, fit, dataparse, metal
import sys, subprocess

expPCS_fileName = sys.argv[1]

# Load the PDB file
prot = protein.load_pdb(PDB_FILE)

# Load the PCS data
rawData = dataparse.read_pcs(expPCS_fileName)

# Associate PCS data with atoms of the PDB
parsedData = prot.parse(rawData)

# Define an initial tensor
mStart = metal.Metal()

# Set the starting position to an atom close to the metal
mStart.position = prot[0]['A'][START_ATOM[0]][START_ATOM[1]].position

# Calculate an initial tensor from an SVD gridsearch
mGuess, calc, qfac = fit.svd_gridsearch_fit_metal_from_pcs(
	[mStart],[parsedData], radius=SVD_RADIUS, points=int(SVD_RADIUS/SVD_DENSITY))

# Refine the tensor using non-linear regression
mFit, calc, qfac = fit.nlr_fit_metal_from_pcs(mGuess, [parsedData], 
	useracs=USE_RACS, userads=USE_RADS)

# Save the fitted tensor to file
mFit[0].save(FITTED_TENSOR_FILE_NAME)

# Save calculated PCS values
back_calc = []
for atom in prot.get_atoms():
	value = mFit[0].atom_pcs(atom, racs=USE_RACS, rads=USE_RADS)
	_, mdl, chn, (_, seq, _), (atm, _) = atom.get_full_id()
	back_calc.append((seq, atm, value))

with open(BACK_CALCULATED_PCS_FILE_NAME, 'w') as o:
	for seq, atm, value in back_calc:
		line = "{0:<3d} {1:3s} {2:8.3f}  0.0\n".format(seq, atm, value)
		o.write(line)


if PYMOL_ISOSURFACE:
	mFit[0].isomap(prot.id, density=1, isoval=1.0)
	subprocess.Popen(["pymol", "isomap.pml"])


if CORRELATION_PLOT:
	#### Plot the correlation ####
	from matplotlib import pyplot as plt
	fig, ax = plt.subplots(figsize=(5,5))

	# Unpack the experimental values
	atoms, experiment, errors = zip(*parsedData)

	# Plot the data
	ax.plot(experiment, calc[0], marker='o', lw=0, ms=3, c='r', 
		label="{0:3d} total PCS\nQ-factor = {1:5.4f}".format(len(atoms), qfac[0]))

	# Plot annotations
	if CORRELATION_PLOT_ANNOTATE:
		for atom, exp, cal in zip(atoms, experiment, calc[0]):
			_, mdl, chn, (_, seq, _), (atm, _) = atom.get_full_id()
			s = "{}{}".format(atm,seq)
			ax.annotate(s, xy=(exp, cal), fontsize=8)

	# Plot a diagonal
	l, h = ax.get_xlim()
	ax.plot([l,h],[l,h],'-k',zorder=0)
	ax.set_xlim(l,h)
	ax.set_ylim(l,h)

	# Make axis labels and save figure
	ax.set_xlabel("Experiment")
	ax.set_ylabel("Calculated")
	ax.legend()
	fig.savefig(CORRELATION_PLOT_FILE_NAME)
	plt.show()
