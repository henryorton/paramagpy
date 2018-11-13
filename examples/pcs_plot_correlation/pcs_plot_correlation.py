# Make two list of experimental and calculated values
xs = []
ys = []
for atom, exp_pcs, error in parsedData:
	xs.append(exp_pcs)
	ys.append(mFit[0].atom_pcs(atom))

# Plot the correlation
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

ax.plot([min(xs+ys),max(xs+ys)],[min(xs+ys),max(xs+ys)], '-k')
ax.plot(xs, ys, marker='o', lw=0, ms=3, c='r')
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
fig.savefig("pcs_plot.png")

# Plot the isosurface to be opened in PyMol
print(prot.id)
mFit[0].isomap(prot.id, density=1)