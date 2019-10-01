from paramagpy import protein, dataparse, metal

# Load the PDB file and get iron centre
prot = protein.load_pdb('../data_files/1bzrH.pdb')
ironAtom = prot[0]['A'][("H_HEM", 154, " ")]['FE']

# Make high and low spin Fe paramagnetic centers
met_cn = metal.Metal(position=ironAtom.position,
                     B0=18.79,
                     temperature=303.0,
                     taur=5.7E-9)
met_f = met_cn.copy()
met_cn.iso = 4.4E-32
met_f.iso = 30.1E-32

# Load experimental data
data_cn = prot.parse(dataparse.read_ccr("../data_files/myoglobin_cn.ccr"))
data_f = prot.parse(dataparse.read_ccr("../data_files/myoglobin_f.ccr"))

# Calculate the cross-correlated realxation
compare_cn = []
for H, N, value, error in data_cn:
    delta = met_cn.atom_ccr(H, N)
    compare_cn.append((value, delta * 0.5))

compare_f = []
for H, N, value, error in data_f:
    delta = met_f.atom_ccr(H, N)
    compare_f.append((value, delta * 0.5))

#### Plot the correlation ####
from matplotlib import pyplot as plt

fig, ax = plt.subplots(figsize=(5, 5))

# Plot the data correlations
ax.scatter(*zip(*compare_cn), s=7, label="myo_cn")
ax.scatter(*zip(*compare_f), s=7, label="myo_f")

# Plot a diagonal
l, h = ax.get_xlim()
ax.plot([l, h], [l, h], '-k', zorder=0)
ax.set_xlim(l, h)
ax.set_ylim(l, h)

# Make axis labels and save figure
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
fig.savefig("ccr_calculate.png")
