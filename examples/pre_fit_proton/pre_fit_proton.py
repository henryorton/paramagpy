from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# Load the PRE data
rawData_er = dataparse.read_pre('../data_files/calbindin_Er_H_R2_600.pre')
rawData_tb = dataparse.read_pre('../data_files/calbindin_Tb_H_R2_800.pre')

# Associate PRE data with atoms of the PDB
data_er = prot.parse(rawData_er)
data_tb = prot.parse(rawData_tb)

# Load PCS fitted tensor
mStart_er = metal.load_tensor('../data_files/calbindin_Er_HN_PCS_tensor.txt')
mStart_tb = metal.load_tensor('../data_files/calbindin_Tb_HN_PCS_tensor.txt')
mStart_er.B0 = 14.1
mStart_tb.B0 = 18.8

# Fit the rotational correlation time using non-linear regression
(m_er,), (cal_er,), qfac = fit.nlr_fit_metal_from_pre(
	[mStart_er], [data_er], params=['taur'], rtypes=['r2'])
(m_tb,), (cal_tb,), qfac = fit.nlr_fit_metal_from_pre(
	[mStart_tb], [data_tb], params=['taur'], rtypes=['r2'])

# # Save the fitted tensor to file
m_er.save('calbindin_Er_H_R2_600_tensor.txt')
m_tb.save('calbindin_Tb_H_R2_800_tensor.txt')

#### Plot the correlation ####
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(5,5))

# Unpack the experimental values
atm_er, exp_er, err_er = zip(*data_er)
atm_tb, exp_tb, err_tb = zip(*data_tb)

# Plot the data
ax.plot(exp_er, cal_er, marker='o', lw=0, ms=3, c='r',
	label="Er: taur = {:3.1f} ns".format(1E9*m_er.taur))
ax.plot(exp_tb, cal_tb, marker='o', lw=0, ms=3, c='g',
	label="Tb: taur = {:3.1f} ns".format(1E9*m_tb.taur))

# Plot a diagonal
l, h = ax.get_ylim()
ax.plot([l,h],[l,h],'-k',zorder=0)
ax.set_xlim(l,h)
ax.set_ylim(l,h)

# Make axis labels and save figure
ax.set_xlabel("Experiment")
ax.set_ylabel("Calculated")
ax.legend()
fig.savefig("pre_fit_proton.png")