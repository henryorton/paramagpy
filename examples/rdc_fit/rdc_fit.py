from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/1ubqH.pdb')

# Load the PCS data
rawData1 = dataparse.read_rdc('../data_files/ubiquitin_a28c_c1_Tb_HN.rdc')
rawData2 = dataparse.read_rdc('../data_files/ubiquitin_s57c_c1_Tb_HN.rdc')

# Associate PCS data with atoms of the PDB
parsedData1 = prot.parse(rawData1)
parsedData2 = prot.parse(rawData2)

# Define an initial tensor
mStart1 = metal.Metal(B0=18.8, temperature=308.0)
mStart2 = metal.Metal(B0=18.8, temperature=308.0)

# Calculate the tensor using SVD
mFit1, calc1, qfac1 = fit.svd_fit_metal_from_rdc(mStart1, parsedData1)
mFit2, calc2, qfac2 = fit.svd_fit_metal_from_rdc(mStart2, parsedData2)

# Save the fitted tensor to file
mFit1.save('ubiquitin_a28c_c1_Tb_tensor.txt')
mFit2.save('ubiquitin_s57c_c1_Tb_tensor.txt')

import numpy as np

vecs = []
gamm = []

for atom1, atom2, value, error in parsedData1:
	vec = atom1.position-atom2.position
	gam = atom1.gamma * atom2.gamma
	vecs.append(vec)
	gamm.append(gam)
	rdc = mFit1.rdc(vec, gam)
	print(rdc)

print(np.array(vecs).shape)

rdcs = mFit1.fast_rdc(np.array(vecs),np.array(gamm))
print(rdcs)