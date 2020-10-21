from paramagpy import protein, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# Load the fitted tensor
met = metal.load_tensor('../data_files/calbindin_Er_HN_PCS_tensor.txt')
met.B0 = 18.8
met.T = 298.0

# Loop over all amide atoms and calulate the RDC value
forFile = []
for atom in prot.get_atoms():
	if atom.name == 'H':
		residue = atom.parent
		seq = residue.id[1]
		if 'N' in residue:
			H = atom
			N = residue['N']
			rdc = met.atom_rdc(H, N)
			line = "{0:2d} {1:^3s} {2:2d} {3:^3s} {4:6.3f} 0.0\n".format(
				seq, H.name, seq, N.name, rdc)
			forFile.append(line)

# Write the calculated RDC values to file
with open("calbindin_Er_RDC_calc.rdc", 'w') as f:
	f.writelines(forFile)