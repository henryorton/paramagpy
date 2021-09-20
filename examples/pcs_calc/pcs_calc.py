from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# Read tensor
met = metal.load_tensor('../data_files/calbindin_Tb_HN_PCS_tensor.txt')

# Open a file to write 
with open("pcs_calc.npc", 'w') as f:
	# Loop over atoms in PDB and calculate PCS
	for atom in prot.get_atoms():
		data = {
			'name':atom.name,
			'seq':atom.parent.id[1],
			'pcs':met.atom_pcs(atom)
		}
		line = "{seq:4d} {name:5s} {pcs:8.5f} 0.0\n".format(**data)
		f.write(line)
