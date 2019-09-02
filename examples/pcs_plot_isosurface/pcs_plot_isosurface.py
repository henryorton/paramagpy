from paramagpy import protein, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# Load the fitted tensor
met = metal.load_tensor('../data_files/calbindin_Er_HN_PCS_tensor.txt')

# Plot the isosurface to be opened in PyMol
met.isomap(prot.id, density=1, isoval=1.0)
