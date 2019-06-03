from paramagpy import protein, fit, dataparse, metal

# Load the PDB file
prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# Load the fitted tensor
met = metal.load_tensor('../data_files/calbindin_Er_HN_PCS_tensor.txt')
tmp = metal.Metal()
tmp.set_lanthanide('Er')
met.mueff = tmp.mueff
met.t1e = tmp.t1e
met.taur = 4.5E-9

H = prot[0]['A'][2]['H']
N = prot[0]['A'][2]['N']
dd = N.dipole_shift_tensor(H.position)
r2 = met.ccr_r2(H.position, H.gamma, dd)


print(r2)
