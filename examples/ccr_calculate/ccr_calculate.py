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

# print(met.info())
# print(r2)
# numres(i),namres(i),namat(i),obs(i),
# *           tolprot(i),wprot(i)

for atom in prot.get_atoms():
	if atom.name == 'H':
		atomH = atom
		res = atom.parent
		atomN = res['N']
		dd = atomN.dipole_shift_tensor(atomH.position)
		r2 = met.ccr_r2(atomH.position, atomH.gamma, dd)

		line = "{0:3d} {1:4s} {2:4s}  {3:7.3f} {4:5.2f} {5:5.2f}".format(res.id[1], res.resname, atomH.name, r2, 0.0, 0.0)
		print(line)
