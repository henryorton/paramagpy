from paramagpy import protein, fit, dataparse, metal
import numpy as np


	
# Load the PDB file
# prot = protein.load_pdb('../data_files/4icbH_mut.pdb')
prot = protein.load_pdb('1bzrH.pdb')



ironAtom = prot[0]['A'][("H_HEM",154," ")]['FE']
met = metal.Metal(position=ironAtom.position)
met.B0 = 18.79
met.T = 273.0 + 30.0

met.iso = 30.1E-32
# met.iso = 4.4E-32
met.taur = 5.7E-9
# met.axrh = (4.24E-32, 0.0)


compare = []
for i, val in vals.items():
# for a in prot.get_atoms():
	# if a.name=='H':
		# try:
			# H = a
			# N = a.parent['N']
		# except KeyError:
			# continue

	H = prot[0]['A'][i]['H']
	N = prot[0]['A'][i]['N']

		# vec = H.position - N.position
		# dd = N.dipole_shift_tensor(H.position)
		# metvec = H.position - met.position

		# theta = np.arccos(vec.dot(metvec)/(np.linalg.norm(vec)*np.linalg.norm(metvec)))

		# pf = -(1.)/(15.)
		# pf*= met.MU0/(np.pi*4)
		# pf*= (met.B0 * H.gamma**2 * N.gamma * met.HBAR)/(np.linalg.norm(vec)**3)
		# pf*= -(3*met.iso)/(4*np.pi*np.linalg.norm(metvec)**3)
		# pf*= (3.*np.cos(theta)**2-1.)/2.
		# pf*= 4*met.taur + (3*met.taur)/(1+(met.taur*met.B0*H.gamma)**2)


	dd = N.dipole_shift_tensor(H.position)
	r2 = met.ccr_r2(H.position, H.gamma, dd)
	compare.append((r2,val))
	# print(r2/pf)



from matplotlib import pyplot as plt
x, y = zip(*compare)
plt.scatter(x, y)
plt.show()


