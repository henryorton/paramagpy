import numpy as np
from metal import Metal
from scipy.optimize import fmin_bfgs
from matplotlib import pyplot as plt
from pprint import pprint


def extract_coords(data):
	atoms, values = zip(*data)
	return np.array([i.coord for i in atoms]), np.array(values)


def sphere_grid(origin, radius, points):
	s = np.linspace(-radius, radius, points)
	mgrid = np.array(np.meshgrid(s, s, s, indexing='ij')).T.reshape(len(s)**3,3)
	norms = np.linalg.norm(mgrid, axis=1)
	sphere_idx = np.where(norms<=radius)
	return mgrid[sphere_idx] + origin


def svd_calc_metal_from_pcs(pos, pcs):
	dist = np.linalg.norm(pos, axis=1)
	x, y, z = pos.T
	a = x**2 - z**2
	b = y**2 - z**2
	c = 2 * x * y
	d = 2 * x * z
	e = 2 * y * z
	mat = (1./(4.*np.pi*dist**5)) * np.array([a,b,c,d,e])
	matinv = np.linalg.pinv(mat)
	sol = matinv.T.dot(pcs)
	calc = mat.T.dot(sol)
	return calc, sol
	

def svd_gridsearch_calc_metal_from_pcs(metals, pcss, 
	origin=None, radius=20.0, points=30):
	if origin is None:
		origin = metals[0].position*1E10

	posarrays = []
	pcsarrays = []
	for pcs in pcss:
		posarray, pcsarray = extract_coords(pcs)
		posarrays.append(posarray*1E-10)
		pcsarrays.append(pcsarray*1E-6)

	sphere = sphere_grid(origin, radius, points)*1E-10

	minscore = 1E30
	print("SVD search started in {} points".format(len(sphere)))
	for pos in sphere:
		score = 0
		for pcsarray, posarray in zip(pcsarrays, posarrays):
			coords = posarray - pos
			calculated, solution = svd_calc_metal_from_pcs(coords, pcsarray)
			score += np.sum((calculated - pcsarray)**2)
		if score<minscore:
			minscore = score
			minpos = pos

	minmetals = [m.copy() for m in metals]
	for pcsarray, posarray, metal in zip(pcsarrays, posarrays, minmetals):
		coords = posarray - minpos
		_, solution = svd_calc_metal_from_pcs(coords, pcsarray)
		metal.position = minpos
		metal.upper_triang = solution

	return minmetals


def nlr_fit_metal_from_pcs(metals, pcss):
	posarrays = []
	pcsarrays = []
	for pcs in pcss:
		posarray, pcsarray = extract_coords(pcs)
		posarrays.append(posarray*1E-10)
		pcsarrays.append(pcsarray)

	def cost(params):
		pos, tensor_params = Metal.unpack_tensor_params(params)
		score = 0.0
		for tensor_param, posarray, pcsarray in zip(tensor_params, posarrays, pcsarrays):
			metal = Metal(position=pos)
			metal.upper_triang = tensor_param
			score += 0.5*np.sum((metal.fast_pcs(posarray) - pcsarray)**2)
		return score

	start_params = Metal.pack_tensor_params(metals)
	opt = fmin_bfgs(cost, start_params)
	fitpos, fit_tensor_params = Metal.unpack_tensor_params(opt)
	fitmetals = [m.copy() for m in metals]
	for fit_tensor_param, m in zip(fit_tensor_params, fitmetals):
		m.position = fitpos
		m.upper_triang = fit_tensor_param

	return fitmetals


def fit_metal_from_pcs(metals, pcss):
	fitmetals = svd_gridsearch_calc_metal_from_pcs(metals, pcss)
	fitmetals = nlr_fit_metal_from_pcs(fitmetals, pcss)
	return fitmetals


def plot_pcs_fit(metals, pcss):
	posarrays = []
	pcsarrays = []
	for pcs in pcss:
		posarray, pcsarray = extract_coords(pcs)
		posarrays.append(posarray*1E-10)
		pcsarrays.append(pcsarray)
	fig = plt.figure(figsize=(5,5))
	ax = fig.add_subplot(111)
	ax.set_xlabel("Experiment")
	ax.set_ylabel("Calculated")
	mini, maxi = min([np.min(i) for i in pcsarrays]), max([np.max(i) for i in pcsarrays])
	ax.plot([mini,maxi], [mini,maxi], '-k', lw=0.5)
	for metal, pos, pcs in zip(metals, posarrays, pcsarrays):
		calc_pcs = metal.fast_pcs(pos)
		ax.plot(pcs, calc_pcs, marker='o', lw=0, ms=3)

	return fig, ax
	

def qfactor(metal, pcs):
	posarray, pcsarray = extract_coords(pcs)
	pcs_calc = metal.fast_pcs(posarray*1E-10)
	numer = np.sum((pcsarray - pcs_calc)**2)
	denom = np.sum(pcsarray**2)
	return (numer/denom)**0.5

def pcs(metal, atom):
	return metal.pcs(atom.coord*1E-10)

def pre(metal, atom, method='dd+dsa+csaxdsa', rtype='r2'):
	pos = atom.coord*1E-10
	gam = atom.gamma
	csa = atom.csa()
	rate = 0.0

	methods = method.split('+')
	if 'dd' in methods and rtype=='r1':
		rate += metal.sbm_r1(pos, gam)
	elif 'dd' in methods and rtype=='r2':
		rate += metal.sbm_r2(pos, gam)

	if 'dsa' in methods and rtype=='r1':
		rate += metal.curie_r1(pos, gam)
	elif 'dsa' in methods and rtype=='r2':
		rate += metal.curie_r2(pos, gam)

	if 'csa' in methods and rtype=='r1':
		rate += metal.curie_r1(pos, gam, csa=csa, ignorePara=True)
	elif 'csa' in methods and rtype=='r2':
		rate += metal.curie_r2(pos, gam, csa=csa, ignorePara=True)

	if ('csaxdsa' in methods or 'dsaxcsa' in methods) and rtype=='r1':
		pre_dsa = metal.curie_r1(pos, gam)
		pre_csa = metal.curie_r1(pos, gam, csa=csa, ignorePara=True)
		pre_cross = metal.curie_r1(pos, gam, csa=csa)
		rate += pre_cross - pre_dsa - pre_csa
	elif ('csaxdsa' in methods or 'dsaxcsa' in methods) and rtype=='r2':
		pre_dsa = metal.curie_r2(pos, gam)
		pre_csa = metal.curie_r2(pos, gam, csa=csa, ignorePara=True)
		pre_cross = metal.curie_r2(pos, gam, csa=csa)
		rate += pre_cross - pre_dsa - pre_csa

	return rate

def rdc(metal, atom1, atom2):
	vector = (atom1.coord - atom2.coord)*1E-10
	return metal.rdc(vector, atom1.gamma, atom2.gamma)
	# vec = (atom1.coord - atom2.coord)*1E-10
	# distance = np.linalg.norm(vec)
	# numer = -metal.HBAR * metal.B0**2 * atom1.gamma * atom2.gamma
	# denom = 120. * metal.K * metal.temperature * np.pi**2
	# preFactor = numer/denom
	# p1 = (1./distance**5)*np.kron(vec,vec).reshape(3,3)
	# p2 = (1./distance**3)*np.identity(3)
	# return preFactor * ((3.*p1 - p2).dot(metal.tensor)).trace()

def ccr(metal, atom):
	pass



from protein import load_pdb
from dataparse import read_pcs

x, y, z, = (25.786,9.515,6.558)
ax, rh = (-8.155,-4.913)
a, b, g = (125.842,142.287,41.759)

m = Metal.make_tensor(x,y,z,ax,rh,a,b,g,'Er')
m.taur = 4.25E-9
m.B0 = 18.8
# m0 = Metal.make_tensor(0,0,0,0,0,0,0,0)

fileName = 'ho4icb43G.pdb'
npcName1 = 'ershifts.npc'
npcName2 = 'ybshifts.npc'

prot = load_pdb(fileName)
# prot = load_pdb('2bcb.pdb')

pcs1 = read_pcs(npcName1)
pcs2 = read_pcs(npcName2)

dat1 = prot.parse(pcs2)

diag = m.tensor_saupe[(0,1,2),(0,1,2)]


atoms, pcss = zip(*dat1)

# for atom1, atom2 in zip(atoms[::2], atoms[1::2]):
# 	pos = atom1.coord*1E-10
# 	print(atom1.parent.id)
# 	print(m.pcs(pos)/ m.pcs2(pos))
	

# mfit = svd_gridsearch_calc_metal_from_pcs([m0], [dat1])
# # mfit = nlr_fit_metal_from_pcs([m, m], [pcs1, pcs2])
# # mfit = fit_metal_from_pcs([m, m], [pcs1, pcs2])
# mfit = fit_metal_from_pcs([m0], [dat1])

# plot_pcs_fit(mfit, [dat1])
# print(mfit[0].info())

# plt.show()
# print(mfit[0].info())
# print(qfactor(prot, mfit[0], pcs1))
# print(mfit[0].info())

# plot_pcs_fit(prot, [m], [pcs])



# print(mfit[0].info())

# a = data[0][0]

# print(a.omega)

# print(p[1], structureBuilder=Structure)

# atom = p.struct[0]['A'][56]['H']

# dat = p.parse_npc(npcName)


# pos = p.struct[0]['A'][56]['H'].coord

# m = Metal.make_tensor(x,y,z,ax,rh,a,b,g)
# m0 = Metal(position=atom.coord*1E-10)


# i = atom.get_full_id()

# l = p[i]

# print(l.coord)
# p[i].set_coord([1,2,3])
# print(l.coord)


# posarray = p.get_coords_by_id(list(dat.keys())[:10])*1E-10
# pcsarray = np.array(list(dat.values())[:10])*1E-6
# m1, m2 = svd_calc_metal_from_pcs(posarray, pcsarray)

# g = fit_metal_from_pcs(p, [m0], [dat])

# print(g[0].info())

# plot_pcs_fit(p, g, [dat])

# plt.show()

# for i in g:
# 	print(i.position)
# 	print(i.axrh)

# print(g.position*1E10)
# print(g.axrh*1E32)

# pos = p.get_coords_by_id(dat.keys())*1E-10
# pcs = np.array(list(dat.values()))*1E-6












