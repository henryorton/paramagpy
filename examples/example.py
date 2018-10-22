from matplotlib import pyplot as plt
import numpy as np
from pprint import pprint
from pylanth import protein, fit, dataparse, metal
import time






prot = protein.load_pdb('ho4icb43G.pdb')
# # prot = protein.load_pdb('2bcb.pdb')

# x, y, z, = (25.786,9.515,6.558)
# m0 = metal.Metal.make_tensor(x,y,z,40,20,0.1,0.2,0.5,'Er')
# m0.B0_MHz = 600.0
# m0.taur = 4.5E-9


# poss = []
# gammas = []
# for a in prot.get_atoms():
# 	if a.name =='H':
# 		poss.append(a.coord*1E-10)
# 		gammas.append(a.gamma)
# poss = np.array(poss)
# gammas = np.array(gammas)

# preer = dataparse.read_pre('erpre.pre')
# pres = prot.parse(preer)

# m = fit.nlr_fit_metal_from_pre([m0], [pres], ['x','y','z','ax','rh'], sumIndices=None, 
# 	rtype='r2', usesbm=True, usedsa=True, usecsa=False, progress=None)

# print(m[0].position)




pcser = dataparse.read_pcs('ershifts_errors.npc')
# # pcser = dataparse.read_rdc('errdcs.rdc')

dater = prot.parse(pcser)

pprint(dater)


# atoms, pcs, errs = zip(*dater)
# summs = np.array([i.serial_number for i in atoms])
# poss = np.array([i.coord*1E-10 for i in atoms])
# pcss = np.array(pcs)*1E-6

x, y, z, = (25.786,9.515,6.558)
# x, y, z, = (6.673,7.289,-3.288)
m0 = metal.Metal.make_tensor(x,y,z,0,0,0,0,0)
# m0.B0 = 18.8
# pos = np.array([x,y,z])*1E-10


# pars = ['x','y','z','ax','rh','a','b','g']

guess = fit.svd_gridsearch_fit_metal_from_pcs([m0],[dater], radius=0, points=1)
# mfit = fit.nlr_fit_metal_from_pcs([m0],[dater],pars)

print(guess[0].info())

# m = mfit[0]
# m.set_utr()

# print(m.info())
# # print(guess[0].info())










# ax    | 1E-32 m^3 :    -6.989
# rh    | 1E-32 m^3 :    -4.246
# x     |   1E-10 m :     6.511
# y     |   1E-10 m :     7.213
# z     |   1E-10 m :    -3.587
# a     |       deg :   162.446
# b     |       deg :    60.780
# g     |       deg :    17.499
# mueff |        Bm :     0.000
# shift |       ppm :     0.000
# B0    |         T :    18.790
# temp  |         K :   298.150
# t1e   |        ps :     0.000


# 2bcb 1,2,3,4 er_errors
# ax    | 1E-32 m^3 :    -7.080
# rh    | 1E-32 m^3 :    -4.469
# x     |   1E-10 m :     6.673
# y     |   1E-10 m :     7.289
# z     |   1E-10 m :    -3.288
# a     |       deg :   164.018
# b     |       deg :    65.619
# g     |       deg :    12.028
# mueff |        Bm :     0.000
# shift |       ppm :     0.000
# B0    |         T :    18.790
# temp  |         K :   298.150
# t1e   |        ps :     0.000



# for a in prot.get_atoms():
# 	if a.name=='H':
# 		res = a.parent
# 		b = res['N']
# 		atoms.append((a, b))

# out = []
# for a, b in atoms:
# 	vec = (b.coord - a.coord)*1E-10
# 	rdc = m.rdc(vec, a.gamma, b.gamma)
# 	_, _, chn, (_, seq, _), (atm1, _) = a.get_full_id()
# 	_, _, chn, (_, seq, _), (atm2, _) = b.get_full_id()
# 	line = "{0:<3} {1:<3} {2:<3} {3:<3} {4:8.3f} {5:3.1f}\n".format(seq, atm1, seq, atm2, rdc, 0.0)
# 	out.append(line)

# with open("errdcs.rdc", 'w') as o:
# 	for line in out:
# 		o.write(line)


# prot = protein.load_pdb('ho4icb43G.pdb')
# rdcer = dataparse.read_rdc('errdcs.rdc')
# dater = prot.parse(rdcer)

# print(rdcer)


# m0 = metal.Metal.make_tensor(0,0,0,0,0,0,0,0)
# m = metal.Metal.make_tensor(0,0,0,-8.688,-4.192,116.011,138.058,-136.508)
# m0.B0 = 18.79

# mfit = fit.svd_fit_metal_from_rdc(m0,dater)

# print(mfit.info())


# vecarray, gamarray, rdcarray, errarray = fit.extract_rdc(dater)
# vecarray *= 1E-10

# dist = np.linalg.norm(vecarray, axis=1)

# pfarray = -3*(m0.MU0 * gamarray * m0.HBAR) / (8 * np.pi**2  * dist**5)

# x, y, z = vecarray.T
# a = x**2 - z**2
# b = y**2 - z**2
# c = 2 * x * y
# d = 2 * x * z
# e = 2 * y * z
# mat = pfarray * np.array([a,b,c,d,e])
# sol = mat.T.dot(m.tensor_saupe[m.upper_coords])

# print(sol/rdcarray)



# print(sol)
# fitMetal.upper_triang_saupe = sol








# m = mfit[0]


# m.axrh[1] = 0.0
# m.eulers[2] = 0.0


# a = np.linalg.eigh(m.tensor)
# print(a)


# m = mfit[0]
# print(m)
# m.axrh[1] = 0.0
# m.eulers[2] = 0.0
# print(m)



# for p in poss[:8]:
# 	rads = mfit[0].rads(p)
# 	print(rads)

# print(mfit[0].fast_rads(poss[:8]))


# print(guess[0].info())
# print(fit.qfactor(guess[0], dater))
# print(mfit[0].info())
# print(fit.qfactor(mfit[0], dater))



# atoms = atoms[:8]
# posss = poss[:8] - m0.position

# csa = prot[0]['A'][5]['N'].csa()
# csat = m0.copy()
# csat.tensor = csa

# numbat_racs = m0.racs(csat.tensor)


# print(m0.rads(poss[2]))
# print(m0.fast_rads(poss[0:10]))


# Tr[S A St U B Ut]
# Tr[S S A B Ut Ut]
# Tr[A B (S Ut)**2]



	


# print(csa)
# print(csat.tensor)
# print(csat.eigenvalues)

# numbat_racs = m0.racs((3./2.)*csat.tensor)


# pf = m0.B0**2 / (5*m0.MU0*m0.K*m0.temperature)

# cp = csat.rotationMatrix.T
# mp = m0.rotationMatrix.T



# def dotty(i,j):
# 	v, w = cp[i], mp[j]
# 	return v.dot(w)

# tmp = []
# for i in range(3):
# 	for j in range(3):
# 		tmp.append(dotty(i,j))

# a = np.array(tmp).reshape(3,3)
# print(a)
# print(cp.dot(mp.T))

# racs = 0.0

# for i in range(3):
# 	for j in range(3):
# 		tmp = csat.eigenvalues[i]*m0.eigenvalues[j]*dotty(i,j)**2
# 		racs += pf * tmp

# print(racs*1E6)
# print(numbat_racs)
# print(racs*1E6/numbat_racs)














# t = m0.fast_pcs_rads(poss[:8])
# print(t)


# distance = np.linalg.norm(posss, axis=1)
# preFactor = 1./(4.*np.pi)
# p1 = (1./distance[:,None,None]**5)*np.einsum('ij,ik->ijk', posss, posss)
# p2 = (1./distance[:,None,None]**3)*np.identity(3)
# tmp = (preFactor * (3.*p1 - p2)).dot(m0.tensor).trace(axis1=1,axis2=2)/3.



# a = np.array([i+10 for i in range(2*2)]).reshape(2,2)
# b = np.array([i for i in range(3*2*2)]).reshape(3,2,2)

# t = np.einsum('jk,ikl->ijl',a,b)

# print(a.dot(b[2]))

# print(t)






# pos = np.array(position, dtype=float) - self.position
# distance = np.linalg.norm(pos)
# preFactor = 1./(4.*np.pi)
# p1 = (1./distance**5)*np.kron(pos,pos).reshape(3,3)
# p2 = (1./distance**3)*np.identity(3)
# return (preFactor * (3.*p1 - p2)).dot(self.tensor)


# pos = np.array(position, dtype=float) - self.position
# distance = np.linalg.norm(pos)
# preFactor = 1./(4.*np.pi)
# p1 = (1./distance**5)*np.kron(pos,pos).reshape(3,3)
# p2 = (1./distance**3)*np.identity(3)
# return (preFactor * (3.*p1 - p2)).dot(self.tensor)



# calc1, sol1 = fit.svd_calc_metal_from_pcs(poss-pos, pcss)
# calc2, sol2 = fit.svd_calc_metal_from_pcs_offset(poss-pos, pcss)



# mfit1 = m0.copy()
# mfit1.upper_triang = sol1

# mfit2 = m0.copy()
# mfit2.upper_triang = sol2[:-1]
# mfit2.shift = sol2[-1]*1E26*1E6
# print(mfit2.info())


# mfit3 = fit.svd_gridsearch_calc_metal_from_pcs([m0], [dater], 
# 	origin=None, radius=0.0, points=1, offset_shift=True)

# print(mfit3[0].info())

# mfiter = fit.fit_metal_from_pcs([m0], [dater])

# print(mfiter[0].info())


# x, y, z, = (25.786,9.515,6.558)
# ax, rh = (-8.155,-4.913)
# a, b, g = (125.842,142.287,41.759)

# m = metal.Metal.make_tensor(x,y,z,ax,rh,a,b,g,'Er')
# m.taur = 4.25E-9
# m.B0 = 18.8

# m0 = metal.Metal.make_tensor(25,10,6,0,0,0,0,0)

# fileName = 'ho4icb43G.pdb'
# npcName1 = 'ershifts.npc'
# npcName2 = 'ybshifts.npc'

# prot = protein.load_pdb(fileName)
# pcs1 = dataparse.read_pcs(npcName1)
# pcs2 = dataparse.read_pcs(npcName2)

# pprint(pcs1)

# dat1 = prot.parse(pcs1)
# dat2 = prot.parse(pcs2)

# mfit1, mfit2 = fit.fit_metal_from_pcs([m0, m0], [dat1, dat2])

# print(mfit2.info())

# fit.t_pcs_fit([mfit2], [dat2])
# .show()
# print(fit.qfactor(mfit2, dat2))


