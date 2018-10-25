import sys
sys.path.append('../')

from paramagpy import protein, fit, dataparse, metal


prot = protein.load_pdb('ho4icb43G.pdb')
rawData = dataparse.read_pcs('calbindin_Er_HN_PCS.npc')

parsedData = prot.parse(rawData)

a = parsedData[0][0]


mStart = metal.Metal()
mStart.set_lanthanide('Er')

mGuess = fit.svd_gridsearch_fit_metal_from_pcs(
	[mStart],[parsedData], radius=20, points=10)



# mFit = fit.nlr_fit_metal_from_pre([m0], [pres], ['x','y','z','ax','rh'], sumIndices=None, 
# 	rtype='r2', usesbm=True, usedsa=True, usecsa=False, progress=None)

# print(m[0].position)



# dater = prot.parse(pcser)

# pprint(dater)


# # atoms, pcs, errs = zip(*dater)
# # summs = np.array([i.serial_number for i in atoms])
# # poss = np.array([i.coord*1E-10 for i in atoms])
# # pcss = np.array(pcs)*1E-6

# x, y, z, = (25.786,9.515,6.558)
# # x, y, z, = (6.673,7.289,-3.288)
# m0 = metal.Metal.make_tensor(x,y,z,0,0,0,0,0)
# # m0.B0 = 18.8
# # pos = np.array([x,y,z])*1E-10


# # pars = ['x','y','z','ax','rh','a','b','g']

# guess = fit.svd_gridsearch_fit_metal_from_pcs([m0],[dater], radius=0, points=1)
# # mfit = fit.nlr_fit_metal_from_pcs([m0],[dater],pars)

# print(guess[0].info())

# m = mfit[0]
# m.set_utr()

# print(m.info())
# # print(guess[0].info())





# x, y, z, = (25.786,9.515,6.558)
# m0 = metal.Metal.make_tensor(x,y,z,40,20,0.1,0.2,0.5,'Er')
# m0.B0_MHz = 600.0
# m0.taur = 4.5E-9

