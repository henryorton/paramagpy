from matplotlib import pyplot as plt
import pylanth as pl


# x, y, z, = (25.786,9.515,6.558)
# ax, rh = (-8.155,-4.913)
# a, b, g = (125.842,142.287,41.759)

# m = pl.metal.Metal.make_tensor(x,y,z,ax,rh,a,b,g,'Er')
# m.taur = 4.25E-9
# m.B0 = 18.8

m0 = pl.metal.Metal.make_tensor(25,10,6,0,0,0,0,0)

fileName = 'ho4icb43G.pdb'
npcName1 = 'ershifts.npc'
npcName2 = 'ybshifts.npc'

prot = pl.protein.load_pdb(fileName)
pcs1 = pl.dataparse.read_pcs(npcName1)
pcs2 = pl.dataparse.read_pcs(npcName2)

dat1 = prot.parse(pcs1)
dat2 = prot.parse(pcs2)

mfit1, mfit2 = pl.fit.fit_metal_from_pcs([m0, m0], [dat1, dat2])

print(mfit2.info())

pl.fit.plot_pcs_fit([mfit2], [dat2])
plt.show()
print(pl.fit.qfactor(mfit2, dat2))


