from paramagpy import protein, fit, dataparse, metal
import numpy as np
from scipy.optimize import fmin_bfgs

# Load the PDB file
calb = protein.load_pdb('./examples/data_files/2bcb.pdb')
ubi = protein.load_pdb('./examples/data_files/1ubqH.pdb')
myo = protein.load_pdb('./examples/data_files/1bzrH.pdb')

# Load the PCS data
rawDataPCS = dataparse.read_pcs('./examples/data_files/calbindin_Er_HN_PCS.npc')
rawDataPRE = dataparse.read_pre('./examples/data_files/calbindin_Er_H_R2_600.pre')
rawDataRDC = dataparse.read_rdc('./examples/data_files/ubiquitin_s57c_c1_Tb_HN.rdc')
rawDataCCR = dataparse.read_ccr('./examples/data_files/myoglobin_cn.ccr')

# Associate PCS data with atoms of the PDB
pcs = calb.parse(rawDataPCS)
pre = calb.parse(rawDataPRE)
rdc = ubi.parse(rawDataRDC)
ccr = myo.parse(rawDataCCR)

m0 = metal.Metal(position=[6.48348859e-10,  7.19131407e-10, -3.56399432e-10])
# m0 = metal.Metal(position=[25.786E-10,   9.515E-10,   6.558E-10])

initMetals = [m0]
dataArrays = [pcs]
ensembleAverage = False

[mfit], [cal] = fit.svd_gridsearch_fit_metal_from_pcs(initMetals, dataArrays, points=1, radius=0.0, ensembleAverage=ensembleAverage)
















































