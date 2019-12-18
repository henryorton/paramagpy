from paramagpy import protein, fit, dataparse, metal
import numpy as np

# Load the PDB file
calb = protein.load_pdb('./examples/data_files/2bcb.pdb')
ubi = protein.load_pdb('./examples/data_files/1ubqH.pdb')
myo = protein.load_pdb('./examples/data_files/1bzrH.pdb')

# Load the PCS data
rawDataPCS = dataparse.read_pcs('./examples/data_files/calbindin_Er_HN_PCS_errors_trunc.npc')
rawDataPRE = dataparse.read_pre('./examples/data_files/calbindin_Er_H_R2_600.pre')
rawDataRDC = dataparse.read_rdc('./examples/data_files/ubiquitin_s57c_c1_Tb_HN.rdc')
rawDataCCR = dataparse.read_ccr('./examples/data_files/myoglobin_cn.ccr')

# Associate PCS data with atoms of the PDB
pcs = calb.parse(rawDataPCS)
pre = calb.parse(rawDataPRE)
rdc = ubi.parse(rawDataRDC)
ccr = myo.parse(rawDataCCR)

a = fit.extract_atom_data(pcs)
b = fit.extract_atom_data(pre)
c = fit.extract_rdc_data(rdc)
d = fit.extract_ccr_data(ccr)

print(c)