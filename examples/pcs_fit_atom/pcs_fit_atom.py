from paramagpy import protein, fit, dataparse, metal
import re

sites = 'A53C-C2', 'N172C-C2', 'S204C-C2', 'N172C-C12'
ions = 'Tm', 'Tb'
pdb_path = "../data_files/5ev6AH.pdb"
prot = protein.load_pdb(pdb_path)
atomns = 'H07', 'H08'
atoms = [prot[0]['A'][28][atomn] for atomn in atomns]
BOOTSTRAP_ITER = 20

# Isosurface colours
surface_colours = {
	(ions[0], atomns[0]): 'teal',
	(ions[1], atomns[0]): 'blue',
	(ions[0], atomns[1]): 'magenta',
	(ions[1], atomns[1]): 'red',
}

# RMSD contour levels
surface_contour = {
	(ions[0], atomns[0]): 0.04,
	(ions[1], atomns[0]): 0.04,
	(ions[0], atomns[1]): 0.016,
	(ions[1], atomns[1]): 0.02,
}

pmlscript = protein.PyMolScript()
pmlscript.add_pdb(path=pdb_path, name='5ev6')

mdata = []
for site in sites: # Loop over sites
	bindingSite = int(re.search("\\d+", site).group()) # Get residue number
	mStart = metal.Metal()
	mStart.position = prot[0]['A'][bindingSite]['CA'].position # Set strating position

	hnpcss = []
	# Assemble exp. PCS data for both ions
	for ion in ions:
		hnpcs_raw = dataparse.read_pcs("../data_files/IMP1_HN_{}_{}_FREE.npc".format(site, ion))
		hnpcs = prot.parse(hnpcs_raw)
		hnpcss.append(hnpcs)

	# Fit the tensor by SVD, then NLR
	mGuess, _ = fit.svd_gridsearch_fit_metal_from_pcs([mStart, mStart], hnpcss)
	mFit, _ = fit.nlr_fit_metal_from_pcs(mGuess, hnpcss)

	# Sample purturbed tensors by bootstrap
	mSamples, mStd = fit.fit_error_bootstrap(
		fit.nlr_fit_metal_from_pcs, 
		BOOTSTRAP_ITER, 
		0.8, 
		initMetals=mFit, 
		dataArrays=hnpcss
	)

	mdata.append(mSamples)


for ion, mSamples in zip(ions, zip(*mdata)):
	trpdata = []
	mdata = []
	# Loop sites with fitted tensors
	for site, mSample in zip(sites, mSamples):
		# Load TRP PCS data
		trppcs_raw = dataparse.read_pcs("../data_files/IMP1_TRP_{}_{}_FREE.npc".format(site, ion))
		trppcs = prot.parse(trppcs_raw)

		# Assemble associated lists of atoms/PCS with tensors
		for atom in atoms:
			dataselect = trppcs[trppcs['atm'] == atom]
			if len(dataselect)>0:
				for m in mSample: 
					trpdata.append(dataselect)
					mdata.append(m)

	gridsVol = fit.gridsearch_fit_atom_from_pcs(mdata, trpdata, mapSize=20.0, mapDensity=2.0)

	for atom in atoms:
		mapname = "{}{}map".format(atom.id, ion)
		dmapFilePath = "{}.ccp4".format(mapname)
		gridsVol[atom].write(dmapFilePath)
		pmlscript.add_map(
			path=dmapFilePath,
			name=mapname,
			isoVals=[surface_contour[(ion, atom.id)]],
			colours=[surface_colours[(ion, atom.id)]],
			surfaceType='isodot',
		)

pmlscript += "set dot_radius, 0.05"
pmlscript += "show sticks, ////28 and sc."
pmlscript += "show sticks, ////28/CA"
pmlscript += "set bg_rgb=[1,1,1]"
pmlscript += "set mesh_width, 0.5"
pmlscript += "zoom ////28/H07\n"
pmlscript += """
set_view (\
     0.505656540,   -0.827194929,   -0.245069817,\
    -0.741597414,   -0.561904311,    0.366465807,\
    -0.440846384,   -0.003562994,   -0.897575319,\
     0.000152570,    0.000080852,  -36.169487000,\
    48.539413452,   83.819839478,   42.674442291,\
    26.907037735,   45.422363281,  -20.000000000 )
"""
pmlscript += "ray 1600"
pmlscript += "png pcs_fit_atom.png"
pmlscript.write("pcs_fit_atom.pml")
