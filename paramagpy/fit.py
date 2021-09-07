import numpy as np
from scipy.optimize import fmin_bfgs
import warnings
from collections import OrderedDict

structdtype = np.dtype([
	('mdl', int          ),
	('atm', object       ),
	('pos', float,  3    ),
	('exp', float        ),
	('cal', float        ),
	('err', float        ),
	('idx', int          ),
	('gam', float        ),
	('xsa', float, (3,3) )])

def extract_atom_data(data, csa=False, separateModels=True):
	"""
	Extract values required for PCS/PRE calculations

	Parameters
	----------
	data : numpy.ndarray
		a numpy structured array containing atomic information
		and experimental data values
		This is returned from.
		:meth:`paramagpy.protein.CustomStructure.parse`
	csa : bool, optional
		when True, calculates the CSA tensor for each atom
		this may be required for RACS and CSAxDSA calculations
	separateModels : bool, optional
		when True, separates data into separate lists 
		by their model number. When False, returns only one list
	Returns
	-------
	arr : numpy.ndarra
		this has fields specified by structdtype
		this array is core to all fitting algorithms
	"""
	arr = np.empty(len(data), dtype=structdtype)

	for name in ('mdl', 'atm', 'exp', 'cal'):
		arr[name] = data[name]
	_, arr['idx'] = np.unique(data['idx'], return_inverse=True)

	if 0.0 in data['err']:
		arr['err'] = np.ones(len(data['err']))
		warnings.warn("0.0 value uncertainty. All values weighted evenly")
	else:
		arr['err'] = data['err']

	arr['pos'] = [a.position for a in data['atm']]
	arr['gam'] = [a.gamma for a in data['atm']]
	if csa:
		arr['xsa'] = [a.csa for a in data['atm']]

	if separateModels:
		return [arr[arr['mdl']==m] for m in np.unique(arr['mdl'])]
	else:
		return [arr]

def extract_rdc_data(data, separateModels=True):
	"""
	Extract values required for RDC calculations

	Parameters
	----------
	data : numpy.ndarray
		a numpy structured array containing atomic information
		and experimental data values
		This is returned from 
		:meth:`paramagpy.protein.CustomStructure.parse`
	separateModels : bool, optional
		when True, separates data into separate lists 
		by their model number. When False, returns only one list
	Returns
	-------
	arr : numpy.ndarra
		this has fields specified by structdtype
		this array is core to all fitting algorithms
	"""
	arr = np.empty(len(data), dtype=structdtype)

	for name in ('mdl', 'atm', 'exp', 'cal'):
		arr[name] = data[name]
	_, arr['idx'] = np.unique(data['idx'], return_inverse=True)

	if 0.0 in data['err']:
		arr['err'] = np.ones(len(data['err']))
		warnings.warn("0.0 value uncertainty. All values weighted evenly")
	else:
		arr['err'] = data['err']

	arr['pos'] = [a2.position - a1.position for a1, a2 in data[['atm','atx']]]
	arr['gam'] = [a1.gamma * a2.gamma for a1, a2 in data[['atm','atx']]]

	if separateModels:
		return [arr[arr['mdl']==m] for m in np.unique(arr['mdl'])]
	else:
		return [arr]

def extract_ccr_data(data, separateModels=True):
	"""
	Extract values required for CCR calculations

	Parameters
	----------
	data : numpy.ndarray
		a numpy structured array containing atomic information
		and experimental data values
		This is returned from 
		:meth:`paramagpy.protein.CustomStructure.parse`
	separateModels : bool, optional
		when True, separates data into separate lists 
		by their model number. When False, returns only one list
	Returns
	-------
	arr : numpy.ndarra
		this has fields specified by structdtype
		this array is core to all fitting algorithms
	"""
	arr = np.empty(len(data), dtype=structdtype)

	for name in ('mdl', 'atm', 'exp', 'cal'):
		arr[name] = data[name]
	_, arr['idx'] = np.unique(data['idx'], return_inverse=True)

	if 0.0 in data['err']:
		arr['err'] = np.ones(len(data['err']))
		warnings.warn("0.0 value uncertainty. All values weighted evenly")
	else:
		arr['err'] = data['err']

	arr['pos'] = [a.position for a in data['atm']]
	arr['gam'] = [a.gamma for a in data['atm']]
	arr['xsa'] = [p.dipole_shift_tensor(a.position) for a,p in data[['atm','atx']]]

	if separateModels:
		return [arr[arr['mdl']==m] for m in np.unique(arr['mdl'])]
	else:
		return [arr]

def sphere_grid(origin, radius, points):
	"""
	Make a grid of cartesian points within a sphere

	Parameters
	----------
	origin : float
		the centre of the sphere
	radius : float
		the radius of the sphere
	points : int
		the number of points per radius

	Returns
	-------
	array : array of [x,y,z] coordinates
		the points within the sphere
	"""
	s = np.linspace(-radius, radius, 2*points-1)
	mgrid = np.array(np.meshgrid(s, s, s, indexing='ij')).T.reshape(len(s)**3,3)
	norms = np.linalg.norm(mgrid, axis=1)
	sphere_idx = np.where(norms<=radius)
	return mgrid[sphere_idx] + origin


def svd_calc_metal_from_pcs(pos, pcs, idx, errors):
	"""
	Solve PCS equation by single value decomposition.
	This function is generally called by higher methods like 
	<svd_gridsearch_fit_metal_from_pcs>

	Parameters
	----------
	pos : array of [x,y,z] floats
		the atomic positions in meters
	pcs : array of floats
		the PCS values in ppm. Note these should be weighted by the
		experimental uncertainties.
	idx : array of ints
		an index assigned to each atom. Common indices determine summation
		between models for ensemble averaging.
	errors : array of floats
		the standard deviation representing experimental uncertainty
		in the measured value

	Returns
	-------
	calc : array of floats
		the calculated PCS values from the fitted tensor
	sol : array of floats
		solution to the linearised PCS equation and 
		consists of the tensor 5 matrix elements
	"""
	floatscale = 1E-24
	dist = np.linalg.norm(pos, axis=1)
	x, y, z = pos.T
	a = x**2 - z**2
	b = y**2 - z**2
	c = 2 * x * y
	d = 2 * x * z
	e = 2 * y * z
	mat = (1./(4.*np.pi*dist**5)) * np.array([a,b,c,d,e]) / errors
	mat = np.array([np.bincount(idx, weights=col) for col in mat])*floatscale
	matinv = np.linalg.pinv(mat)
	sol = matinv.T.dot(pcs*1E-6)
	calc = mat.T.dot(sol)*1E6
	return calc, sol*floatscale


def svd_calc_metal_from_pcs_offset(pos, pcs, idx, errors):
	"""
	Solve PCS equation by single value decomposition with offset.
	An offset arising from referencing errors between diamagnetic
	and paramagnetic datasets can be accounted for using this method.
	This function is generally called by higher methods like 
	<svd_gridsearch_fit_metal_from_pcs>

	NOTE: the factor of 1E26 is required for floating point error mitigation

	Parameters
	----------
	pos : array of [x,y,z] floats
		the atomic positions in meters
	pcs : array of floats
		the PCS values in ppm. Note these should be weighted by the
		experimental uncertainties.
	idx : array of ints
		an index assigned to each atom. Common indices determine summation
		between models for ensemble averaging.
	errors : array of floats
		the standard deviation representing experimental uncertainty
		in the measured value

	Returns
	-------
	calc : array of floats
		the calculated PCS values from the fitted tensor
	sol : array of floats
		solution to the linearised PCS equation and 
		consists of the tensor 5 matrix elements and offset values
	"""
	floatscale = 1E-24
	dist = np.linalg.norm(pos, axis=1)
	x, y, z = pos.T
	a = x**2 - z**2
	b = y**2 - z**2
	c = 2 * x * y
	d = 2 * x * z
	e = 2 * y * z
	scale = 1./(4.*np.pi*dist**5)
	mat = scale * np.array([a,b,c,d,e,1E26/scale]) / errors
	mat = np.array([np.bincount(idx, weights=col) for col in mat])*floatscale
	matinv = np.linalg.pinv(mat)
	sol = matinv.T.dot(pcs*1E-6)
	calc = mat.T.dot(sol)*1E6
	sol[-1] *= 1E26
	return calc, sol*floatscale


def svd_gridsearch_fit_metal_from_pcs(initMetals, dataArrays, ensembleAverage=False,
	origin=None, radius=20.0, points=16, offsetShift=False, progress=None):
	"""
	Fit deltaChi tensor to PCS values using Single Value Decomposition over
	a grid of points in a sphere.
	Note this uses a weighted SVD fit which takes into account 
	experimental errors

	Parameters
	----------
	initMetals : list of Metal objects
		a list of metals used as starting points for fitting. 
		a list must always be provided, but may also contain 
		only one element. If multiple metals are provided, each metal
		is fitted to their respective PCS dataset by index in <dataArrays, 
		but all are fitted to a common position.
	dataArrays : list of PCS dataArray
		each PCS dataArray must correspond to an associated metal for fitting.
		each PCS dataArray has structure determined by 
		:meth:`paramagpy.protein.CustomStructure.parse`.
	ensembleAverage : bool, optional
		when False, each model of the structure is fit independently.
		The parameters for each fitted tensor are then averaged before 
		returning the final averaged tensor.
		When True, the structure models are treated as an ensemble and 
		ensemble averaging of calculated PCS/PRE/RDC/CCR values is 
		conducted at all stages of fitting to fit a single tensor to all
		models simultaneously. The 'idx' column of the dataArray
		determines the ensemble averaging behaviour with common indices 
		for atoms between models resulting in their summation.
	origin : float, optional
		the centre of the gridsearch of positions in Angstroms.
		If None, the position of the first metal is used
	radius : float, optional
		the radius of the gridsearch in Angstroms.
	points : int, optional
		the number of points per radius in the gridsearch
	offsetShift : bool, optional
		if True, an offset value added to all PCS values is included in
		the SVD fitting. This may arise due to a referencing error between
		diamagnetic and paramagnetic PCS datasets and may be used when
		many data points are available.
		Default False, no offset is included in the fitting.
	progress : object, optional
		to keep track of the calculation, progress.set(x) is called each
		iteration and varies from 0.0 -> 1.0 when the calculation is complete.

	Returns
	-------
	fitMetals : list of metals
		a list of the fitted tensors.
	dataArrays : list of dataArray
		each dataArray is copy of the original dataArray with
		the 'cal' column populated with back-calculated values from the 
		fitted tensor.
	"""
	if len(initMetals)!=len(dataArrays):
		raise ValueError("<metals> and <dataArrays> must have same length")

	if offsetShift:
		svd_func = svd_calc_metal_from_pcs_offset
	else:
		svd_func = svd_calc_metal_from_pcs

	datas = {}
	metalAvgs = []
	for metal, dataArray in zip(initMetals, dataArrays):
		metalAvg = []
		if ensembleAverage:
			d = extract_atom_data(dataArray, csa=False, 
				separateModels=False)[0]
			m = metal.copy()
			m.par['weightSum'] = np.bincount(d['idx'], weights=d['exp']/d['err'])
			metalAvg.append(m)
			if 0 not in datas:
				datas[0] = []
			datas[0].append((m, d))
		else:
			for d in extract_atom_data(dataArray, csa=False, 
				separateModels=True):
				mdl = d['mdl'][0]
				if mdl not in datas:
					datas[mdl] = []
				m = metal.copy()
				m.par['weightSum'] = np.bincount(d['idx'], weights=d['exp']/d['err'])
				metalAvg.append(m)
				datas[mdl].append((m, d))
		metalAvgs.append(metalAvg)

	if origin is None:
		origin = initMetals[0].position*1E10

	sphere = sphere_grid(origin, radius, points)*1E-10
	tot = len(sphere)*len(datas)
	prog = 0.0
	if tot<1:
		raise ValueError("Zero grid points selected for SVD search")
	print("SVD gridsearch started in {} points".format(len(sphere)))

	for mdl in datas:
		minscore = 1E308
		data = datas[mdl]
		for pos in sphere:
			if progress:
				prog += 1
				progress.set(prog/tot)
			score = 0.0
			sols = []
			for m, d in data:
				calculated, solution = svd_func(d['pos']-pos, m.par['weightSum'], d['idx'], d['err'])
				sols.append(solution)
				score += np.sum((calculated - m.par['weightSum'])**2)
			if score<minscore:
				minscore = score
				for sol, (m, _) in zip(sols, data):
					m.position = pos
					if offsetShift:
						m.upper_triang = sol[:-1]
						m.shift = sol[-1]*1E6
					else:
						m.upper_triang = sol

	fitMetals = []
	for metalAvg in metalAvgs:
		mAvg = metalAvg[0].copy()
		mAvg.average(metalAvg)
		mAvg.set_utr()
		mAvg.par['mav'] = metalAvg
		fitMetals.append(mAvg)

	newDataArrays = []
	for m, dataArray in zip(fitMetals, dataArrays):
		d = extract_atom_data(dataArray, 
			separateModels=False)[0]
		newDataArray = dataArray.copy()
		newDataArray['cal'] = m.fast_pcs(d['pos'])
		newDataArrays.append(newDataArray)

	if progress:
		progress.set(1.0)

	return fitMetals, newDataArrays


def nlr_fit_metal_from_pcs(initMetals, dataArrays, 
	params=('x','y','z','ax','rh','a','b','g'), ensembleAverage=False,
	userads=False, useracs=False, progress=None):
	"""
	Fit deltaChi tensor to PCS values using non-linear regression.

	Parameters
	----------
	initMetals : list of Metal objects
		a list of metals used as starting points for fitting. 
		a list must always be provided, but may also contain 
		only one element. If multiple metals are provided, each metal
		is fitted to their respective PCS dataset by index in <dataArrays, 
		but all are fitted to a common position.
	dataArrays : list of PCS dataArray
		each PCS dataArray must correspond to an associated metal for fitting.
		each PCS dataArray has structure determined by 
		:meth:`paramagpy.protein.CustomStructure.parse`.
	params : list of str
		the parameters to be fit. 
		For example ['x','y','z','ax','rh','a','b','g','shift']
	ensembleAverage : bool, optional
		when False, each model of the structure is fit independently.
		The parameters for each fitted tensor are then averaged before 
		returning the final averaged tensor.
		When True, the structure models are treated as an ensemble and 
		ensemble averaging of calculated PCS/PRE/RDC/CCR values is 
		conducted at all stages of fitting to fit a single tensor to all
		models simultaneously. The 'idx' column of the dataArray
		determines the ensemble averaging behaviour with common indices 
		for atoms between models resulting in their summation.
	userads : bool, optional
		include residual anisotropic dipolar shielding (RADS) during fitting
	useracs : bool, optional
		include residual anisotropic chemical shielding (RACS) during fitting.
		CSA tensors are taken using the <csa> method of atoms.
	progress : object, optional
		to keep track of the calculation, progress.set(x) is called each
		iteration and varies from 0.0 -> 1.0 when the calculation is complete.

	Returns
	-------
	fitMetals : list of metals
		a list of the fitted tensors.
	dataArrays : list of dataArray
		each dataArray is copy of the original dataArray with
		the 'cal' column populated with back-calculated values from the 
		fitted tensor.
	"""
	if len(initMetals)!=len(dataArrays):
		raise ValueError("initMetals and dataArrays must have same length")

	datas = {}
	metalAvgs = []
	for metal, dataArray in zip(initMetals, dataArrays):
		metalAvg = []
		if ensembleAverage:
			d = extract_atom_data(dataArray, csa=useracs, 
				separateModels=False)[0]
			m = metal.copy()
			metalAvg.append(m)
			if 0 not in datas:
				datas[0] = []
			datas[0].append((m, d))
		else:
			for d in extract_atom_data(dataArray, csa=useracs, 
				separateModels=True):
				mdl = d['mdl'][0]				
				if mdl not in datas:
					datas[mdl] = []
				m = metal.copy()
				metalAvg.append(m)
				datas[mdl].append((m, d))
		metalAvgs.append(metalAvg)

	params = set(params)
	pospars = tuple(params & set(['x','y','z']))
	otherpars = tuple(params - set(['x','y','z']))

	for mdl in datas:
		data = datas[mdl]
		for i, (m, d) in enumerate(data):
			m.par['pos'] = slice(0, len(pospars))
			m.par['oth'] = slice(len(pospars) + i*len(otherpars), 
								  len(pospars) + (i+1)*len(otherpars))

	tot = len(datas)
	prog = 0.0
	for mdl in datas:
		data = datas[mdl]
		startpars = data[0][0].get_params(pospars)
		for m, _ in data:
			startpars += m.get_params(otherpars)
		def cost(args):
			score = 0.0
			for m, d in data:
				m.set_params(zip(pospars, args[m.par['pos']]))
				m.set_params(zip(otherpars, args[m.par['oth']]))
				cal = m.fast_pcs(d['pos'])
				if userads:
					cal += m.fast_rads(d['pos'])
				if useracs:
					cal += m.fast_racs(d['xsa'])
				diff = (cal - d['exp']) / d['err']
				selectiveSum = np.bincount(d['idx'], weights=diff)
				score += np.sum(selectiveSum**2)

			return score

		fmin_bfgs(cost, startpars, disp=False)
		if progress:
			prog += 1
			progress.set(prog/tot)

	fitMetals = []
	for metalAvg in metalAvgs:
		mAvg = metalAvg[0].copy()
		mAvg.average(metalAvg)
		mAvg.set_utr()
		mAvg.par['mav'] = metalAvg
		fitMetals.append(mAvg)

	newDataArrays = []
	for m, dataArray in zip(fitMetals, dataArrays):
		d = extract_atom_data(dataArray, csa=useracs, 
			separateModels=False)[0]
		newDataArray = dataArray.copy()
		newDataArray['cal'] = m.fast_pcs(d['pos'])
		if userads:
			newDataArray['cal'] += m.fast_rads(d['pos'])
		if useracs:
			newDataArray['cal'] += m.fast_racs(d['xsa'])
		newDataArrays.append(newDataArray)

	if progress:
		progress.set(1.0)

	return fitMetals, newDataArrays




def nlr_fit_metal_from_pre(initMetals, dataArrays, rtypes, params=('x','y','z'), 
	usesbm=True, usegsbm=False, usedsa=True, 
	usecsa=False, ensembleAverage=False, progress=None):
	"""
	Fit Chi tensor to PRE values using non-linear regression.

	Parameters
	----------
	initMetals : list of Metal objects
		a list of metals used as starting points for fitting. 
		a list must always be provided, but may also contain 
		only one element. If multiple metals are provided, each metal
		is fitted to their respective PRE dataset by index in <dataArrays, 
		but all are fitted to a common position.
	dataArrays : list of PRE dataArray
		each PRE dataArray must correspond to an associated metal for fitting.
		each PRE dataArray has structure determined by 
		:meth:`paramagpy.protein.CustomStructure.parse`.
	rtypes : list of str, optional
		the relaxtion type, either 'r1' or 'r2'. A list must be provided with
		each element corresponding to an associated dataset.
		Defaults to 'r2' for all datasets of None is specified.
	params : list of str
		the parameters to be fit. 
		For example ['x','y','z','ax','rh','a','b','g','shift']
	usesbm : bool, optional
		include Solomon-Bloemenbergen-Morgan (Dipole-dipole) relaxation theory.
		default is True
	usegsbm : bool, optional
		include anisotropic dipolar relaxation theory.
		note that the g-tensor must be set for this 
		default is False
	usedsa : bool, optional
		include Dipolar-Shielding-Anisotropy (Curie Spin) relaxation theory.
		default is True
	usecsa : bool, optional
		include Chemical-Shift-Anisotropy cross-correlated realxation theory.
		default is False
	ensembleAverage : bool, optional
		when False, each model of the structure is fit independently.
		The parameters for each fitted tensor are then averaged before 
		returning the final averaged tensor.
		When True, the structure models are treated as an ensemble and 
		ensemble averaging of calculated PCS/PRE/RDC/CCR values is 
		conducted at all stages of fitting to fit a single tensor to all
		models simultaneously. The 'idx' column of the dataArray
		determines the ensemble averaging behaviour with common indices 
		for atoms between models resulting in their summation.
	progress : object, optional
		to keep track of the calculation, progress.set(x) is called each
		iteration and varies from 0.0 -> 1.0 when the calculation is complete.

	Returns
	-------
	fitMetals : list of metals
		a list of the fitted tensors.
	dataArrays : list of dataArray
		each dataArray is copy of the original dataArray with
		the 'cal' column populated with back-calculated values from the 
		fitted tensor.
	"""
	if len(initMetals)!=len(dataArrays)!=len(rtypes):
		raise ValueError("initMetals, dataArrays and rtypes must have same length")

	if len(set(rtypes) | set(['r1','r2']))>2:
		raise TypeError("rtype must be a list with values 'r1' and 'r2' only")

	datas = {}
	metalAvgs = []
	for metal, dataArray, rtype in zip(initMetals, dataArrays, rtypes):
		metalAvg = []
		if ensembleAverage:
			d = extract_atom_data(dataArray, csa=usecsa, 
				separateModels=False)[0]
			m = metal.copy()
			m.par['rtp'] = rtype
			metalAvg.append(m)
			if 0 not in datas:
				datas[0] = []
			datas[0].append((m, d))
		else:
			for d in extract_atom_data(dataArray, csa=usecsa, 
				separateModels=True):
				mdl = d['mdl'][0]				
				if mdl not in datas:
					datas[mdl] = []
				m = metal.copy()
				m.par['rtp'] = rtype
				metalAvg.append(m)
				datas[mdl].append((m, d))
		metalAvgs.append(metalAvg)

	params = set(params)
	pospars = tuple(params & set(['x','y','z']))
	otherpars = tuple(params - set(['x','y','z']))

	for mdl in datas:
		data = datas[mdl]
		for i, (m, d) in enumerate(data):
			m.par['pos'] = slice(0, len(pospars))
			m.par['oth'] = slice(len(pospars) + i*len(otherpars), 
								  len(pospars) + (i+1)*len(otherpars))

	for mdl in datas:
		data = datas[mdl]
		startpars = data[0][0].get_params(pospars)
		for m, _ in data:
			startpars += m.get_params(otherpars)

		def cost(args):
			score = 0.0
			for m, d in data:
				m.set_params(zip(pospars, args[m.par['pos']]))
				m.set_params(zip(otherpars, args[m.par['oth']]))
				cal = m.fast_pre(d['pos'], d['gam'], m.par['rtp'], 
				dsa=usedsa, sbm=usesbm, gsbm=usegsbm, csaarray=d['xsa'])
				diff = (cal - d['exp']) / d['err']
				selectiveSum = np.bincount(d['idx'], weights=diff)
				score += np.sum(selectiveSum**2)
			return score

		fmin_bfgs(cost, startpars, disp=False)

	fitMetals = []
	for metalAvg in metalAvgs:
		mAvg = metalAvg[0].copy()
		mAvg.par['rtp'] = metalAvg[0].par['rtp']
		mAvg.par['mav'] = metalAvg
		mAvg.average(metalAvg)
		mAvg.set_utr()
		mAvg.metalAvg = metalAvg
		fitMetals.append(mAvg)

	newDataArrays = []
	for m, dataArray in zip(fitMetals, dataArrays):
		d = extract_atom_data(dataArray, 
			csa=usecsa, separateModels=False)[0]
		newDataArray = dataArray.copy()
		newDataArray['cal'] = m.fast_pre(d['pos'], d['gam'], m.par['rtp'], 
				dsa=usedsa, sbm=usesbm, gsbm=usegsbm, csaarray=d['xsa'])
		newDataArrays.append(newDataArray)

	if progress:
		progress.set(1.0)

	return fitMetals, newDataArrays


def svd_calc_metal_from_rdc(vec, rdc_parameterised, idx, errors):
	"""
	Solve RDC equation by single value decomposition.
	This function is generally called by higher methods like 
	<svd_fit_metal_from_rdc>

	Parameters
	----------
	vec : array of [x,y,z] floats
		the internuclear vectors in meters
	rdc_parameterised : array of floats
		the experimental RDC values, normalised by a prefactor
	idx : array of ints
		an index assigned to each atom. Common indices determine summation
		between models for ensemble averaging.
	errors : array of floats
		the standard deviation representing experimental uncertainty
		in the measured value

	Returns
	-------
	calc : array of floats
		the calculated RDC values from the fitted tensor
	sol : array of floats
		sol is the solution to the linearised PCS equation and 
		consists of the tensor matrix elements
	"""
	dist = np.linalg.norm(vec, axis=1)
	x, y, z = vec.T
	a = x**2 - z**2
	b = y**2 - z**2
	c = 2 * x * y
	d = 2 * x * z
	e = 2 * y * z
	mat = (1./dist**5) * np.array([a,b,c,d,e]) / errors
	matSum = np.array([np.bincount(idx, weights=col) for col in mat])
	matinv = np.linalg.pinv(matSum)
	sol = matinv.T.dot(rdc_parameterised)
	calc = matSum.T.dot(sol)
	return calc, sol


def svd_fit_metal_from_rdc(initMetals, dataArrays,
	params=('ax','rh','a','b','g'), ensembleAverage=False, progress=None):
	"""
	Fit deltaChi tensor to RDC values using Single Value Decomposition.
	Note this is a weighted SVD calculation which takes into account
	experimental errors.

	Parameters
	----------
	initMetals : list of Metal objects
		a list of metals used as starting points for fitting. 
		a list must always be provided, but may also contain 
		only one element. If multiple metals are provided, each metal
		is fitted to their respective RDC dataset by index in <dataArrays>.
	dataArrays : list of PRE dataArray
		each RDC dataArray must correspond to an associated metal for fitting.
		each RDC dataArray has structure determined by 
		:meth:`paramagpy.protein.CustomStructure.parse`.
	params : list of str
		the parameters to be fit. 
		NOTE: This is a dummy argument and does not influence the fitting.
		The default parameters ('ax','rh','a','b','g') are the only option.
	ensembleAverage : bool, optional
		when False, each model of the structure is fit independently.
		The parameters for each fitted tensor are then averaged before 
		returning the final averaged tensor.
		When True, the structure models are treated as an ensemble and 
		ensemble averaging of calculated PCS/PRE/RDC/CCR values is 
		conducted at all stages of fitting to fit a single tensor to all
		models simultaneously. The 'idx' column of the dataArray
		determines the ensemble averaging behaviour with common indices 
		for atoms between models resulting in their summation.
	progress : object, optional
		to keep track of the calculation, progress.set(x) is called each
		iteration and varies from 0.0 -> 1.0 when the calculation is complete.

	Returns
	-------
	fitMetals : list of metals
		a list of the fitted tensors.
	dataArrays : list of dataArray
		each dataArray is copy of the original dataArray with
		the 'cal' column populated with back-calculated values from the 
		fitted tensor.
	"""
	if len(initMetals)!=len(dataArrays):
		raise ValueError("initMetals and dataArrays must have same length")

	datas = {}
	metalAvgs = []
	for metal, dataArray in zip(initMetals, dataArrays):
		metalAvg = []
		if ensembleAverage:
			d = extract_rdc_data(dataArray, separateModels=False)[0]
			m = metal.copy()
			metalAvg.append(m)
			if 0 not in datas:
				datas[0] = []
			datas[0].append((m, d))
		else:
			for d in extract_rdc_data(dataArray, separateModels=True):
				mdl = d['mdl'][0]				
				if mdl not in datas:
					datas[mdl] = []
				m = metal.copy()
				metalAvg.append(m)
				datas[mdl].append((m, d))
		metalAvgs.append(metalAvg)

	tot = len(datas)
	prog = 0.0
	for mdl in datas:
		data = datas[mdl]
		for m, d in data:
			pfarray = -3*(m.MU0 * d['gam'] * m.HBAR) / (8 * np.pi**2)
			rdc_parameterised = np.bincount(d['idx'], weights=d['exp'] / (pfarray * d['err']))
			calculated, solution = svd_calc_metal_from_rdc(d['pos'], rdc_parameterised, d['idx'], d['err'])
			m.upper_triang_alignment = solution

		if progress:
			prog += 1
			progress.set(prog/tot)

	fitMetals = []
	for metalAvg in metalAvgs:
		mAvg = metalAvg[0].copy()
		mAvg.average(metalAvg)
		mAvg.set_utr()
		mAvg.par['mav'] = metalAvg
		fitMetals.append(mAvg)

	newDataArrays = []
	for m, dataArray in zip(fitMetals, dataArrays):
		d = extract_rdc_data(dataArray, separateModels=False)[0]
		newDataArray = dataArray.copy()
		newDataArray['cal'] = m.fast_rdc(d['pos'], d['gam'])
		newDataArrays.append(newDataArray)

	if progress:
		progress.set(1.0)

	return fitMetals, newDataArrays


def nlr_fit_metal_from_ccr(initMetals, dataArrays, params=('x','y','z'), 
	ensembleAverage=False, progress=None):
	"""
	Fit Chi tensor to CCR values using non-linear regression.
	This algorithm applies to CSA/Curie spin cross-correlated relaxation
	for R2 differential line broadening.

	Parameters
	----------
	initMetals : list of Metal objects
		a list of metals used as starting points for fitting. 
		a list must always be provided, but may also contain 
		only one element. If multiple metals are provided, each metal
		is fitted to their respective PRE dataset by index in <dataArrays, 
		but all are fitted to a common position.
	dataArrays : list of PRE dataArray
		each PRE dataArray must correspond to an associated metal for fitting.
		each PRE dataArray has structure determined by 
		:meth:`paramagpy.protein.CustomStructure.parse`.
	params : list of str
		the parameters to be fit. 
		For example ['x','y','z','ax','rh','a','b','g']
	ensembleAverage : bool, optional
		when False, each model of the structure is fit independently.
		The parameters for each fitted tensor are then averaged before 
		returning the final averaged tensor.
		When True, the structure models are treated as an ensemble and 
		ensemble averaging of calculated PCS/PRE/RDC/CCR values is 
		conducted at all stages of fitting to fit a single tensor to all
		models simultaneously. The 'idx' column of the dataArray
		determines the ensemble averaging behaviour with common indices 
		for atoms between models resulting in their summation.
	progress : object, optional
		to keep track of the calculation, progress.set(x) is called each
		iteration and varies from 0.0 -> 1.0 when the calculation is complete.

	Returns
	-------
	fitMetals : list of metals
		a list of the fitted tensors.
	dataArrays : list of dataArray
		each dataArray is copy of the original dataArray with
		the 'cal' column populated with back-calculated values from the 
		fitted tensor.
	"""
	if len(initMetals)!=len(dataArrays):
		raise ValueError("initMetals and dataArrays must have same length")

	datas = {}
	metalAvgs = []
	for metal, dataArray in zip(initMetals, dataArrays):
		metalAvg = []
		if ensembleAverage:
			d = extract_ccr_data(dataArray, separateModels=False)[0]
			m = metal.copy()
			metalAvg.append(m)
			if 0 not in datas:
				datas[0] = []
			datas[0].append((m, d))
		else:
			for d in extract_ccr_data(dataArray, separateModels=True):
				mdl = d['mdl'][0]				
				if mdl not in datas:
					datas[mdl] = []
				m = metal.copy()
				metalAvg.append(m)
				datas[mdl].append((m, d))
		metalAvgs.append(metalAvg)

	params = set(params)
	pospars = tuple(params & set(['x','y','z']))
	otherpars = tuple(params - set(['x','y','z']))

	for mdl in datas:
		data = datas[mdl]
		for i, (m, d) in enumerate(data):
			m.par['pos'] = slice(0, len(pospars))
			m.par['oth'] = slice(len(pospars) + i*len(otherpars), 
								  len(pospars) + (i+1)*len(otherpars))

	tot = len(datas)
	prog = 0.0
	for mdl in datas:
		data = datas[mdl]
		startpars = data[0][0].get_params(pospars)
		for m, _ in data:
			startpars += m.get_params(otherpars)
		def cost(args):
			score = 0.0
			for m, d in data:
				m.set_params(zip(pospars, args[m.par['pos']]))
				m.set_params(zip(otherpars, args[m.par['oth']]))
				cal = m.fast_ccr(d['pos'], d['gam'], d['xsa'])
				diff = (cal - d['exp']) / d['err']
				selectiveSum = np.bincount(d['idx'], weights=diff)
				score += np.sum(selectiveSum**2)

			return score

		fmin_bfgs(cost, startpars, disp=False)
		if progress:
			prog += 1
			progress.set(prog/tot)

	fitMetals = []
	for metalAvg in metalAvgs:
		mAvg = metalAvg[0].copy()
		mAvg.average(metalAvg)
		mAvg.set_utr()
		mAvg.par['mav'] = metalAvg
		fitMetals.append(mAvg)

	newDataArrays = []
	for m, dataArray in zip(fitMetals, dataArrays):
		d = extract_ccr_data(dataArray, separateModels=False)[0]
		newDataArray = dataArray.copy()
		newDataArray['cal'] = m.fast_ccr(d['pos'], d['gam'], d['xsa'])
		newDataArrays.append(newDataArray)

	if progress:
		progress.set(1.0)

	return fitMetals, newDataArrays


def metal_standard_deviation(metals, params):
	"""
	Calculate the standard deviation in parameters <params> for a
	list of metal objects <metals>.

	Parameters
	----------
	metals : list of Metal objects
		the metals for which the standard deviation in parameters
		will be calculated
	params : list of str
		the parameters for the standard deviation calculation. 
		For example ['x','y','z','ax','rh','a','b','g','shift']

	Returns
	-------
	std_metal : metal object
		the returned metal object has attributes equal to the
		standard deviation in the given parameter.
		All other attributes are zero.
	"""
	std_metals = []
	all_param_values = []
	for metal in metals:
		all_param_values.append(metal.get_params(params))

	std_params = {}
	for param, values in zip(params, zip(*all_param_values)):
		std_params[param] = np.std(values)

	std_metal = metal.__class__(temperature=0.0, B0=0.0)
	std_metal.set_params(std_params.items())
	return std_metal


def fit_error_models(fittingFunction, **kwargs):
	"""
	Perform uncertainty analysis sourcing noise from cooridinates as defined
	by models of the PDB structure.
	This function takes a fitting routine <fittingFunction> and repeats it for
	each model. The standard deviation in the fitted parameters is then returned.

	Parameters
	----------
	fittingFunction : function
		the fitting routine to be used. 
		This could be 'nlr_fit_metal_from_ccr' for example
	kwargs : dict
		all key-word arguments will be bundled into this variable and
		parsed to the fittingFunction.

	Returns
	-------
	sample_metals : list of list of metals
		the metals fitted to the data with noise at each iteration
	std_metals : list of metals
		the standard deviation in fitted parameters over all iterations of the
		Monte Carlo simulation.
		These are stored within the metal object. All unfitted parameters 
		are zero.
	"""
	if kwargs.get('ensembleAverage'):
		raise ValueError("""You cannot define ensemble averaging 
			when sourcing noise from models of the structure as 
			they require separate fits""")

	fitMetals, _ = fittingFunction(**kwargs)

	if kwargs.get('progress'):
			kwargs['progress'].set(1.0)

	if 'params' not in kwargs:
		kwargs['params'] = fittingFunction.__defaults__[0]

	sampleMetals = [m.par['mav'] for m in fitMetals]
	stdMetals = []
	for metals in sampleMetals:
		stdMetals.append(metal_standard_deviation(metals, kwargs['params']))

	return sampleMetals, stdMetals



def fit_error_monte_carlo(fittingFunction, iterations, **kwargs):
	"""
	Perform uncertainty analysis sourcing noise from experimental uncertainties
	This function takes a fitting routine <fittingFunction> and repeats it for
	the specified iterations in a Monte-Carlo approach. With each iteration,
	random noise sourced from a uniform distribution scaled by the experimental
	uncertainties is added to the experimental values. The standard deviation in 
	the fitted parameters is then returned.

	NOTE: the 'err' column of the dataArrays must be set to non-zero values for
	this method to work.

	Parameters
	----------
	fittingFunction : function
		the fitting routine to be used. 
		This could be 'nlr_fit_metal_from_ccr' for example
	iterations : int
		the number of iterations for the Monte-Carlo simulation
	kwargs : dict
		all key-word arguments will be bundled into this variable and
		parsed to the fittingFunction.

	Returns
	-------
	sample_metals : list of list of metals
		the metals fitted to the data with noise at each iteration
	std_metals : list of metals
		the standard deviation in fitted parameters over all iterations of the
		Monte Carlo simulation.
		These are stored within the metal object. All unfitted parameters 
		are zero.
	"""
	initMetals = kwargs['initMetals']
	dataArrays = kwargs['dataArrays']
	sampleMetals = []
	for i in range(iterations):
		newInitMetals = []
		newDataArrays = []
		for m, d in zip(initMetals, dataArrays):
			newInitMetals.append(m.copy())
			d['exp'] += d['err'] * np.random.uniform(low=-1, high=1, size=len(d))
			newDataArrays.append(d)

		kwargs['initMetals'] = newInitMetals
		kwargs['dataArrays'] = newDataArrays

		metals, _ = fittingFunction(**kwargs)

		sampleMetals.append(metals)

		if kwargs.get('progress'):
			kwargs['progress'].set(float(i+1)/iterations)

	sampleMetals = list(zip(*sampleMetals))
	stdMetals = []

	if 'params' not in kwargs:
		kwargs['params'] = fittingFunction.__defaults__[0]

	for metals in sampleMetals:
		stdMetals.append(metal_standard_deviation(metals, kwargs['params']))

	return sampleMetals, stdMetals


def fit_error_bootstrap(fittingFunction, iterations, fraction, **kwargs):
	"""
	Perform uncertainty analysis sourcing noise from fractioning the 
	experimental data. This function takes a fitting routine <fittingFunction> 
	and repeats it for the specified iterations in a Bootstrap approach. 
	With each iteration, a random subset of the experimental data is sampled as
	specified by the <fraction> argument. The standard deviation in the fitted
	parameters is then returned.

	Parameters
	----------
	fittingFunction : function
		the fitting routine to be used. 
		This could be 'nlr_fit_metal_from_ccr' for example
	iterations : int
		the number of iterations for the Monte-Carlo simulation
	fraction : float
		the proportion of data to be samples. Must be between 0 and 1.0
	kwargs : dict
		all key-word arguments will be bundled into this variable and
		parsed to the fittingFunction.

	Returns
	-------
	sample_metals : list of list of metals
		the metals fitted to the data with noise at each iteration
	std_metals : list of metals
		the standard deviation in fitted parameters over all iterations of the
		Monte Carlo simulation.
		These are stored within the metal object. All unfitted parameters 
		are zero.
	"""
	if not (0.0<fraction<1.0):
		raise ValueError("The bootstrap sample fraction must be between 0 and 1")

	initMetals = kwargs['initMetals']
	dataArrays = kwargs['dataArrays']
	sampleMetals = []
	for i in range(iterations):
		newInitMetals = []
		newDataArrays = []
		for m, d in zip(initMetals, dataArrays):
			newInitMetals.append(m.copy())
			d = np.random.choice(d, int(fraction*len(d)), replace=False)
			newDataArrays.append(d)

		kwargs['initMetals'] = newInitMetals
		kwargs['dataArrays'] = newDataArrays

		metals, _ = fittingFunction(**kwargs)

		sampleMetals.append(metals)

		if kwargs.get('progress'):
			kwargs['progress'].set(float(i+1)/iterations)

	if 'params' not in kwargs:
		kwargs['params'] = fittingFunction.__defaults__[0]

	sampleMetals = list(zip(*sampleMetals))
	stdMetals = []
	for metals in sampleMetals:
		stdMetals.append(metal_standard_deviation(metals, kwargs['params']))

	return sampleMetals, stdMetals


def qfactor(dataArray, ensembleAverage=False, calDenominator=False):
	"""
	Calculate the Q-factor to judge tensor fit quality

	A lower value indicates a better fit. The Q-factor is calculated using
	the following equation:

	.. math::
		Q = \\sqrt{
			\\frac{\\sum_i\\left[\\left(\\sum_m\\left[
			PCS^{exp}_{m,i}-PCS^{calc}_{m,i}\\right]\\right)^2\\right]}
			{\\sum_i\\left[
			\\left(\\sum_m\\left[PCS^{exp}_{m,i}\\right]\\right)^2\\right]}
		}

	where :math:`m` and :math:`i` are usually indexed over models and atoms
	respectively.

	Parameters
	----------
	dataArray : numpy array
		the dataArray must contain the columns 'exp', 'cal' and 'idx'
		corresponding to the experimenta, calculated and index values
		respectively. The index value determines the ensemble averaging
		behaviour, and can be ignored if the argument <ensembleAverage>
		is False.
	ensembleAverage : bool, optional
		when False, the q-factor calculation squares each difference
		independently.
		When True, the q-factor calculates an ensemble average before
		taking the square of differences. The 'idx' column of the dataArray
		determines the ensemble averaging behaviour with common indices 
		for atoms between models resulting in their summation.
	calDenominator : bool, optional
		when False, the standard Q-factor is calculated with only the sum
		of squares for the experimental values used in the denominator
		when True, the Q-factor established by Ubbink et al. is 
		calculated which has a sum of absolute values of exp and cal
		values squared in the denominator.

	Returns
	-------
	qfactor : float
		the Q-factor
	"""
	if len(dataArray)==0:
		return np.nan
	diff = dataArray['exp'] - dataArray['cal']
	if ensembleAverage:
		numer = np.sum(np.bincount(dataArray['idx'], weights=diff)**2)
		if calDenominator:
			tmp = np.abs(dataArray['exp']) + np.abs(dataArray['cal'])
		else:
			tmp = dataArray['exp']
		denom = np.sum(np.bincount(dataArray['idx'], weights=tmp)**2)
	else:
		numer = np.sum(diff**2)
		if calDenominator:
			tmp = np.abs(dataArray['exp']) + np.abs(dataArray['cal'])
		else:
			tmp = dataArray['exp']
		denom = np.sum(tmp**2)
	return (numer/denom)**0.5


def ensemble_average(dataArray):
	"""
	Calculate the ensemble average for the calculated values in the 
	column 'cal' of the argument <dataArray> over models of the
	PDB file.
	Ensemble averaging behaviour is determined by the column 'idx'
	of the input array.

	Parameters
	----------
	dataArray : numpy array
		the input array for ensemble averaging

	Returns
	-------
	data : numpy array
		a smaller dataArray with ensemble averaged values
	"""
	mdl = sorted(set(dataArray['mdl']))[0]
	data = dataArray.copy()[dataArray['mdl']==mdl]
	for row in data:
		d = dataArray[dataArray['idx']==row['idx']]
		row['cal'] = np.mean(d['cal'])

	return data




class DensityMap(object):
	"""
	A class to help with calculations on a grid.
	The grid can be setup by specifying the origin (centre)
	edge size and density.
	Arguments are given in Angstrom, but stored in metres
	The grid positions are available for calculations
	The density can be written out as a .ccp4 electron density map
	Use this class for any contour plots to view in PyMol
	"""
	CCP4_HEADER_DTYPE = np.dtype([
	    ('np', 'i4', 3),       # Points in each axis
	    ('mode', 'i4'),        # 2 for 32-bit float
	    ('ns', 'i4', 3),       # Start number for each axis
	    ('ni', 'i4', 3),       # Number of intervals along each axis
	    ('dim', 'f4', 3),      # Cell dimensions in Angstrom
	    ('ang', 'f4', 3),      # Cell angles in degrees
	    ('map', 'i4', 3),      # Axis mappings, 1=x,2=y,3=z
	    ('amin', 'f4'),        # Minimum density value
	    ('amax', 'f4'),        # Maximum density value
	    ('amean', 'f4'),       # Mean density value
	    ('ispg', 'i4'),        # space group number
	    ('padd', 'V932')       # additional headers
	])
	def __init__(self, origin, size, density):
		"""
		Make a density map over a cubic grid of coordinates

		Parameters
		----------
		origin : np.ndarray of floats
			[x,y,z] position in Angstrom defining the centre of the cubic grid
		size : float
			the grid edge size in Angstrom
		density : float
			grid density defined as points per angstrom
		"""
		origin_vertex = np.asarray(
			density * (origin - size/2.0), dtype=int)
		low = origin_vertex / float(density)
		high = low + size
		points = np.array([int(density*size)]*3) + 1
		domains = [1E-10*np.linspace(*i) for i in zip(low, high, points)]
		
		self.shape = points
		self.mesh = np.array(np.meshgrid(*domains, indexing='ij')).T
		self.positions = self.mesh.reshape(np.prod(points),3)
		self.density = np.zeros(np.prod(points), dtype=np.float64)

		h = np.zeros(1, dtype=self.CCP4_HEADER_DTYPE)
		h['np'] = points
		h['mode'] = 2
		h['ns'] = origin_vertex
		h['ni'] = points - 1
		h['dim'] = high - low
		h['ang'] = [90.0, 90.0, 90.0]
		h['map'] = [1, 2, 3]
		h['ispg'] = 1

		self.header = h

	def minpos(self):
		"""
		Fetch the grid point position with minimum density value

		Returns
		-------
		position : np.ndarray
			[x,y,z] position in metres of the point
		"""
		return self.positions[np.argmin(self.density)]

	def boundary_min(self):
		"""
		Does the minimum value lie on the boundary of the grid?

		Return
		------
		truth: bool
			return true if the minimum density is located
			on the boundary of the grid
		"""
		minarg = np.argmin(self.density.reshape(*self.shape))
		idx = np.array(np.unravel_index(minarg, self.shape))
		for i, s in zip(idx, self.shape):
			if i==0 or i==(s-1):
				return True
		return False

	def get_positions_below_density(self, cutoffValue):
		"""
		Returns all positions in the grid that are below the 
		argument <cutoffValue>

		Parameters
		----------
		cutoffValue : float
			the threshold below which points will be taken

		Returns
		-------
		positions : numpy.ndarray
			a list of positions which are below the cutoff value
		"""
		return self.positions[np.where(self.density < cutoffValue)]

	def write(self, fileName):
		"""
		Write density values to file in .ccp4 format
		Headers are accounted for from the given geometry of the grid

		Parameters
		----------
		fileName : string
			the name of the .ccp4 file
		"""
		self.header['amin'] = np.min(self.density)
		self.header['amax'] = np.max(self.density)
		self.header['amean'] = np.mean(self.density)
		with open(fileName, 'wb') as f:
			f.write(self.header.tobytes())
			f.write(np.asarray(self.density, dtype=np.float32).tobytes())



def gridsearch_fit_atom_from_pcs(metals, dataArrays, mapSize=10.0, mapDensity=1.0):
	"""
	Calculate likely regions for an atom on a grid using
	an experimental PCS value and multiple delta-chi-tensors.
	The calculation returns a grid of PCS RMSD values for each atom
	The smallest values on the grid indicated the likely positions
	This function returns a dictionary with key/value pairs
	associating the atoms/grids.

	Parameters
	----------
	metals : list of Metal objects
		a list of metals used for calculating the theoretical
		PCS values used during the RMSD calculation 
		a list must always be provided, but may also contain 
		only one element.
	dataArrays : list of PCS dataArray
		each PCS dataArray must correspond to an associated metal.
		each PCS dataArray has structure determined by 
		:meth:`paramagpy.protein.CustomStructure.parse`.
	mapSize : float, optional
		the edge length of the grid in angstrom
		defaults to 10 Angstrom
	mapDensity : float, optional
		the density of points in the grid
		defaults to 1 point per Angstrom

	Returns
	-------
	positions : dict of :class:`paramagpy.fit.DensityMap`
		a dictionary of density maps. Each key is an atom
		as defined in <dataArrays> and corresponds to a value
		which is the density map.
		a density map defines the RMSD calculation on a grid
	"""
	data = {}
	# Sort data for separate calculations for each atom
	for metal, dataArray in zip(metals, dataArrays):
		for row in dataArray:
			if row['atm'] not in data:
				data[row['atm']] = []
			data[row['atm']].append((metal, row['exp']))

	out = {}
	# Loop over each atom
	for atom in data:
		dm = DensityMap(atom.position*1E10, mapSize, mapDensity)
		# Loop over each metal and calculate sum of squares
		for metal, exp in data[atom]:
			dm.density += (metal.fast_pcs(dm.positions) - exp)**2

		# Average and square root to finish RMSD calculation
		dm.density = (dm.density / len(data[atom]))**0.5 

		# Display a warning if minimum point is on the grid boundary
		# Might be necessary to change box definition for grid
		if dm.boundary_min():
			warnings.warn("Atom {} has minimum position solved on grid boundary.".format(atom))

		out[atom] = dm

	return out


def gridsearch_fit_atom_restrain_distance(densityMapA, densityMapB, distUpper, distLower, number):
	"""
	Given two RMSD density maps, this function will compare
	all points pairwise and return only those within the bounds of
	a given distance cutoff and within a certain number of points 
	that are that are sorted by RMSD value.
	This might be useful if two density maps for separate atoms
	in a molecule are known to be constrained w.r.t. one another
	and you would like to use that restraint to further restrict 
	the space of PCS RMSD points.
	The calculation first sorts the RMSD points and takes the bottom
	<number> of points and then compares each point pariwise to 
	fulfill the distance condition. It then returns those points
	from both maps. Unfortunately there is no correlation data
	available between these two maps.

	Parameters
	----------
	densityMapA: :class:`paramagpy.fit.DensityMap`
		a density map of PCS RMSD values.
	densityMapB: :class:`paramagpy.fit.DensityMap`
		a second density map of PCS RMSD values.
	distUpper : float
		The upper limit for distance.
		Any pairwise distances larger than this value
		will be rejected from the final space of points
	distLower : float
		The lower distance limit for distance
		Any pairwise distance smaller than this value
		will be rejectet from the final space of points
	number : int
		The total number of positions to be considered
		for the pairwise distance comparison.
		This calculation first sorts points by RMSD
		and takes <number> of points with minimum RMSD
		and uses these for the pairwise distance calculation.
		Note that the total number of points returned
		could be significnalty more than this value after
		the pairwise comparison

	Returns
	-------
	tuple : np.ndarray of position coordinates
		two lists of [x,y,z] coordinates are returned associated
		with the inputs <densityMapA> and <densityMapB>.
		The returned coordinates are taken from the original grids
		and represent points that have another associated point in
		the other grid which is within the distance bounds and 
		contained with an RMSD low enough to be within the lowest
		<number> of points
	"""
	idxSortA = np.argsort(densityMapA.density)[0:number]
	idxSortB = np.argsort(densityMapB.density)[0:number]
	posA = densityMapA.positions[idxSortA]
	posB = densityMapB.positions[idxSortB]
	dist = np.linalg.norm(posA[:,None] - posB, axis=2)
	idxA, idxB = np.where(np.logical_and(dist < distUpper, dist > distLower))
	idxSort = np.argsort(densityMapA.density[idxA] + densityMapB.density[idxB])
	return posA[idxA[idxSort]], posB[idxB[idxSort]]


def gridsearch_fit_atom_restrain_distance_cutoff(densityMapA, densityMapB, distUpper, distLower, cutoffValue):
	"""
	Given two RMSD density maps, this function will compare
	all points pairwise and return only those within the bounds of
	a given distance cutoff and within a given RMSD cutoff.
	This might be useful if two density maps for separate atoms
	in a molecule are known to be constrained w.r.t. one another
	and you would like to use that restraint to further restrict 
	the space of PCS RMSD points.
	The calculation first selects points with an RMSD less than
	the given cutoff value and then compares each point pariwise to 
	fulfill the distance condition. It then returns those points
	from both maps. Unfortunately there is no correlation data
	available between these two maps.

	Parameters
	----------
	densityMapA: :class:`paramagpy.fit.DensityMap`
		a density map of PCS RMSD values.
	densityMapB: :class:`paramagpy.fit.DensityMap`
		a second density map of PCS RMSD values.
	distUpper : float
		The upper limit for distance.
		Any pairwise distances larger than this value
		will be rejected from the final space of points
	distLower : float
		The lower distance limit for distance
		Any pairwise distance smaller than this value
		will be rejectet from the final space of points
	cutoffValue : float
		The RMSD threshold for which points are taken
		for the pairwise distance comparison.
		Note that the total number of points returned
		is significnalty influenced by this parameter
		and should be chosen carefully

	Returns
	-------
	tuple : np.ndarray of position coordinates
		two lists of [x,y,z] coordinates are returned associated
		with the inputs <densityMapA> and <densityMapB>.
		The returned coordinates are taken from the original grids
		and represent points that have another associated point in
		the other grid which is within the distance bounds and 
		have an RMSD below the cutoff threshold.
	"""
	idxCutoffA = np.where(densityMapA.density < cutoffValue)
	idxCutoffB = np.where(densityMapB.density < cutoffValue)
	posA = densityMapA.positions[idxCutoffA]
	posB = densityMapB.positions[idxCutoffB]
	dist = np.linalg.norm(posA[:,None] - posB, axis=2)
	idxA, idxB = np.where(np.logical_and(dist < distUpper, dist > distLower))
	idxSort = np.argsort(densityMapA.density[idxA] + densityMapB.density[idxB])
	return posA[idxA[idxSort]], posB[idxB[idxSort]]



def pcs_gradient_orthogonality_dot(metals, position):
	"""
	An experimental metric for calculating the likelihood
	of a particular nuclear position being well localised
	from multiple tensors.
	This algorithm calculates the normalised PCS graident
	vector arising at a given position from each tensor
	provided. It then calculates the pairwise dot product
	of each vector and averages their absolute value.

	Parameters
	----------
	metals : list of Metal objects
		a list of paramagpy metal objects which define the 
		tensors.
	position : numpy.ndarray of floats
		the [x,y,z] coordinate to calculated the
		orthogonality score

	Returns
	-------
	score : float
		the orthogonality score. This is necessarily
		between 0 and 1 for this metric. A lower
		score defies a more favourable orthogonality 
		between PCS gradient vectors.
		Note that a value of zero is not possible when
		more than 4 metals are provided.
	"""
	grads = np.array([metal.pcs_gradient(position) for metal in metals])
	norms = grads / np.linalg.norm(grads, axis=1)[:,None]
	dots = norms.dot(norms.T)
	mask = ~np.eye(dots.shape[0], dtype=bool)
	score = np.sum(np.abs(dots[mask])) / (len(metals)**2 - len(metals))
	return score


def pcs_gradient_orthogonality_cross(metals, position):
	"""
	An experimental metric for calculating the likelihood
	of a particular nuclear position being well localised
	from multiple tensors.
	The metric scores the orthogonality of the PCS
	gradient vectors and accounts for their magnitude.
	This algorithm calculates the PCS graident
	vector arising at a given position from each tensor
	provided. It then calculates the pairwise cross product
	of each vector and averages their norm (length).

	Parameters
	----------
	metals : list of Metal objects
		a list of paramagpy metal objects which define the 
		tensors.
	position : numpy.ndarray of floats
		the [x,y,z] coordinate to calculated the
		orthogonality score

	Returns
	-------
	score : float
		the orthogonality score. This is necessarily
		between greater than zero. A larger
		score defines a more favourable orthogonality 
		between PCS gradient vectors and also score larger
		PCS gradients more highly.
	"""
	grads = np.array([metal.pcs_gradient(position) for metal in metals])*1E-10
	vals = []
	for g1 in grads:
		for g2 in grads:
			cross = np.cross(g1, g2)
			vals.append(np.linalg.norm(cross)**0.5)

	score = np.sum(vals) / (len(metals)**2 - len(metals))
	return score





