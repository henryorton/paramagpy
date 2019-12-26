import numpy as np
from scipy.optimize import fmin_bfgs
from pprint import pprint
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
	Extract values required for PCS calculations

	Parameters
	----------
	data : 

	Returns
	-------
	
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
	Extract values required for PCS calculations

	Parameters
	----------
	data : 

	Returns
	-------
	
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
	Extract values required for PCS calculations

	Parameters
	----------
	data : 

	Returns
	-------
	
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
		the PCS values in ppm
	idx : array of ints
		an index assigned to each atom. Common indices determine summation
		between models for ensemble averaging.
	errors : array of floats
		the standard deviation representing experimental uncertainty
		in the measured value

	Returns
	-------
	tuple : (calc, sol)
		calc are the calculated PCS values from the fitted tensor
		sol is the solution to the linearised PCS equation and 
		consists of the tensor matrix elements
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
		the PCS values in ppm
	idx : array of ints
		an index assigned to each atom. Common indices determine summation
		between models for ensemble averaging.
	errors : array of floats
		the standard deviation representing experimental uncertainty
		in the measured value

	Returns
	-------
	tuple : (calc, sol)
		calc are the calculated PCS values from the fitted tensor
		sol is the solution to the linearised PCS equation and 
		consists of the tensor matrix elements and offset
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
	Ensemble averaging is determined by atom number.

	Parameters
	----------
	metals : list of Metal objects
		a list of metals used as starting points for fitting. 
		a list must always be provided, but may also contain 
		only one element. If multiple metals are provided, each metal
		is fitted to their respective PCS dataset by index, but all are 
		fitted to a common position.
	pcss : list of PCS datasets
		each PCS dataset must correspond to an associated metal for fitting.
		each PCS dataset has structure [Atom, value, error], where Atom is 
		an Atom object, value is the PCS/RDC/PRE value
		and error is the uncertainty
	sumIndices : list of arrays of ints, optional
		each index list must correspond to an associated pcs dataset.
		each index list contains an index assigned to each atom. 
		Common indices determine summation between models 
		for ensemble averaging.
		If None, defaults to atom serial number to determine summation 
		between models.
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
	minmetals : list of metals
		the metals fitted by SVD to the PCS data provided
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
		fitMetals.append(mAvg)

	for m, data in zip(fitMetals, dataArrays):
		d = extract_atom_data(data, csa=False, separateModels=False)[0]
		data['cal'] = m.fast_pcs(d['pos'])

	if progress:
		progress.set(1.0)

	return fitMetals, dataArrays


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
		is fitted to their respective PCS dataset by index, but all are 
		fitted to a common position.
	pcss : list of PCS datasets
		each PCS dataset must correspond to an associated metal for fitting.
		each PCS dataset has structure [Atom, value, error], where Atom is 
		an Atom object, value is the PCS/RDC/PRE value
		and error is the uncertainty
	params : list of str
		the parameters to be fit. 
		For example ['x','y','z','ax','rh','a','b','g','shift']
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
	metals : list of metals
		the metals fitted by NLR to the PCS data provided
	calc_pcss : list of lists of floats
		the calculated PCS values
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
		fitMetals.append(mAvg)

	for m, data in zip(fitMetals, dataArrays):
		d = extract_atom_data(data, csa=useracs, separateModels=False)[0]
		data['cal'] = m.fast_pcs(d['pos'])
		if userads:
			data['cal'] += m.fast_rads(d['pos'])
		if useracs:
			data['cal'] += m.fast_racs(d['xsa'])

	if progress:
		progress.set(1.0)

	return fitMetals, dataArrays


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


def pcs_fit_error_monte_carlo(initMetals, dataArrays, iterations,
	params=('x','y','z','ax','rh','a','b','g'),
	ensembleAverage=False, userads=False, useracs=False, progress=None):
	"""
	Analyse uncertainty of PCS fit by Monte-Carlo simulation
	This repeatedly adds noise to experimental PCS data and fits the tensor.
	The standard deviation of the fitted parameters across each iteration
	is then reported.

	Parameters
	----------
	initMetals : list of Metal objects
		a list of metals used as starting points for fitting. 
		a list must always be provided, but may also contain 
		only one element. If multiple metals are provided, each metal
		is fitted to their respective PCS dataset by index, but all are 
		fitted to a common position.
	pcss : list of PCS datasets
		each PCS dataset must correspond to an associated metal for fitting.
		each PCS dataset has structure [Atom, value, error], where Atom is 
		an Atom object, value is the PCS/RDC/PRE value
		and error is the uncertainty
	iterations : int
		the number of Monte Carlo iterations to perform
	params : list of str
		the parameters to be fit. 
		For example ['x','y','z','ax','rh','a','b','g','shift']
	sumIndices : list of arrays of ints, optional
		each index list must correspond to an associated pcs dataset.
		each index list contains an index assigned to each atom. 
		Common indices determine summation between models 
		for ensemble averaging.
		If None, defaults to atom serial number to determine summation 
		between models.
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
	sample_metals : list of list of metals
		the metals fitted by NLR to the PCS data with noise at each iteration
	std_metals : list of metals
		the standard deviation in fitted parameters over all iterations of the
		Monte Carlo simulation.
		These are stored within the metal object. All unfitted parameters 
		are zero.
	"""

	sampleMetals = []
	for i in range(iterations):
		newInitMetals = []
		newDataArrays = []
		for m, d in zip(initMetals, dataArrays):
			newInitMetals.append(m.copy())
			d['exp'] += d['err'] * np.random.uniform(low=-1, high=1, size=len(d))
			newDataArrays.append(d)

		metals, _ = nlr_fit_metal_from_pcs(newInitMetals, newDataArrays, params=params, 
			ensembleAverage=ensembleAverage, userads=userads, useracs=useracs)

		sampleMetals.append(metals)

		if progress:
			progress.set(float(i+1)/iterations)

	sampleMetals = list(zip(*sampleMetals))
	stdMetals = []
	for metals in sampleMetals:
		stdMetals.append(metal_standard_deviation(metals, params))

	return sampleMetals, stdMetals


def pcs_fit_error_bootstrap(initMetals, dataArrays, iterations, fraction,
	params=('x','y','z','ax','rh','a','b','g'), 
	ensembleAverage=False, userads=False, useracs=False, progress=None):
	"""
	Analyse uncertainty of PCS fit by Bootstrap methods.
	This repeats the tensor fitting, but each time samples a fraction
	of the data randomly. The standard deviation in fitted parameters
	over each iteration is then reported.

	Parameters
	----------
	initMetals : list of Metal objects
		a list of metals used as starting points for fitting. 
		a list must always be provided, but may also contain 
		only one element. If multiple metals are provided, each metal
		is fitted to their respective PCS dataset by index, but all are 
		fitted to a common position.
	pcss : list of PCS datasets
		each PCS dataset must correspond to an associated metal for fitting.
		each PCS dataset has structure [Atom, value, error], where Atom is 
		an Atom object, value is the PCS/RDC/PRE value
		and error is the uncertainty
	iterations : int
		the number of Monte Carlo iterations to perform
	fraction : float
		must be between 0 and 1
		the proportion of data to be sample for fitting with each iteration
		of the bootstrap method.
	params : list of str
		the parameters to be fit. 
		For example ['x','y','z','ax','rh','a','b','g','shift']
	sumIndices : list of arrays of ints, optional
		each index list must correspond to an associated pcs dataset.
		each index list contains an index assigned to each atom. 
		Common indices determine summation between models 
		for ensemble averaging.
		If None, defaults to atom serial number to determine summation 
		between models.
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
	sample_metals : list of list of metals
		the metals fitted by NLR to the PCS data with noise at each iteration
	std_metals : list of metals
		the standard deviation in fitted parameters over all iterations of the
		Monte Carlo simulation.
		These are stored within the metal object. All unfitted parameters 
		are zero.
	"""
	if not (0.0<fraction<1.0):
		raise ValueError("The bootstrap sample fraction must be between 0 and 1")

	sampleMetals = []
	for i in range(iterations):
		newInitMetals = []
		newDataArrays = []
		for m, d in zip(initMetals, dataArrays):
			newInitMetals.append(m.copy())
			d = np.random.choice(d, int(fraction*len(d)), replace=False)
			newDataArrays.append(d)

		metals, _ = nlr_fit_metal_from_pcs(newInitMetals, newDataArrays, params=params, 
			ensembleAverage=ensembleAverage, userads=userads, useracs=useracs)

		sampleMetals.append(metals)

		if progress:
			progress.set(float(i+1)/iterations)

	sampleMetals = list(zip(*sampleMetals))
	stdMetals = []
	for metals in sampleMetals:
		stdMetals.append(metal_standard_deviation(metals, params))

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
	experiment : list of floats
		the experimental values
	calculated : list of floats
		the corresponding caluclated values from the fitted model
	sumIndices : list ints, optional
		Common indices determine summation between models 
		for ensemble averaging.
		If None, no ensemble averaging is conducted

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
		denom = np.sum(np.bincount(dataArray['idx'], weights=tmp)**2)
	return (numer/denom)**0.5


def nlr_fit_metal_from_pre(initMetals, dataArrays, rtypes, params=('x','y','z'), 
	usesbm=True, usegsbm=False, usedsa=True, 
	usecsa=False, ensembleAverage=False, progress=None):
	"""
	Fit deltaChi tensor to PRE values using non-linear regression.

	Parameters
	----------
	initMetals : list of Metal objects
		a list of metals used as starting points for fitting. 
		a list must always be provided, but may also contain 
		only one element. If multiple metals are provided, each metal
		is fitted to their respective PRE dataset by index, but all are 
		fitted to a common position.
	pres : list of PRE datasets
		each PRE dataset must correspond to an associated metal for fitting.
		each PRE dataset has structure [Atom, value, error], where Atom is 
		an Atom object, value is the PCS/RDC/PRE value
		and error is the uncertainty
	params : list of str
		the parameters to be fit. 
		For example ['x','y','z','ax','rh','a','b','g','iso','taur','t1e']
	sumIndices : list of arrays of ints, optional
		each index list must correspond to an associated pcs dataset.
		each index list contains an index assigned to each atom. 
		Common indices determine summation between models 
		for ensemble averaging.
		If None, defaults to atom serial number to determine summation 
		between models.
	rtypes : list of str, optional
		the relaxtion type, either 'r1' or 'r2'. A list must be provided with
		each element corresponding to an associated dataset.
		Defaults to 'r2' for all datasets of None is specified.
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
	progress : object, optional
		to keep track of the calculation, progress.set(x) is called each
		iteration and varies from 0.0 -> 1.0 when the calculation is complete.

	Returns
	-------
	metals : list of metals
		the metals fitted by NLR to the PRE data provided
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
		mAvg.average(metalAvg)
		mAvg.set_utr()
		mAvg.metalAvg = metalAvg
		fitMetals.append(mAvg)

	for m, data in zip(fitMetals, dataArrays):
		d = extract_atom_data(data, csa=usecsa, separateModels=False)[0]
		data['cal'] = m.fast_pre(d['pos'], d['gam'], m.par['rtp'], 
				dsa=usedsa, sbm=usesbm, gsbm=usegsbm, csaarray=d['xsa'])

	if progress:
		progress.set(1.0)

	return fitMetals, dataArrays






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


def svd_fit_metal_from_rdc(initMetal, dataArray, ensembleAverage=False):
	"""
	Fit deltaChi tensor to RDC values using SVD algorithm.
	Note this is a weighted SVD calculation which takes into account
	experimental errors.
	Ensemble averaging defaults to atom numbers.

	Parameters
	----------
	metal : Metal object
		the starting metal for fitting
	rdc : the RDC dataset
		each RDC dataset has structure [(Atom1, Atom2), value, error], 
		where Atom is an Atom object, value is the RDC value
		and error is the uncertainty
	sumIndices : array of ints, optional
		the list contains an index assigned to each atom. 
		Common indices determine summation between models 
		for ensemble averaging.
		If None, defaults to atom serial number to determine summation 
		between models.

	Returns
	-------
	fitMetal : Metal object
		the fitted metal by NLR to the RDC data provided
	calculated : array of floats
		the calculated RDC values
	qfac : float
		the qfactor judging the fit quality
	"""
	datas = {}
	metalAvgs = []
	if ensembleAverage:
		d = extract_rdc_data(dataArray, separateModels=False)[0]
		m = initMetal.copy()
		metalAvgs.append(m)
		datas[0] = m, d
	else:
		for d in extract_rdc_data(dataArray, separateModels=True):
			mdl = d['mdl'][0]
			m = initMetal.copy()
			metalAvgs.append(m)
			datas[mdl] = m, d

	for mdl in datas:
		m, d = datas[mdl]
		pfarray = -3*(m.MU0 * d['gam'] * m.HBAR) / (8 * np.pi**2)
		rdc_parameterised = np.bincount(d['idx'], weights=d['exp'] / (pfarray * d['err']))
		calculated, solution = svd_calc_metal_from_rdc(d['pos'], rdc_parameterised, d['idx'], d['err'])
		m.upper_triang_alignment = solution

	fitMetal = metalAvgs[0].copy()
	fitMetal.average(metalAvgs)
	fitMetal.set_utr()

	d = extract_rdc_data(dataArray, separateModels=False)[0]
	dataArray['cal'] = m.fast_rdc(d['pos'], d['gam'])

	return fitMetal, dataArray


def nlr_fit_metal_from_ccr(initMetals, dataArrays, params=('x','y','z'), 
	ensembleAverage=False, progress=None):
	"""
	Fit deltaChi tensor to CCR values using non-linear regression.

	Parameters
	----------
	initMetals : list of Metal objects
		a list of metals used as starting points for fitting. 
		a list must always be provided, but may also contain 
		only one element. If multiple metals are provided, each metal
		is fitted to their respective CCR dataset by index, but all are 
		fitted to a common position.
	dataArrays : list of CCR datasets
		each CCR dataset must correspond to an associated metal for fitting.
		each CCR dataset has structure [Atom, value, error], where Atom is 
		an Atom object, value is the PCS/RDC/PRE/CCR value
		and error is the uncertainty
	params : list of str
		the parameters to be fit. 
		For example ['x','y','z','ax','rh','a','b','g','shift']
		This defaults to ['x','y','z']
	sumIndices : list of arrays of ints, optional
		each index list must correspond to an associated ccr dataset.
		each index list contains an index assigned to each atom. 
		Common indices determine summation between models 
		for ensemble averaging.
		If None, defaults to atom serial number to determine summation 
		between models.
	progress : object, optional
		to keep track of the calculation, progress.set(x) is called each
		iteration and varies from 0.0 -> 1.0 when the calculation is complete.

	Returns
	-------
	metals : list of metals
		the metals fitted by NLR to the CCR data provided
	calc_ccrs : list of lists of floats
		the calculated CCR values
	qfactors : list
	"""
	if len(initMetals)!=len(dataArrays):
		raise ValueError("initMetals and dataArrays must have same length")

	datas = {}
	metalAvgs = []
	for metal, dataArray in zip(initMetals, dataArrays):
		metalAvg = []
		if ensembleAverage:
			tmp = fit.extract_ccr_data(dataArray, separateModels=False)[0]
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

	startpars = initMetals[0].get_params(pospars)
	for mdl in datas:
		data = datas[mdl]
		for i, (m, d) in enumerate(data):
			m.par['pos'] = slice(0, len(pospars))
			m.par['oth'] = slice(len(pospars) + i*len(otherpars), 
								  len(pospars) + (i+1)*len(otherpars))
			startpars += m.get_params(otherpars)

	def cost(args):
		score = 0.0
		for mdl in datas:
			data = datas[mdl]
			for m, d in data:
				m.set_params(zip(pospars, args[m.par['pos']]))
				m.set_params(zip(otherpars, args[m.par['oth']]))
				d['cal'] = m.fast_ccr(d['pos'], d['gam'], d['xsa'])
				diff = (d['cal'] - d['exp']) / d['err']
				selectiveSum = np.bincount(d['idx'], weights=diff)
				score += np.sum(selectiveSum**2)

		return score

	fmin_bfgs(cost, startpars, disp=False)

	fitMetals = []
	for metalAvg in metalAvgs:
		mAvg = metalAvg[0].copy()
		mAvg.average(metalAvg)
		mAvg.set_utr()
		fitMetals.append(mAvg)

	if progress:
		progress.set(1.0)

	return fitMetals, dataArrays

















