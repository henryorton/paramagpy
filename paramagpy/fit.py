import numpy as np
from scipy.optimize import fmin_bfgs
from pprint import pprint
import warnings
from collections import OrderedDict

def unique_pairing(a, b):
	"""
	Bijectively map two integers to a single integer.
	The mapped space is minimum size.
	The input is symmetric.
	see `Bijective mapping f:ZxZ->N <https://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way/>`_.

	Parameters
	----------
	a : int
	b : int

	Returns
	-------
	c : int
		bijective symmetric mapping (a, b) | (b, a) -> c
	"""
	c = a * b + (abs(a - b) - 1)**2 // 4
	return c

def cantor_pairing(a, b):
	"""
	Map two integers to a single integer.
	The mapped space is minimum size.
	Ordering matters in this case.
	see `Bijective mapping f:ZxZ->N <https://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way/>`_.

	Parameters
	----------
	a : int
	b : int

	Returns
	-------
	c : int
		bijective mapping (a, b) -> c
	"""
	c = (a + b)*(a + b + 1)//2 + b
	return c

def clean_indices(indices):
	"""
	Uniquely map a list of integers to their smallest size.
	For example: [7,4,7,9,9,10,1] -> [4 2 4 0 0 1 3]

	Parameters
	----------
	indices : array-like integers
		a list of integers

	Returns
	-------
	new_indices : array-like integers
		the mapped integers with smallest size
	"""
	translation = {idx:i for i, idx in enumerate(set(indices))}
	new_indices = [translation[idx] for idx in indices]
	return np.array(new_indices)

def extract_pcs(data):
	"""
	Extract values required for PCS calculations

	Parameters
	----------
	data : list of lists
		A list with elements [Atom, value, error], where Atom is 
		an Atom object, value is the PCS value, and error is the uncertainty

	Returns
	-------
	tuple : (atom coordinates, PCS values, PCS errors, atom indices)
		all information required for PCS calculations
	"""
	atoms, values, errors = zip(*data)
	coords = np.array([i.position for i in atoms])
	values = np.array(values)
	errors = np.array(errors)
	if 0.0 in errors:
		errors = np.ones(len(errors))
		warnings.warn("0.0 value uncertainty. All values weighted evenly")
	idxs = clean_indices([i.serial_number for i in atoms])
	return (coords, values, errors, idxs)

def extract_pre(data):
	"""
	Extract values required for PRE calculations

	Parameters
	----------
	data : list of lists
		A list with elements [Atom, value, error], where Atom is 
		an Atom object, value is the PRE value, and error is the uncertainty

	Returns
	-------
	tuple : (atom coordinates, PRE values, PRE errors, atom indices)
		all information required for PRE calculations
	"""
	atoms, values, errors = zip(*data)
	coords = np.array([i.position for i in atoms])
	gammas = np.array([i.gamma for i in atoms])
	values = np.array(values)
	errors = np.array(errors)
	if 0.0 in errors:
		errors = np.ones(len(errors))
		warnings.warn("0.0 value uncertainty. All values weighted evenly")
	idxs = clean_indices([i.serial_number for i in atoms])
	return (coords, gammas, values, errors, idxs)

def extract_csa(data):
	"""
	Extract CSA tensors from atoms

	Parameters
	----------
	data : list of lists
		A list with elements [Atom, value, error], where Atom is 
		an Atom object, value is the PCS/RDC/PRE value, 
		and error is the uncertainty

	Returns
	-------
	csas : array of 3x3 arrays
		an array of each CSA tensor
	"""
	atoms, values, errors = zip(*data)
	return np.array([i.csa for i in atoms])

def extract_rdc(data):
	"""
	Extract values required for RDC calculations

	Parameters
	----------
	data : list of lists
		A list with elements [Atom1, Atom2, value, error], where Atom is 
		an Atom object, value is the RDC value, and error is the uncertainty

	Returns
	-------
	tuple : (inter-atomic vector, gamma values, RDC values, 
			RDC errors, atom indices)
		all information required for RDC calculations and fitting
	"""
	atoms1, atoms2, values, errors = zip(*data)
	vectors = [j.position - i.position for i, j in zip(atoms1, atoms2)]
	gammas = [i.gamma * j.gamma for i, j in zip(atoms1, atoms2)]
	if 0.0 in errors:
		errors = np.ones(len(errors))
		warnings.warn("0.0 value uncertainty. All values weighted evenly")
	idxs = clean_indices([unique_pairing(i.serial_number, 
		j.serial_number) for i, j in zip(atoms1, atoms2)])
	return map(np.array, [vectors, gammas, values, errors, idxs])

def extract_ccr(data):
	"""
	Extract values required for CCR calculations

	Parameters
	----------
	data : list of lists
		A list with elements [Atom1, Atom2, value, error], where Atom is 
		an Atom object, value is the CCR value, and error is the uncertainty

	Returns
	-------
	tuple : dict
		all information required for CCR calculations
	"""
	atoms, atomsPartner, values, errors = zip(*data)
	if 0.0 in errors:
		errors = np.ones(len(errors))
		warnings.warn("0.0 value uncertainty. All values weighted evenly")

	d = {
		'pos':np.array([i.position for i in atoms]),
		'gam':np.array([i.gamma for i in atoms]),
		'dst':np.array([j.dipole_shift_tensor(i.position) 
						for i, j in zip(atoms, atomsPartner)]),
		'val':np.array(values),
		'err':np.array(errors),
		'idx':clean_indices([cantor_pairing(i.serial_number, j.serial_number) 
							 for i, j in zip(atoms, atomsPartner)])
	}
	return d

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
	mat = np.array([np.bincount(idx, weights=col) for col in mat])*1E-24
	matinv = np.linalg.pinv(mat)
	sol = matinv.T.dot(pcs*1E-6)*1E-24
	calc = mat.T.dot(sol)
	return calc, sol


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
	dist = np.linalg.norm(pos, axis=1)
	x, y, z = pos.T
	a = x**2 - z**2
	b = y**2 - z**2
	c = 2 * x * y
	d = 2 * x * z
	e = 2 * y * z
	scale = 1./(4.*np.pi*dist**5)
	mat = scale * np.array([a,b,c,d,e,1E26/scale]) / errors
	mat = np.array([np.bincount(idx, weights=col) for col in mat])*1E-24
	matinv = np.linalg.pinv(mat)
	sol = matinv.T.dot(pcs*1E-6)*1E-24
	sol[-1] *= 1E26
	calc = mat.T.dot(sol)
	return calc, sol


def svd_gridsearch_fit_metal_from_pcs(metals, pcss, sumIndices=None,
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
	assert len(metals)==len(pcss)

	if origin is None:
		origin = metals[0].position*1E10

	if offsetShift:
		svd_func = svd_calc_metal_from_pcs_offset
	else:
		svd_func = svd_calc_metal_from_pcs

	posarrays = []
	pcsarrays = []
	errarrays = []
	idxarrays = []
	for pcs in pcss:
		posarray, pcsarray, errarray, idxarray = extract_pcs(pcs)
		posarrays.append(posarray)
		pcsarrays.append(pcsarray)
		errarrays.append(errarray)
		idxarrays.append(idxarray)
	sphere = sphere_grid(origin, radius, points)*1E-10

	if sumIndices is not None:
		idxarrays = sumIndices

	pcsarrays_eavg = []
	for pcsarray, idxarray, errarray in zip(pcsarrays, idxarrays, errarrays):
		pcsarrays_eavg.append(np.bincount(idxarray, weights=pcsarray/errarray))

	minscore = 1E50
	print("SVD search started in {} points".format(len(sphere)))
	tot = len(sphere)
	prog = 0.0
	for pos in sphere:
		if progress:
			prog += 1
			progress.set(prog/tot)
		score = 0.0
		sols = []
		zipped = zip(pcsarrays_eavg, posarrays, idxarrays, errarrays)
		for pcsarray_eavg, posarray, idxarray, errarray in zipped:
			coords = posarray - pos
			calculated, solution = svd_func(coords, 
				pcsarray_eavg, idxarray, errarray)
			sols.append(solution)
			score += np.sum((calculated - pcsarray_eavg)**2)
		if score<minscore:
			minscore = score
			minpos = pos
			minsols = sols

	minmetals = [m.copy() for m in metals]
	calc_pcss = []
	qfactors = []
	for pcsarray, posarray, idxarray, metal, sol in zip(pcsarrays, posarrays, 
		idxarrays, minmetals, minsols):
		metal.position = minpos
		if offsetShift:
			metal.upper_triang = sol[:-1]
			metal.shift = sol[-1]*1E6
		else:
			metal.upper_triang = sol
		calculated = metal.fast_pcs(posarray)
		calc_pcss.append(calculated)
		qfac = qfactor(pcsarray, calculated, idxarray)
		qfactors.append(qfac)

	return minmetals, calc_pcss, qfactors


def nlr_fit_metal_from_pcs(initMetals, pcss, 
	params=('x','y','z','ax','rh','a','b','g'), 
	sumIndices=None, userads=False, useracs=False, progress=None):
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
	metals : list of metals
		the metals fitted by NLR to the PCS data provided
	calc_pcss : list of lists of floats
		the calculated PCS values
	"""
	posarrays = []
	csaarrays = []
	pcsarrays = []
	errarrays = []
	idxarrays = []
	for pcs in pcss:
		posarray, pcsarray, errarray, idxarray = extract_pcs(pcs)
		posarrays.append(posarray)
		pcsarrays.append(pcsarray)
		errarrays.append(errarray)
		idxarrays.append(idxarray)
		if useracs:
			csas = extract_csa(pcs)
			csaarrays.append(csas)
		else:
			csaarrays.append(None)

	metals = [metal.copy() for metal in initMetals]
	pospars = [param for param in params if param in ['x','y','z']]
	otherpars = [param for param in params if param not in ['x','y','z']]

	if sumIndices is not None:
		idxarrays = sumIndices

	def cost(args):
		pos = args[:len(pospars)]
		allother = args[len(pospars):]
		for i, metal in enumerate(metals):
			other = allother[len(otherpars)*i:len(otherpars)*(i+1)]
			metal.set_params(zip(pospars, pos))
			metal.set_params(zip(otherpars, other))
		score = 0.0
		zipped = zip(metals, posarrays, csaarrays, pcsarrays, idxarrays, errarrays)
		for metal, posarray, csaarray, pcsarray, idxarray, errarray in zipped:
			calcpcs = metal.fast_pcs(posarray)
			if userads:
				calcpcs += metal.fast_rads(posarray)
			if useracs:
				calcpcs += metal.fast_racs(csaarray)
			diff = (calcpcs - pcsarray) / errarray
			selectiveSum = np.bincount(idxarray, weights=diff)
			score += np.sum(selectiveSum**2)
		return score

	startpars = metals[0].get_params(pospars)

	for metal in metals:
		pars = metal.get_params(otherpars)
		startpars += pars
	fmin_bfgs(cost, startpars, disp=False)
	calc_pcss = []
	qfactors = []
	zipped = zip(metals, posarrays, csaarrays, pcsarrays, idxarrays)
	for metal, posarray, csaarray, pcsarray, idxarray in zipped:
		metal.set_utr()
		calculated = metal.fast_pcs(posarray)
		if userads:
			calculated += metal.fast_rads(posarray)
		if useracs:
			calculated += metal.fast_racs(csaarray)
		calc_pcss.append(calculated)
		qfac = qfactor(pcsarray, calculated, idxarray)
		qfactors.append(qfac)

	if progress:
		progress.set(1.0)
	return metals, calc_pcss, qfactors


def pcs_fit_error_monte_carlo(initMetals, pcss, params, iterations,
	sumIndices=None, userads=False, useracs=False, progress=None):
	stds = [[] for metal in initMetals]
	atmarrays = []
	pcsarrays = []
	errarrays = []
	for data in pcss:
		atmarray, pcsarray, errarray = zip(*data)
		atmarrays.append(atmarray)
		pcsarrays.append(pcsarray)
		errarrays.append(errarray)

	for i in range(iterations):
		data = []
		for atm, pcs, err in zip(atmarrays, pcsarrays, errarrays):
			# noisey = pcs + np.random.normal(scale=err)
			noisey = pcs + (np.random.uniform(low=-1, high=-1, 
				size=len(err)))*err
			tmp = zip(atm, noisey, err)
			data.append(list(tmp))
		metals, _, _ = nlr_fit_metal_from_pcs(initMetals, data, params, 
			sumIndices, userads, useracs, progress=None)
		for metal, std in zip(metals, stds):
			std.append(metal.get_params(params))
		if progress:
			progress.set(float(i+1)/iterations)

	return [dict(zip(params, zip(*std))) for std in stds]


def pcs_fit_error_bootstrap(initMetals, pcss, params, iterations, 
	fraction_removed, sumIndices=None, userads=False, useracs=False, 
	progress=None):
	assert 0.0<fraction_removed<1.0
	stds = [[] for metal in initMetals]
	atmarrays = []
	pcsarrays = []
	errarrays = []
	sumarrays = []
	for data in pcss:
		atmarray, pcsarray, errarray = zip(*data)
		atmarrays.append(atmarray)
		pcsarrays.append(pcsarray)
		errarrays.append(errarray)
		sumarrays.append(clean_indices([a.serial_number for a in atmarray]))

	if sumIndices is None:
		sumIndices = sumarrays

	for i in range(iterations):
		datas_trunc = []
		sumIndices_trunc = []
		for atm, pcs, err, idx in zip(atmarrays, pcsarrays, errarrays, sumarrays):
			unique_idx = np.unique(idx)
			chosen_idx = np.random.choice(unique_idx, 
				int(len(unique_idx)*(1-fraction_removed)), replace=False)
			mask = np.isin(idx, chosen_idx)
			idxs_trunc = idx[mask]
			sumIndices_trunc.append(idxs_trunc)
			data_trunc = np.array(list(zip(atm, pcs, err)))[mask]
			datas_trunc.append(data_trunc.tolist())
		metals, _, _ = nlr_fit_metal_from_pcs(initMetals, datas_trunc, params, 
			sumIndices_trunc, userads, useracs, progress=None)
		for metal, std in zip(metals, stds):
			std.append(metal.get_params(params))
		if progress:
			progress.set(float(i+1)/iterations)

	return [dict(zip(params, zip(*std))) for std in stds]


def qfactor(experiment, calculated, sumIndices=None):
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
	experiment = np.array(experiment)
	calculated = np.array(calculated)

	if len(experiment)==0:
		return np.nan
	if sumIndices is None:
		idxs = np.arange(len(experiment))
	else:
		idxs = clean_indices(sumIndices)
	diff = experiment - calculated
	numer = np.sum(np.bincount(idxs, weights=diff)**2)
	denom = np.sum(np.bincount(idxs, weights=experiment)**2)
	return (numer/denom)**0.5


def nlr_fit_metal_from_pre(initMetals, pres, params=('x','y','z'), 
	sumIndices=None, rtypes=None, usesbm=True, usedsa=True, usecsa=False, 
	progress=None):
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
	if rtypes is None:
		rtypes = ['r2']*len(initMetals)

	posarrays = []
	csaarrays = []
	prearrays = []
	gamarrays = []
	idxarrays = []
	errarrays = []
	for pre in pres:
		posarray, gamarray, prearray, errarray, idxarray = extract_pre(pre)
		posarrays.append(posarray)
		gamarrays.append(gamarray)
		prearrays.append(prearray)
		idxarrays.append(idxarray)
		errarrays.append(errarray)
		if usecsa:
			csas = extract_csa(pre)
			csaarrays.append(csas)
		else:
			csaarrays.append(0.0)

	metals = [metal.copy() for metal in initMetals]
	pospars = [param for param in params if param in ['x','y','z']]
	otherpars = [param for param in params if param not in ['x','y','z']]

	if sumIndices is not None:
		idxarrays = sumIndices

	def cost(args):
		pos = args[:len(pospars)]
		allother = args[len(pospars):]
		for i, metal in enumerate(metals):
			other = allother[len(otherpars)*i:len(otherpars)*(i+1)]
			metal.set_params(zip(pospars, pos))
			metal.set_params(zip(otherpars, other))
		score = 0.0
		zipped = zip(metals, posarrays, gamarrays, csaarrays, 
						prearrays, idxarrays, errarrays, rtypes)
		for metal, posarr, gamarr, csaarr, prearr, idxarr, errarr, rtype in zipped:
			calcpre = metal.fast_pre(posarr, gamarr, rtype, 
				dsa=usedsa, sbm=usesbm, csaarray=csaarr)
			diff = (calcpre - prearr) / errarr
			selectiveSum = np.bincount(idxarr, weights=diff)
			score += np.sum(selectiveSum**2)
		return score

	startpars = metals[0].get_params(pospars)
	for metal in metals:
		pars = metal.get_params(otherpars)
		startpars += pars
	fmin_bfgs(cost, startpars, disp=False)
	calc_pres = []
	qfactors = []
	zipped = zip(metals, posarrays, gamarrays, csaarrays, 
						prearrays, idxarrays, errarrays, rtypes)
	for metal, posarr, gamarr, csaarr, prearr, idxarr, errarr, rtype in zipped:
		metal.set_utr()
		calculated = metal.fast_pre(posarr, gamarr, rtype, 
			dsa=usedsa, sbm=usesbm, csaarray=csaarr)
		calc_pres.append(calculated)
		qfac = qfactor(prearr, calculated, idxarr)
		qfactors.append(qfac)

	if progress:
		progress.set(1.0)
	return metals, calc_pres, qfactors

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


def svd_fit_metal_from_rdc(metal, rdc, sumIndices=None):
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
	vecarray, gamarray, rdcarray, errarray, idxarray = extract_rdc(rdc)
	if sumIndices is None:
		sumIndices = idxarray
	pfarray = -3*(metal.MU0 * gamarray * metal.HBAR) / (8 * np.pi**2)
	rdc_parameterised = np.bincount(idxarray, 
		weights=rdcarray / (pfarray * errarray))
	fitMetal = metal.copy()
	_, sol = svd_calc_metal_from_rdc(vecarray, rdc_parameterised, 
		idxarray, errarray)
	fitMetal.upper_triang_alignment = sol
	calculated = fitMetal.fast_rdc(vecarray, gamarray)
	qfac = qfactor(rdcarray, calculated, idxarray)
	return fitMetal, calculated, qfac


def nlr_fit_metal_from_rdc(metal, rdc, params=('ax','rh','a','b','g'), 
	sumIndices=None, progress=None):
	"""
	Fit deltaChi tensor to RDC values using non-linear regression.

	Parameters
	----------
	metal : Metal object
		the starting metal for fitting
	rdc : the RDC dataset
		each RDC dataset has structure [Atom, value, error], where Atom is 
		an Atom object, value is the PCS/RDC/PRE value
		and error is the uncertainty
	params : list of str, optional
		the parameters to be fit. 
		this defaults to ['ax','rh','a','b','g']
	sumIndices : array of ints, optional
		the list contains an index assigned to each atom. 
		Common indices determine summation between models 
		for ensemble averaging.
		If None, defaults to atom serial number to determine summation 
		between models.
	progress : object, optional
		to keep track of the calculation, progress.set(x) is called each
		iteration and varies from 0.0 -> 1.0 when the calculation is complete.

	Returns
	-------
	fitMetal : Metal object
		the fitted metal by NLR to the RDC data provided
	calculated : array of floats
		the calculated RDC values
	qfac : float
		the qfactor judging the fit quality
	"""
	vecarray, gamarray, rdcarray, errarray, idxarray = extract_rdc(rdc)
	if sumIndices is None:
		sumIndices = idxarray
	fitMetal = metal.copy()

	def cost(args):
		fitMetal.set_params(zip(params, args))
		calrdc = fitMetal.fast_rdc(vecarray, gamarray)
		diff = (calrdc - rdcarray) / errarray
		selectiveSum = np.bincount(idxarray, weights=diff)
		score = np.sum(selectiveSum**2)
		return score

	startpars = fitMetal.get_params(params)
	fmin_bfgs(cost, startpars, disp=False)
	fitMetal.set_utr()
	calculated = fitMetal.fast_rdc(vecarray, gamarray)
	qfac = qfactor(rdcarray, calculated, idxarray)

	if progress:
		progress.set(1.0)

	return fitMetal, calculated, qfac


def nlr_fit_metal_from_ccr(initMetals, ccrs, params=('x','y','z'), 
	sumIndices=None, progress=None):
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
	ccrs : list of CCR datasets
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
	datas = []
	for metal, ccr in zip(initMetals, ccrs):
		d = extract_ccr(ccr)
		d['met'] = metal.copy()
		datas.append(d)

	params = set(params)
	pospars = tuple(params & set(['x','y','z']))
	otherpars = tuple(params - set(['x','y','z']))

	startpars = initMetals[0].get_params(pospars)
	for i, d in enumerate(datas):
		d['pospars'] = slice(0, len(pospars))
		d['othpars'] = slice(len(pospars) + i*len(otherpars), 
						  len(pospars) + (i+1)*len(otherpars))
		startpars += d['met'].get_params(otherpars)

	if sumIndices is not None:
		idxarrays = sumIndices

	def cost(args):
		score = 0.0
		for d in datas:
			d['met'].set_params(zip(pospars, args[d['pospars']]))
			d['met'].set_params(zip(otherpars, args[d['othpars']]))
			d['cal'] = d['met'].fast_ccr(d['pos'], d['gam'], d['dst'])
			diff = (d['cal'] - d['val']) / d['err']
			selectiveSum = np.bincount(d['idx'], weights=diff)
			score += np.sum(selectiveSum**2)

		return score

	fmin_bfgs(cost, startpars, disp=False)
	fitmetals = []
	calc_ccrs = []
	qfactors = []
	for d in datas:
		d['met'].set_utr()
		fitmetals.append(d['met'])
		calc_ccrs.append(d['cal'])
		qfac = qfactor(d['val'], d['cal'], d['idx'])
		qfactors.append(qfac)

	if progress:
		progress.set(1.0)
	return fitmetals, calc_ccrs, qfactors





def ensemble_average(atoms, *values):
	if type(atoms[0]) in (list, tuple):
		dtype = 'RDC'
	else:
		dtype = 'other'

	d = {}
	if dtype=='RDC':
		for x in zip(atoms, *values):
			key = unique_pairing(x[0][0].serial_number, x[0][1].serial_number)
			if key not in d:
				d[key] = []
			d[key].append(x[1:])

	else:
		for x in zip(atoms, *values):
			key = x[0].serial_number
			if key not in d:
				d[key] = []
			d[key].append(x[1:])

	out = []
	for x in d.values():
		vals = [sum(i)/float(len(i)) for i in zip(*x)]
		out.append(vals)

	return list(zip(*out))













