import numpy as np
import struct
from collections import OrderedDict
import ntpath

def euler_to_matrix(eulers):
	"""
	Calculate a rotation matrix from euler angles using ZYZ convention

	Parameters
	----------
	eulers : array of floats
		the euler angles [alpha,beta,gamma] in radians 
		by ZYZ convention.

	Returns
	-------
	matrix : 3x3 numpy ndarray
		the rotation matrix

	Examples
	--------
	>>> eulers = np.array([0.5,1.2,0.8])
	>>> euler_to_matrix(eulers)
	array([[-0.1223669 , -0.5621374 ,  0.81794125],
	       [ 0.75057357,  0.486796  ,  0.44684334],
	       [-0.64935788,  0.66860392,  0.36235775]])
	"""
	c1, c2, c3 = np.cos(eulers)
	s1, s2, s3 = np.sin(eulers)
	M = np.identity(3)
	M[0,0] = c1*c2*c3 - s1*s3
	M[0,1] = -c3*s1 -c1*c2*s3
	M[0,2] = c1*s2
	M[1,0] = c1*s3 + c2*c3*s1
	M[1,1] = c1*c3 - c2*s1*s3
	M[1,2] = s1*s2
	M[2,0] = -c3*s2
	M[2,1] = s2*s3
	M[2,2] = c2
	return M


def matrix_to_euler(M):
	"""
	Calculate Euler angles from a rotation matrix using ZYZ convention

	Parameters
	----------
	M : 3x3 numpy ndarray
		a rotation matrix

	Returns
	-------
	eulers : array of floats
		the euler angles [alpha,beta,gamma] in radians
		by ZYZ convention

	Examples
	--------
	>>> matrix = array([[-0.1223669 , -0.5621374 ,  0.81794125],
	                    [ 0.75057357,  0.486796  ,  0.44684334],
	                    [-0.64935788,  0.66860392,  0.36235775]])
	>>> matrix_to_euler(matrix)
	np.array([0.5,1.2,0.8])
	"""
	if M[2,2]<1.0:
		if M[2,2] > -1.0:
			alpha = np.arctan2(M[1,2],M[0,2])
			beta  = np.arccos(M[2,2])
			gamma = np.arctan2(M[2,1],-M[2,0])
		else:
			alpha = -np.arctan2(M[1,0],M[1,1])
			beta  = np.pi
			gamma = 0.0
	else:
		alpha = np.arctan2(M[1,0],M[1,1])
		beta =  0.0
		gamma = 0.0
	return np.array([alpha,beta,gamma])


def unique_eulers(eulers):
	"""
	Calculate Euler angles in unique tensor representation.

	Given general Euler angles by ZYZ convention, this function accounts for 
	the symmetry of	a second rank symmetric tensor to map all angles within 
	the range [0, pi].

	Parameters
	----------
	eulers : array of float
		the three Euler angles in radians

	Returns
	-------
	eulers_utr : array of floats
		the euler angles [alpha,beta,gamma] in radians
		by ZYZ convention

	Examples
	--------
	>>> eulers = np.array([-5.2,10.3,0.1])
	>>> unique_eulers(eulers)
	np.array([1.08318531 0.87522204 3.04159265])
	"""
	def normalise(x):
		mod = x % (2*np.pi)
		if mod < np.pi:
			return mod
		else:
			return mod - 2*np.pi

	a, b, g = map(normalise, eulers)

	if   a>=0 and b>=0 and g>=0:
		alpha = a
		beta  = b
		gamma = g
	elif a<0 and b>=0 and g>=0:
		alpha = a + np.pi
		beta  = np.pi - b
		gamma = np.pi - g
	elif a>=0 and b<0 and g>=0:
		alpha = a
		beta  = b + np.pi
		gamma = np.pi - g
	elif a>=0 and b>=0 and g<0:
		alpha = a
		beta  = b
		gamma = g + np.pi
	elif a<0 and b<0 and g>=0:
		alpha = a + np.pi
		beta  = -b
		gamma = g
	elif a<0 and b>=0 and g<0:
		alpha = a + np.pi
		beta  = np.pi - b
		gamma = -g
	elif a>=0 and b<0 and g<0:
		alpha = a
		beta  = b + np.pi
		gamma = -g
	elif a<0 and b<0 and g<0:
		alpha = a + np.pi
		beta  = -b
		gamma = g + np.pi
	else:
		alpha = a
		beta  = b
		gamma = g

	eulers_utr = np.array([alpha, beta, gamma])
	return eulers_utr


def anisotropy_to_eigenvalues(axial_rhombic):
	"""
	Calculate [dx,dy,dz] eigenvalues from axial and rhombic
	tensor anisotropies (axial and rhombic parameters).

	Calculations assume traceless tensor.

	Parameters
	----------
	axial_rhombic : array of floats
		two values defining the axial and rhombic anisotropy 
		of the tensor respectively

	Returns
	-------
	eigenvalues : array of floats
		the three eigenvalues defining the magnitude of the priciple axes

	Examples
	--------
	>>> ax = 10.5
	>>> rh = 3.0
	>>> anisotropy_to_eigenvalues([ax, rh])
	[-2. -5.  7.]
	"""
	axial, rhombic = axial_rhombic
	dx =  rhombic/2. - axial/3.
	dy = -rhombic/2. - axial/3.
	dz = (axial*2.)/3.
	# return np.array(sorted([dx,dy,dz], key = lambda v: abs(v)))
	return np.array([dx,dy,dz])


def eigenvalues_to_anisotropy(eigenvalues):
	"""
	Calculate axial and rhombic tensor anisotropies from 
	eigenvalues dx,dy,dz

	Parameters
	----------
	eigenvalues : array of floats
		the three eigenvalues of the tensor.
		These are the principle axis magnitudes

	Returns
	-------
	axial_rhombic : tuple of floats
		the tensor axial/rhombic anisotropies respectively

	Examples
	--------
	>>> eigenvalues = [-2.0, -5.0, 7.0]
	>>> eigenvalues_to_anisotropy(eigenvalues)
	[10.5  3. ]
	"""
	dx, dy, dz = eigenvalues
	axial = dz - (dx + dy) / 2.
	rhombic = dx - dy
	return np.array([axial, rhombic])



class Metal(object):
	"""
	An object for paramagnetic chi tensors and delta-chi tensors.
	This class has basic attributes that specify position, 
	axiality/rhombicity, isotropy and euler angles.
	It also has methods for calculating PCS, RDC, PRE and CCR values.
	"""

	# Gyromagnetic ratio of an electron and other constants
	MU0 = 4*np.pi*1E-7
	MUB = 9.274E-24
	K = 1.381E-23
	HBAR = 1.0546E-34
	GAMMA = 1.760859644E11

	# Stored values get scaled by this amount for the fitting algorithm
	# This ensures against floating point errors during least-squares
	fit_scaling = {
		'x': 1E10,
		'y': 1E10,
		'z': 1E10,
		'a': 180.0/np.pi,
		'b': 180.0/np.pi,
		'g': 180.0/np.pi,
		'ax': 1E32,
		'rh': 1E32,
		'iso': 1E32,
		't1e': 1E12,
		'taur': 1E9,
		'mueff': 1./MUB,
		'gax': 1E60,
		'grh': 1E60
	}

	# J, g, T1e values for lanthanide series
	lanth_lib = OrderedDict([
		('Zero',(0.   , 0.    , 0.       )),
		('Ce', ( 5./2., 6./7. , 0.133E-12)),
		('Pr', ( 4./1., 4./5. , 0.054E-12)),
		('Nd', ( 9./2., 8./11., 0.210E-12)),
		('Pm', ( 4./1., 3./5. , 0.0      )),
		('Sm', ( 5./2., 2./7. , 0.074E-12)),
		('Eu', ( 2./1., 3./2. , 0.015E-12)),
		('Gd', ( 7./2., 2./1. , 1E-7    )),
		('Tb', ( 6./1., 3./2. , 0.251E-12)),
		('Dy', (15./2., 4./3. , 0.240E-12)),
		('Ho', ( 8./1., 5./4. , 0.209E-12)),
		('Er', (15./2., 6./5. , 0.189E-12)),
		('Tm', ( 6./1., 7./6. , 0.268E-12)),
		('Yb', ( 7./2., 8./7. , 0.157E-12))]
	)

	# Template anisotropies [axial, rhombic] from Bertini
	lanth_axrh = OrderedDict([
		('Zero',( 0.0,   0.0 )),
		('Ce', (  2.08,  0.71)),
		('Pr', (  3.40,  2.11)),
		('Nd', (  1.74,  0.46)),
		('Pm', (  0.0,   0.0 )),
		('Sm', (  0.19,  0.08)),
		('Eu', ( -2.34, -1.63)),
		('Gd', (  0.0,   0.0 )),
		('Tb', ( 42.1,  11.2 )),
		('Dy', ( 34.7,  20.3 )),
		('Ho', ( 18.5,   5.79)),
		('Er', (-11.6,  -8.58)),
		('Tm', (-21.9,  -20.1)),
		('Yb', ( -8.26, -5.84))]
	)

	fundamental_attributes = (
		'position', 
		'eulers', 
		'axrh', 
		'mueff', 
		'g_axrh', 
		't1e',
		'shift', 
		'temperature', 
		'B0', 
		'taur'
		)

	# Indices defining 5 unique elements of 3x3 tensor anisotropy
	upper_coords = ((0,1,0,0,1),(0,1,1,2,2))
	lower_coords = ((0,1,1,2,2),(0,1,0,0,1))

	def __init__(self, position=(0,0,0), eulers=(0,0,0), 
		axrh=(0,0), mueff=0.0, 
		g_axrh=(0,0), t1e=0.0,
		shift=0.0, temperature=298.15, 
		B0=18.79, taur=0.0):
		"""
		Instantiate ChiTensor object

		Parameters
		----------
		position : array of floats, optional
			the (x,y,z) position in meters. Default is (0,0,0)
			stored as a np.matrix object.
		eulers : array of floats, optional
			the euler angles [alpha,beta,gamma] in radians 
			by ZYZ convention. Defualt is (0,0,0)
		axrh : array of floats, optional
			the axial and rhombic values defining the magnetic susceptibility
			anisotropy
		g_axrh : array of floats, optional
			the axial and rhombic values defining the power spectral density
			tensor
		mueff : float
			the effective magnetic moment in units of A.m^2
		shift : float
			a bulk shift value applied to all PCS calculations.
			This is a correction parameter that may arise due to an offset
			between diamagnetic and paramagnetic PCS datasets.
		temperature : float
			the temperature in Kelvin
		t1e : float
			the longitudinal electronic relaxation time
		B0 : float
			the magnetic field in Telsa
		taur : float
			the rotational correlation time in seconds
		"""
		self.position = np.array(position, dtype=float)
		self.eulers = np.array(eulers, dtype=float)
		self.axrh = np.array(axrh, dtype=float)
		self.g_axrh = np.array(g_axrh, dtype=float)
		self.mueff = mueff
		self.shift = shift
		self.temperature = temperature
		self.t1e = t1e
		self.B0 = B0
		self.taur = taur

		# A dictionary of volatile parameters used during fitting
		self.par = {}

	def __str__(self):
		return str(self.tensor)

	def __abs__(self):
		return sum(self.eigenvalues)/3.

	@property
	def tauc(self):
		"""
		The effective rotational correlation time.

		This is calculated by combining the rotational correaltion time
		and the electronic relaxation time:

		.. math::
			\\tau_c = \\frac{1}{\\frac{1}{\\tau_r}+\\frac{1}{T_{1e}}}
		"""
		if self.taur==0.0 or self.t1e==0.0:
			raise ValueError("You need to set both taur and t1e attributes")
		return 1./(1./self.taur + 1./self.t1e)

	def copy(self):
		"""
		Copy the current Metal object to a new instance

		Returns
		-------
		new_tensor : Metal object
			a new Metal instance with the same parameters
		"""
		return self.__class__(
			position=tuple(self.position), 
			eulers=tuple(self.eulers),
			axrh=tuple(self.axrh),
			g_axrh=tuple(self.g_axrh),
			mueff=self.mueff, 
			shift=self.shift, 
			temperature=self.temperature, 
			t1e=self.t1e,
			B0=self.B0, 
			taur=self.taur)

	def set_lanthanide(self, lanthanide, set_dchi=True):
		"""
		Set the anisotropy, isotropy and T1e parameters from
		literature values

		Parameters
		----------
		lanthanide : str
			one of 'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb',
			'Dy','Ho','Er','Tm','Yb'
		set_dichi : bool (optional)
			if True (default), the tensor anisotropy is set.
			Otherwise only the isotropy and T1e values are set 
		"""
		J, g, t1e = self.lanth_lib[lanthanide]
		self.t1e = t1e
		self.set_Jg(J, g)
		if set_dchi:
			ax, rh = self.lanth_axrh[lanthanide]
			self.axrh = np.array([ax,rh])*1E-32

	def set_Jg(self, J, g):
		"""
		Set the magnetic susceptibility absolute magnitude from J/g.

		This is achieved using the following formula:

		.. math::
			\\mu_{eff}=g\\mu_B\\sqrt{J(J+1)}

		Parameters
		----------
		J : str
			the total spin angular momentum quantum number
		g : bool, optional
			the Lande g-factor
		"""
		self.mueff = g * self.MUB * (J*(J+1))**0.5

	def info(self, comment=True):
		"""
		Get basic information about the Metal object

		This is returned as a string in human readable units
		This is also the file format for saving the tensor
		
		Parameters
		----------
		comment : bool (optional)
			if True, each line has a '#' placed at the front

		Returns
		-------
		information : str
			a string containing basic information about the Metal

		Examples
		--------
		>>> metal = Metal()
		>>> metal.set_lanthanide('Er')
		>>> metal.info()
		# ax    | 1E-32 m^3 :   -11.600
		# rh    | 1E-32 m^3 :    -8.600
		# x     |   1E-10 m :     0.000
		# y     |   1E-10 m :     0.000
		# z     |   1E-10 m :     0.000
		# a     |       deg :     0.000
		# b     |       deg :     0.000
		# g     |       deg :     0.000
		# mueff |        Bm :     9.581
		# shift |       ppm :     0.000
		# B0    |         T :    18.790
		# temp  |         K :   298.150
		# t1e   |        ps :     0.189
		# taur  |        ns :     0.000
		"""
		l = "{0:<6}| {1:>9} : {2:9.3f}\n"
		if comment:
			l = '# ' + l 
		i  = l.format('ax','1E-32 m^3',self.axrh[0]*1E32)
		i += l.format('rh','1E-32 m^3',self.axrh[1]*1E32)
		i += l.format('x','1E-10 m',self.position[0]*1E10)
		i += l.format('y','1E-10 m',self.position[1]*1E10)
		i += l.format('z','1E-10 m',self.position[2]*1E10)
		i += l.format('a','deg',self.eulers[0]*(180.0/np.pi))
		i += l.format('b','deg',self.eulers[1]*(180.0/np.pi))
		i += l.format('g','deg',self.eulers[2]*(180.0/np.pi))
		i += l.format('mueff','Bm',self.mueff/self.MUB)
		i += l.format('shift','ppm',self.shift)
		i += l.format('B0','T',self.B0)
		i += l.format('temp','K',self.temperature)
		if self.t1e:
			i += l.format('t1e','ps',self.t1e*1E12)
		else:
			i += l.format('t1e','ps',0.0)
		if self.taur:
			i += l.format('taur','ns',self.taur*1E9)
		else:
			i += l.format('taur','ns',0.0)
		return i

	def get_params(self, params):
		"""
		Get tensor parameters that have been scaled appropriately

		This is often used to get parameter values during fitting where
		floating point errors would otherwise occur on the small values
		encountered.

		Parameters
		----------
		params : list of str
			each element of the list is a string that corresponds to 
			an attribute of the Metal to be retrieved.

		Returns
		-------
		scaled_params : list
			a list with respective scaled parameter values from the input.

		Examples
		--------
		>>> metal = Metal(axrh=[20E-32, 3E-32],position=[0.0,10E-10,-5E-10])
		>>> metal.get_params(['ax','rh','x','y','z'])
		[20.0, 3.0, 0.0, 10.0, -5.0]
		"""
		pars = []
		for param in params:
			scale = self.fit_scaling.get(param, 1.0)
			pars.append(scale * getattr(self, param))
		return pars

	def set_params(self, paramValues):
		"""
		Set tensor parameters that have been scaled appropriately

		This is the inverse of the method <get_params>

		Parameters
		----------
		paramValues : list of tuple
			each element is a tuple (variable, value) where 'variable'
			is the string indentifying the attribute to be set, and 'value'
			is the corresponding value

		Examples
		--------
		>>> metal = Metal()
		>>> metal.set_params([('ax',20.0),('rh',3.0)])
		>>> metal.axrh
		[2.e-31 3.e-32]
		"""
		for param, value in paramValues:
			scale = self.fit_scaling.get(param, 1.0)
			setattr(self, param, value/scale)

	@property
	def x(self):
		"""x coordinate"""
		return self.position[0]
	@x.setter
	def x(self, value):
		self.position[0] = value
	@property
	def y(self):
		"""y coordinate"""
		return self.position[1]
	@y.setter
	def y(self, value):
		self.position[1] = value
	@property
	def z(self):
		"""z coordinate"""
		return self.position[2]
	@z.setter
	def z(self, value):
		self.position[2] = value
	@property
	def a(self):
		"""alpha euler anglue"""
		return self.eulers[0]
	@a.setter
	def a(self, value):
		self.eulers[0] = value
	@property
	def b(self):
		"""beta euler anglue"""
		return self.eulers[1]
	@b.setter
	def b(self, value):
		self.eulers[1] = value
	@property
	def g(self):
		"""gamma euler anglue"""
		return self.eulers[2]
	@g.setter
	def g(self, value):
		self.eulers[2] = value
	@property
	def ax(self):
		"""axiality"""
		return self.axrh[0]
	@ax.setter
	def ax(self, value):
		self.axrh[0] = value
	@property
	def rh(self):
		"""rhombicity"""
		return self.axrh[1]
	@rh.setter
	def rh(self, value):
		self.axrh[1] = value
	@property
	def iso(self):
		"""isotropy"""
		return self.isotropy
	@iso.setter
	def iso(self, value):
		self.isotropy = value
	@property
	def B0_MHz(self):
		"""1H NMR frequency for the given field in MHz"""
		return self.B0 * 42.57747892
	@B0_MHz.setter
	def B0_MHz(self, value):
		self.B0 = value / 42.57747892
	@property
	def gax(self):
		"""axial componenet of spectral power density tensor"""
		return self.g_axrh[0]
	@gax.setter
	def gax(self, value):
		self.g_axrh[0] = value
	@property
	def grh(self):
		"""axial componenet of spectral power density tensor"""
		return self.g_axrh[1]
	@gax.setter
	def grh(self, value):
		self.g_axrh[1] = value

	@property
	def eigenvalues(self):
		"""The eigenvalues defining the magnitude of the principle axes"""
		return anisotropy_to_eigenvalues(self.axrh) + self.isotropy

	@eigenvalues.setter
	def eigenvalues(self, newEigenvalues):
		self.axrh = eigenvalues_to_anisotropy(newEigenvalues)
		# self.isotropy = np.round(np.sum(newEigenvalues)/3., 40)
		self.isotropy = np.sum(newEigenvalues)/3.

	@property
	def isotropy(self):
		"""The magnidue of the isotropic component of the tensor"""
		return (self.MU0 * self.mueff**2) / (3*self.K*self.temperature)

	@isotropy.setter
	def isotropy(self, newIsotropy):
		if newIsotropy<=0:
			newIsotropy = 0.0
		self.mueff = ((newIsotropy*3*self.K*self.temperature) / self.MU0)**0.5

	@property
	def g_eigenvalues(self):
		"""The eigenvalues defining the magnitude of the principle axes"""
		return anisotropy_to_eigenvalues(self.g_axrh) + self.g_isotropy

	@g_eigenvalues.setter
	def g_eigenvalues(self, newEigenvalues):
		self.g_axrh = eigenvalues_to_anisotropy(newEigenvalues)
		# self.isotropy = np.round(np.sum(newEigenvalues)/3., 40)
		self.g_isotropy = np.sum(newEigenvalues)/3.

	@property
	def g_isotropy(self):
		"""Estimate of the spectral power density tensor isotropy"""
		return (self.t1e * self.isotropy * self.K * self.temperature) / self.MU0

	@g_isotropy.setter
	def g_isotropy(self, newIsotropy):
		if newIsotropy<=0 or self.isotropy<=0:
			self.t1e = 0.0
		else:	
			self.t1e = (newIsotropy * self.MU0) / (
				self.isotropy * self.K * self.temperature)

	@property
	def rotationMatrix(self):
		"""The rotation matrix as defined by the euler angles"""
		return euler_to_matrix(self.eulers)

	@rotationMatrix.setter
	def rotationMatrix(self, newMatrix):
		self.eulers = unique_eulers(matrix_to_euler(newMatrix))

	@property
	def tensor(self):
		"""The magnetic susceptibility tensor matrix representation"""
		R = self.rotationMatrix
		return R.dot(np.diag(self.eigenvalues)).dot(R.T)

	@tensor.setter
	def tensor(self, newTensor):
		eigenvals, eigenvecs = np.linalg.eigh(newTensor)
		eigs = zip(eigenvals, np.array(eigenvecs).T)
		iso = np.sum(eigenvals)/3.
		eigenvals, (x, y, z) = zip(*sorted(eigs, key=lambda v: abs(v[0]-iso)))
		eigenvecs = x * z.dot(np.cross(x,y)), y, z
		rotationMatrix = np.vstack(eigenvecs).T
		eulers = unique_eulers(matrix_to_euler(rotationMatrix))
		self.eulers = np.array(eulers, dtype=float)
		self.eigenvalues = eigenvals

	@property
	def tensor_traceless(self):
		"""The traceless magnetic susceptibility tensor matrix representation"""
		return self.tensor - np.identity(3)*self.isotropy

	@property
	def alignment_factor(self):
		"""Factor for conversion between magnetic susceptibility 
		and alignment tensors"""
		return (self.B0**2) / (15 * self.MU0 * self.K * self.temperature)

	@property
	def saupe_factor(self):
		"""Factor for conversion between magnetic susceptibility 
		and saupe tensors"""
		return (3./2.)*self.alignment_factor

	@property
	def tensor_alignment(self):
		"""The alignment tensor matrix representation"""
		return self.alignment_factor * self.tensor_traceless

	@tensor_alignment.setter
	def tensor_alignment(self, new_tensor_alignment):
		old_mueff = self.mueff
		self.tensor = new_tensor_alignment / self.alignment_factor
		self.mueff = old_mueff

	@property
	def tensor_saupe(self):
		"""The saupe tensor matrix representation"""
		return self.saupe_factor * self.tensor_traceless

	@tensor_saupe.setter
	def tensor_saupe(self, new_tensor_saupe):
		old_mueff = self.mueff
		self.tensor = new_tensor_saupe / self.saupe_factor
		self.mueff = old_mueff

	@property
	def upper_triang(self):
		"""Fetch 5 unique matrix element defining the magnetic 
		susceptibility tensor"""
		return self.tensor_traceless[self.upper_coords]

	@upper_triang.setter
	def upper_triang(self, elements):
		old_mueff = self.mueff
		newTensor = np.zeros(9).reshape(3,3)
		newTensor[self.upper_coords] = elements
		newTensor[self.lower_coords] = elements
		newTensor[2,2] = - elements[0] - elements[1]
		self.tensor = newTensor
		self.mueff = old_mueff

	@property
	def upper_triang_alignment(self):
		"""Fetch 5 unique matrix element defining the alignment tensor"""
		return self.tensor_alignment[self.upper_coords]

	@upper_triang_alignment.setter
	def upper_triang_alignment(self, elements):
		newTensor = np.zeros(9).reshape(3,3)
		newTensor[self.upper_coords] = elements
		newTensor[self.lower_coords] = elements
		newTensor[2,2] = - elements[0] - elements[1]
		self.tensor_alignment = newTensor

	@property
	def upper_triang_saupe(self):
		"""Fetch 5 unique matrix element defining the saupe tensor"""
		return self.tensor_saupe[self.upper_coords]

	@upper_triang.setter
	def upper_triang_saupe(self, elements):
		newTensor = np.zeros(9).reshape(3,3)
		newTensor[self.upper_coords] = elements
		newTensor[self.lower_coords] = elements
		newTensor[2,2] = - elements[0] - elements[1]
		self.tensor_saupe = newTensor

	@property
	def g_tensor(self):
		"""The magnetic susceptibility tensor matrix representation"""
		R = self.rotationMatrix
		return R.dot(np.diag(self.g_eigenvalues)).dot(R.T)

	@g_tensor.setter
	def g_tensor(self, newTensor):
		eigenvals, eigenvecs = np.linalg.eigh(newTensor)
		eigs = zip(eigenvals, np.array(eigenvecs).T)
		iso = np.sum(eigenvals)/3.
		eigenvals, (x, y, z) = zip(*sorted(eigs, key=lambda x: abs(x[0]-iso)))
		eigenvecs = x * z.dot(np.cross(x,y)), y, z
		rotationMatrix = np.vstack(eigenvecs).T
		eulers = unique_eulers(matrix_to_euler(rotationMatrix))
		self.eulers = np.array(eulers, dtype=float)
		self.g_eigenvalues = eigenvals

	def set_utr(self):
		"""
		Modify current tensor parameters to unique tensor representation (UTR)

		Note that multiple axial/rhombic and euler angles can give congruent
		tensors. 
		This method ensures that identical tensors may always be compared
		by using Numbat style representation.
		"""
		self.tensor = self.tensor

	def atom_set_position(self, atom):
		"""
		Set the position of the Metal object to that of an atom

		Parameters
		----------
		atom : biopython atom object
			must have 'position' attribute
		"""
		self.position = atom.position

	def average(self, metals):
		"""
		Set the attributes of the current instance to the average
		of a list of provided tensor objects

		WARNING: averging is unstable for spectral power density <g_tensor>

		Parameters
		----------
		metals : a list of Metal objects
			the average of attributes of this list will be taken
		"""
		self.g_tensor = sum([m.g_tensor for m in metals])/len(metals)
		self.tensor = sum([m.tensor for m in metals])/len(metals)
		self.position = sum(m.position for m in metals)/len(metals)
		self.taur = sum([m.taur for m in metals])/len(metals)

	################################
	# Methods for PCS calculations #
	################################

	def dipole_shift_tensor(self, position):
		"""
		Calculate the chemical shift tensor at the given postition

		This arises due to the paramagnetic dipole tensor field

		Parameters
		----------
		position : array floats
			the position (x, y, z) in meters

		Returns
		-------
		dipole_shift_tensor : 3x3 array
			the tensor describing chemical shift at the nuclear position
		"""
		pos = np.array(position, dtype=float) - self.position
		distance = np.linalg.norm(pos)
		preFactor = 1./(4.*np.pi)
		p1 = (1./distance**5)*np.kron(pos,pos).reshape(3,3)
		p2 = (1./distance**3)*np.identity(3)
		return (preFactor * (3.*p1 - p2)).dot(self.tensor)

	def fast_dipole_shift_tensor(self, posarray):
		"""
		A vectorised version of 
		:meth:`paramagpy.metal.Metal.dipole_shift_tensor`

		This is generally used for fast calculations.

		Parameters
		----------
		posarray : array
			an array of positions with shape (n,3)

		Returns
		-------
		dipole_shift_tensor_array : array
			and array of dipole shift tensors at corresponding positions.
			This has shape (n,3,3)
		"""
		pos = posarray - self.position
		distance = np.linalg.norm(pos, axis=1)[:,None,None]
		preFactor = 1./(4.*np.pi)
		p1 = (1./distance**5)*np.einsum('ij,ik->ijk', pos, pos)
		p2 = (1./distance**3)*np.identity(3)
		ds = preFactor*(3.*p1 - p2).dot(self.tensor)
		return ds

	def pcs(self, position):
		"""
		Calculate the psuedo-contact shift at the given postition

		Parameters
		----------
		position : array floats
			the position (x, y, z) in meters

		Returns
		-------
		pcs : float
			the pseudo-contact shift in parts-per-million (ppm)

		Examples
		--------
		>>> metal = Metal()
		>>> metal.set_lanthanide('Er')
		>>> metal.pcs([0.,0.,10E-10])
		-6.153991132886608
		"""
		val = self.dipole_shift_tensor(position).trace()/3.
		return 1E6*val + self.shift

	def atom_pcs(self, atom, racs=False, rads=False):
		"""
		Calculate the psuedo-contact shift at the given atom

		Parameters
		----------
		atom : biopython atom object
			must have 'position' attribute
		racs : bool (optional)
			when True, RACS (residual anisotropic chemical shielding) 
			correction is included. Default is False
		rads : bool (optional)
			when True, RADS (residual anisotropic dipolar shielding) 
			correction is included. Defualt is False

		Returns
		-------
		pcs : float
			the pseudo-contact shift in parts-per-million (ppm)
		"""
		value = self.pcs(atom.position)
		if racs:
			value += self.racs(atom.csa)
		if rads:
			value += self.rads(atom.position)
		return value

	def fast_pcs(self, posarray):
		"""
		A vectorised version of :meth:`paramagpy.metal.Metal.pcs`

		This efficient algorithm calculates the PCSs for an array of 
		positions and is best used where speed is required for fitting.

		Parameters
		----------
		posarray : array with shape (n,3)
			array of 'n' positions (x, y, z) in meters

		Returns
		-------
		pcs : array of floats with shape (n,1)
			the peudo-contact shift in parts-per-million (ppm)
		"""
		pos = posarray - self.position
		dist = np.linalg.norm(pos, axis=1)
		dot1 = np.einsum('ij,jk->ik', pos, self.tensor_traceless)
		dot2 = np.einsum('ij,ij->i', pos, dot1)
		return 1E6*(1./(4.*np.pi))*(dot2/dist**5) + self.shift

	def rads(self, position):
		"""
		Calculate the residual anisotropic dipolar shift at the 
		given postition. 

		The partial alignment induced by an anisotropic 
		magnetic susecptiblity causes the dipole shift tensor at a nuclear
		position to average to a value different to the PCS.

		Parameters
		----------
		position : array floats
			the position (x, y, z) in meters

		Returns
		-------
		rads : float
			the residual anisotropic dipole shift in parts-per-million (ppm)
		"""
		ds = self.dipole_shift_tensor(position)
		rads = self.tensor_saupe.dot(ds).trace()/3.
		return 1E6*rads

	def fast_rads(self, posarray):
		"""
		A vectorised version of :meth:`paramagpy.metal.Metal.rads`

		This is generally used when speed is required for fitting

		Parameters
		----------
		posarray : array with shape (n,3)
			an array of 'n' positions (x, y, z) in meters

		Returns
		-------
		rads_array : array of floats with shape (n,1)
			the residual anisotropic dipole shift in parts-per-million (ppm)
		"""
		ds = self.fast_dipole_shift_tensor(posarray)
		rads = np.einsum('jk,ikl->ijl',self.tensor_saupe,ds)
		return 1E6*rads.trace(axis1=1,axis2=2)/3.

	def racs(self, csa):
		"""
		Calculate the residual anisotropic chemical shift at the 
		given postition. 

		The partial alignment induced by an anisotropic 
		magnetic susecptiblity causes the chemical shift tensor at a nuclear
		position to average to a value different to the isotropic value.

		Parameters
		----------
		csa : 3 x 3 array
			the chemical shift anisotropy tensor

		Returns
		-------
		racs : float
			the residual anisotropic chemical shift in parts-per-million (ppm)
		"""
		racs = self.tensor_alignment.dot(csa).trace()
		return 1E6*racs

	def fast_racs(self, csaarray):
		"""
		A vectorised version of :meth:`paramagpy.metal.Metal.racs`

		This is generally used when speed is required for fitting

		Parameters
		----------
		csaarray : array with shape (n,3,3)
			array of chemical shift anisotropy tensors

		Returns
		-------
		racs_array : array of floats with shape (n,1)
			the residual anisotropic chemical shift in parts-per-million (ppm)
		"""
		racs = 2*np.einsum('jk,ikl->ijl',self.tensor_alignment,csaarray)
		return 1E6*racs.trace(axis1=1,axis2=2)/3.

	def pcs_gradient(self, position):
		"""
		Calculate the gradient of the psuedo-contact shift 
		at the given postition.
		This equation uses analytic partial derivatives calculated in 
		Mathematica.

		Parameters
		----------
		position : array floats
			the position (x, y, z) in metres

		Returns
		-------
		gradient : array of floats
			the [x,y,z] gradient vector for the pseudo-contact shift 
			in parts-per-million (ppm) per metre
		"""
		pos = position - self.position
		r = np.linalg.norm(pos)
		x, y, z = pos
		xx, yy, xy, xz, yz = self.upper_triang
		dnt = (x**2-z**2)*xx + (y**2-z**2)*yy + 2*x*y*xy + 2*x*z*xz + 2*y*z*yz
		dxt = 2*x*xx + 2*y*xy + 2*z*xz
		dx = dxt / (4*np.pi*r**5) - 5*x*dnt / (4*np.pi*r**7)
		dyt = 2*x*xy + 2*y*yy + 2*z*yz
		dy = dyt / (4*np.pi*r**5) - 5*y*dnt / (4*np.pi*r**7)
		dzt = -2*z*xx - 2*z*yy + 2*x*xz +2*y*yz
		dz = dzt / (4*np.pi*r**5) - 5*z*dnt / (4*np.pi*r**7)
		return 1E6 * np.array([dx,dy,dz])

	################################
	# Methods for PRE calculations #
	################################

	@staticmethod
	def spec_dens(tau, omega):
		"""
		A spectral density function with Lorentzian shape:

		.. math::
			\\mathbb{J}(\\tau,\\omega)=\\frac{\\tau}{1+(\\omega\\tau)^2}

		Parameters
		----------
		tau : float
			correaltion time
		omega : float
			frequency

		Returns
		-------
		value : float
			the value of the spectral denstiy
		"""
		return tau/(1+(omega*tau)**2)

	@staticmethod
	def first_invariant_squared(t):
		"""
		Calculate the antisymmetric contribution to relaxation via the
		first invariant of a tensor.

		This is required for PRE calculations using the shilding tensor

		Parameters
		----------
		tensor : 3x3 matrix
			a second rank tensor

		Returns
		-------
		firstInvariantSquared : float
			the first invariant squared of the shift tensor
		"""
		xy, yx = t[0,1], t[1,0]
		xz, zx = t[0,2], t[2,0]
		yz, zy = t[1,2], t[2,1]
		sis = (xy-yx)**2 + (xz-zx)**2 + (yz-zy)**2
		return sis

	@staticmethod
	def second_invariant_squared(t):
		"""
		Calculate the second invariant squared of a tensor.

		This is required for PRE calculations using the shilding tensor

		Parameters
		----------
		tensor : 3x3 matrix
			a second rank tensor

		Returns
		-------
		secondInvariantSquared : float
			the second invariant squared of the shift tensor
		"""
		xx, yy, zz = t[0,0], t[1,1], t[2,2]
		xy, yx = t[0,1], t[1,0]
		xz, zx = t[0,2], t[2,0]
		yz, zy = t[1,2], t[2,1]
		sis = xx**2 + yy**2 + zz**2 -xx*yy - xx*zz - yy*zz
		sis += 0.75*((xy+yx)**2 + (xz+zx)**2 + (yz+zy)**2)
		return sis

	@staticmethod
	def fast_first_invariant_squared(t):
		"""
		Vectorised version of 
		:meth:`paramagpy.metal.Metal.first_invariant_squared`

		This is generally used for speed in fitting PRE data

		Parameters
		----------
		tensorarray : array with shape (n,3,3)
			array of shielding tensors

		Returns
		-------
		firstInvariantSquared : array with shape (n,1)
			the first invariants squared of the tensors
		"""
		xy, yx = t[:,0,1], t[:,1,0]
		xz, zx = t[:,0,2], t[:,2,0]
		yz, zy = t[:,1,2], t[:,2,1]
		sis = (xy-yx)**2 + (xz-zx)**2 + (yz-zy)**2
		return sis

	@staticmethod
	def fast_second_invariant_squared(t):
		"""
		Vectorised version of 
		:meth:`paramagpy.metal.Metal.second_invariant_squared`

		This is generally used for speed in fitting PRE data

		Parameters
		----------
		tensorarray : array with shape (n,3,3)
			array of shielding tensors

		Returns
		-------
		secondInvariantSquared : array with shape (n,1)
			the second invariants squared of the tensors
		"""
		xx, yy, zz = t[:,0,0], t[:,1,1], t[:,2,2]
		xy, yx = t[:,0,1], t[:,1,0]
		xz, zx = t[:,0,2], t[:,2,0]
		yz, zy = t[:,1,2], t[:,2,1]
		sis = xx**2 + yy**2 + zz**2 -xx*yy - xx*zz - yy*zz
		sis += 0.75*((xy+yx)**2 + (xz+zx)**2 + (yz+zy)**2)
		return sis

	def dsa_r1(self, position, gamma, csa=0.0):
		"""
		Calculate R1 relaxation due to Curie Spin

		If the metal has an anisotropic magnetic susceptibility, this is
		taken into account, resulting in orientation dependent PRE as
		predicted by Vega and Fiat. CSA cross-correlated relaxation may
		be included by providing an appropriate CSA tensor.

		Parameters
		----------
		position : array of floats
			three coordinates (x,y,z) in meters
		gamma : float
			the gyromagnetic ratio of the spin
		csa : 3x3 matrix (optional)
			the CSA tensor of the given spin.
			This defualts to 0.0, meaning CSAxDSA crosscorrelation is
			not accounted for.

		Returns
		-------
		value : float
			The R1 relaxation rate in /s
		"""
		ds = self.dipole_shift_tensor(position)
		fis = self.first_invariant_squared(ds + csa)
		sis = self.second_invariant_squared(ds + csa)
		if isinstance(csa, np.ndarray):
			fis -= self.first_invariant_squared(csa)
			sis -= self.second_invariant_squared(csa)
		omega = self.B0 * gamma
		pfasym = (1./2. ) * fis * omega**2
		pfsymm = (2./15.) * sis * omega**2
		rate  = pfasym * self.spec_dens(self.taur, 3*omega)
		rate += pfsymm * self.spec_dens(self.taur,   omega)
		return rate

	def fast_dsa_r1(self, posarray, gammaarray, csaarray=0.0):
		"""
		Vectorised version of :meth:`paramagpy.metal.Metal.dsa_r1`

		This is generally used for speed in fitting PRE data

		Parameters
		----------
		posarray : array with shape (n,3)
			array of positions in meters
		gammaarray : array with shape (n,3)
			array of gyromagnetic ratios of the spins
		csaarray : array with shape (m,3,3) (optional)
			array of CSA tensors of the spins.
			This defualts to 0.0, meaning CSAxDSA crosscorrelation is
			not accounted for.

		Returns
		-------
		rates : array with shape (n,1)
			The R1 relaxation rates in /s
		"""
		ds = self.fast_dipole_shift_tensor(posarray)
		fis = self.fast_first_invariant_squared(ds + csaarray)
		sis = self.fast_second_invariant_squared(ds + csaarray)
		if isinstance(csaarray, np.ndarray):
			fis -= self.fast_first_invariant_squared(csaarray)
			sis -= self.fast_second_invariant_squared(csaarray)
		omegas = self.B0 * gammaarray
		pfasym = (1./2. ) * fis * omegas**2
		pfsymm = (2./15.) * sis * omegas**2
		rate  = pfasym * self.spec_dens(self.taur, 3*omegas)
		rate += pfsymm * self.spec_dens(self.taur,   omegas)
		return rate

	def dsa_r2(self, position, gamma, csa=0.0):
		"""
		Calculate R2 relaxation due to Curie Spin

		If the metal has an anisotropic magnetic susceptibility, this is
		taken into account, resulting in orientation dependent PRE as
		predicted by Vega and Fiat. CSA cross-correlated relaxation may
		be included by providing an appropriate CSA tensor.

		Parameters
		----------
		position : array of floats
			three coordinates (x,y,z)
		gamma : float
			the gyromagnetic ratio of the spin
		csa : 3x3 matrix (optional)
			the CSA tensor of the given spin.
			This defualts to 0.0, meaning CSAxDSA crosscorrelation is
			not accounted for.

		Returns
		-------
		value : float
			The R2 relaxation rate in /s
		"""
		ds = self.dipole_shift_tensor(position)
		fis = self.first_invariant_squared(ds + csa)
		sis = self.second_invariant_squared(ds + csa)
		if isinstance(csa, np.ndarray):
			fis -= self.first_invariant_squared(csa)
			sis -= self.second_invariant_squared(csa)
		omega = self.B0 * gamma
		pfasym = (1./4. ) * fis * omega**2
		pfsymm = (1./45.) * sis * omega**2
		rate  = pfasym * self.spec_dens(self.taur, 3*omega)
		rate += pfsymm *(4*self.spec_dens(self.taur, 0.   ) +
						 3*self.spec_dens(self.taur, omega))
		return rate

	def fast_dsa_r2(self, posarray, gammaarray, csaarray=0.0):
		"""
		Vectorised version of :meth:`paramagpy.metal.Metal.dsa_r2`

		This is generally used for speed in fitting PRE data.

		Parameters
		----------
		posarray : array with shape (n,3)
			array of positions in meters
		gammaarray : array with shape (n,3)
			array of gyromagnetic ratios of the spins
		csaarray : array with shape (m,3,3) (optional)
			array of CSA tensors of the spins.
			This defualts to 0.0, meaning CSAxDSA crosscorrelation is
			not accounted for.

		Returns
		-------
		rates : array with shape (n,1)
			The R2 relaxation rates in /s
		"""
		ds = self.fast_dipole_shift_tensor(posarray)
		fis = self.fast_first_invariant_squared(ds + csaarray)
		sis = self.fast_second_invariant_squared(ds + csaarray)
		if isinstance(csaarray, np.ndarray):
			fis -= self.fast_first_invariant_squared(csaarray)
			sis -= self.fast_second_invariant_squared(csaarray)
		omegas = self.B0 * gammaarray
		pfasym = (1./4. ) * fis * omegas**2
		pfsymm = (1./45.) * sis * omegas**2
		rate  = pfasym * self.spec_dens(self.taur, 3*omegas)
		rate += pfsymm *(4*self.spec_dens(self.taur, 0.   ) +
						 3*self.spec_dens(self.taur, omegas))
		return rate

	def sbm_r1(self, position, gamma):
		"""
		Calculate R1 relaxation due to Solomon-Bloembergen-Morgan theory

		Parameters
		----------
		position : array of floats
			three coordinates (x,y,z)
		gamma : float
			the gyromagnetic ratio of the spin

		Returns
		-------
		value : float
			The R1 relaxation rate in /s
		"""
		distance = np.linalg.norm(position-self.position)
		p1 = self.MU0/(4.*np.pi)
		p2 = (gamma * self.mueff)/distance**3
		rate = (2./15.)*(p1*p2)**2 * (
			3*self.spec_dens(self.tauc, self.B0*gamma) + 
			7*self.spec_dens(self.tauc, self.B0*self.GAMMA))
		return rate

	def fast_sbm_r1(self, posarray, gammaarray):
		"""
		Vectorised version of :meth:`paramagpy.metal.Metal.sbm_r1`

		This is generally used for speed in fitting PRE data

		Parameters
		----------
		posarray : array with shape (n,3)
			array of positions in meters
		gammaarray : array with shape (n,3)
			array of gyromagnetic ratios of the spins

		Returns
		-------
		rates : array with shape (n,1)
			The R1 relaxation rates in /s
		"""
		distance = np.linalg.norm(posarray-self.position, axis=1)
		p1 = self.MU0/(4.*np.pi)
		p2 = (gammaarray * self.mueff)/distance**3
		rate = (2./15.)*(p1*p2)**2 * (
			3*self.spec_dens(self.tauc, self.B0*gammaarray) + 
			7*self.spec_dens(self.tauc, self.B0*self.GAMMA))
		return rate

	def sbm_r2(self, position, gamma):
		"""
		Calculate R2 relaxation due to Solomon-Bloembergen-Morgan theory

		Parameters
		----------
		position : array of floats
			three coordinates (x,y,z)
		gamma : float
			the gyromagnetic ratio of the spin

		Returns
		-------
		value : float
			The R2 relaxation rate in /s
		"""
		distance = np.linalg.norm(position-self.position)
		p1 = self.MU0/(4.*np.pi)
		p2 = (gamma * self.mueff)/distance**3
		rate = (1./15.)*(p1*p2)**2 * (
			 4*self.spec_dens(self.tauc, 0.) +
			 3*self.spec_dens(self.tauc, self.B0*gamma) + 
			13*self.spec_dens(self.tauc, self.B0*self.GAMMA))
		return rate

	def fast_sbm_r2(self, posarray, gammaarray):
		"""
		Vectorised version of :meth:`paramagpy.metal.Metal.sbm_r2`

		This is generally used for speed in fitting PRE data

		Parameters
		----------
		posarray : array with shape (n,3)
			array of positions in meters
		gammaarray : array with shape (n,3)
			array of gyromagnetic ratios of the spins

		Returns
		-------
		rates : array with shape (n,1)
			The R2 relaxation rates in /s
		"""
		distance = np.linalg.norm(posarray-self.position, axis=1)
		p1 = self.MU0/(4.*np.pi)
		p2 = (gammaarray * self.mueff)/distance**3
		rate = (1./15.)*(p1*p2)**2 * (
			 4*self.spec_dens(self.tauc, 0.) +
			 3*self.spec_dens(self.tauc, self.B0*gammaarray) + 
			13*self.spec_dens(self.tauc, self.B0*self.GAMMA))
		return rate

	def g_sbm_r1(self, position, gamma):
		"""
		Calculate R1 relaxation due to Solomon-Bloembergen-Morgan theory
		from anisotropic power spectral density tensor

		Parameters
		----------
		position : array of floats
			three coordinates (x,y,z)
		gamma : float
			the gyromagnetic ratio of the spin

		Returns
		-------
		value : float
			The R1 relaxation rate in /s
		"""
		pos = np.array(position, dtype=float) - self.position
		distance = np.linalg.norm(pos)
		pos_unit = pos / distance
		preFactor = (2./3.) * (self.MU0 / (4*np.pi))**2
		preFactor *= gamma**2 / distance**6
		D = 3*np.kron(pos_unit,pos_unit).reshape(3,3) - np.identity(3)
		return preFactor * (D.dot(D)).dot(self.g_tensor).trace()

	def fast_g_sbm_r1(self, posarray, gammaarray):
		"""
		Vectorised version of :meth:`paramagpy.metal.Metal.g_sbm_r1`

		This is generally used for speed in fitting PRE data

		Parameters
		----------
		posarray : array with shape (n,3)
			array of positions in meters
		gammaarray : array with shape (n,3)
			array of gyromagnetic ratios of the spins

		Returns
		-------
		rates : array with shape (n,1)
			The R1 relaxation rates in /s
		"""
		n = len(gammaarray)
		pos = posarray - self.position
		distance = np.linalg.norm(pos, axis=1)
		pos_unit = (pos.T / distance).T
		preFactor = (2./3.) * (self.MU0 / (4*np.pi))**2
		preFactor *= gammaarray**2 / distance**6
		D1 = np.einsum('ij,ik->ijk', pos_unit, pos_unit)
		D2 = np.tile(np.identity(3), n).T.reshape(n,3,3)
		D = 3*D1 - D2
		D_squared = np.einsum('ijk,ilk->ijl', D,D)
		tmp = np.einsum('ijk,lk->ijl', D_squared, self.g_tensor)
		return preFactor * np.einsum('ijj->i', tmp)


	def pre(self, position, gamma, rtype, dsa=True, sbm=True, 
			gsbm=False, csa=0.0):
		"""
		Calculate the PRE for a set of spins using Curie and or SBM theory

		Parameters
		----------
		position : array of floats
			position in meters
		gamma : float
			gyromagnetic ratio of the spin
		rtype : str
			either 'r1' or 'r2', the relaxation type
		dsa : bool (optional)
			when True (defualt), DSA or Curie spin relaxation is included
		sbm : bool (optional)
			when True (defualt), SBM spin relaxation is included
		gsbm : bool (optional)
			when True (default=False), anisotropic dipolar relaxation is 
			included using the spectral power density gensor <g_tensor>
			NOTE: when true, ignores relaxation of type SBM
			NOTE: only implemented for R1 relaxation calculations
		csa : array with shape (3,3) (optional)
			CSA tensor of the spin.
			This defualts to 0.0, meaning CSAxDSA crosscorrelation is
			not accounted for. 

		Returns
		-------
		rate : float
			The PRE rate in /s
		"""
		if gsbm:
			sbm = False
			if rtype=='r2':
				raise NotImplementedError(
					"Anisotropic dipolar relaxation has not been implement for R2 caluculations yet.")
		rate = 0.0
		if rtype=='r1':
			if dsa:
				rate += self.dsa_r1(position, gamma, csa)
			if sbm:
				rate += self.sbm_r1(position, gamma)
			if gsbm:
				rate += self.g_sbm_r1(position, gamma)
		elif rtype=='r2':
			if dsa:
				rate += self.dsa_r2(position, gamma, csa)
			if sbm:
				rate += self.sbm_r2(position, gamma)
		return rate

	def fast_pre(self, posarray, gammaarray, rtype, 
		dsa=True, sbm=True, gsbm=False, csaarray=0.0):
		"""
		Calculate the PRE for a set of spins using Curie and or SBM theory

		Parameters
		----------
		posarray : array with shape (n,3)
			array of positions in meters
		gammaarray : array with shape (n,3)
			array of gyromagnetic ratios of the spins
		rtype : str
			either 'r1' or 'r2', the relaxation type
		dsa : bool (optional)
			when True (defualt), DSA or Curie spin relaxation is included
		sbm : bool (optional)
			when True (defualt), SBM spin relaxation is included
		gsbm : bool (optional)
			when True (default=False), anisotropic dipolar relaxation is 
			included using the spectral power density gensor <g_tensor>
			NOTE: when true, ignores relaxation of type SBM
			NOTE: only implemented for R1 relaxation calculations
		csaarray : array with shape (m,3,3) (optional)
			array of CSA tensors of the spins.
			This defualts to 0.0, meaning CSAxDSA crosscorrelation is
			not accounted for. 

		Returns
		-------
		rates : array with shape (n,1)
			The PRE rates in /s
		"""
		if gsbm:
			sbm = False
		rates = 0.0
		if rtype=='r1':
			if dsa:
				rates += self.fast_dsa_r1(posarray, gammaarray, csaarray)
			if sbm:
				rates += self.fast_sbm_r1(posarray, gammaarray)
			if gsbm:
				rates += self.fast_g_sbm_r1(posarray, gammaarray)
		elif rtype=='r2':
			if dsa:
				rates += self.fast_dsa_r2(posarray, gammaarray, csaarray)
			if sbm:
				rates += self.fast_sbm_r2(posarray, gammaarray)
		return rates

	def atom_pre(self, atom, rtype='r2', dsa=True, sbm=True, csa=0.0):
		"""
		Calculate the PRE for an atom

		Parameters
		----------
		atom : paramagpy.protein.CustomAtom
			the active nuclear spin for which relaxation will be calculated
			must have attributes 'position' and 'gamma'
		rtype : str
			either 'r1' or 'r2', the relaxation type
		dsa : bool (optional)
			when True (defualt), DSA or Curie spin relaxation is included
		sbm : bool (optional)
			when True (defualt), SBM spin relaxation is included 
		csa : array with shape (3,3) (optional)
			CSA tensor of the spin.
			This defualts to 0.0, meaning CSAxDSA crosscorrelation is
			not accounted for. 

		Returns
		-------
		rate : float
			The PRE rate in /s
		"""
		return self.pre(atom.position, atom.gamma, 
			rtype, dsa=dsa, sbm=sbm, csa=csa)


	################################
	# Methods for CCR calculations #
	################################

	def ccr(self, position, gamma, dipole_shift_tensor):
		"""
		Calculate R2 cross-corelated relaxation due to DDxDSA 

		If the metal has an anisotropic magnetic susceptibility, this is
		taken into account.

		Parameters
		----------
		position : array of floats
			three coordinates (x,y,z)
			this is the position of the nuclear spin
		gamma : float
			the gyromagnetic ratio of the relaxing spin
		dipole_shift_tensor : 3x3 array of floats
			this is the dipole shift tensor arising from 
			the nuclear spin of the coupling partner

		Returns
		-------
		value : float
			The R2 differential line broadening rate in /s
		"""
		# NOTE: B0 factor expresses shielding tensor in ppm
		shield = dipole_shift_tensor/self.B0
		spinUp   = self.dsa_r2(position, gamma, csa= shield)
		spinDown = self.dsa_r2(position, gamma, csa=-shield)
		return spinUp - spinDown

	def fast_ccr(self, posarray, gammaarray, dstarray):
		"""
		Vectorised version of :meth:`paramagpy.metal.Metal.ccr`

		This is generally used for speed in fitting DDxDSA data

		If the metal has an anisotropic magnetic susceptibility, this is
		taken into account.

		Parameters
		----------
		posarray : array with shape (n,3)
			array of positions in meters
		gammaarray : array with shape (n,3)
			array of gyromagnetic ratios of the spins
		dstarray : array with shape (n,3,3)
			array of nuclear dipole shift tensors arising from
			the coupling partners

		Returns
		-------
		rates : array with shape (n,1)
			The R2 differential line broadening rates in /s
		"""
		# NOTE: B0 factor expresses shielding tensor in ppm
		shield = dstarray/self.B0
		spinUp   = self.fast_dsa_r2(posarray, gammaarray, 
			csaarray= shield)
		spinDown = self.fast_dsa_r2(posarray, gammaarray, 
			csaarray=-shield)
		return spinUp - spinDown

	def atom_ccr(self, atom, atomPartner):
		"""
		Calculate R2 cross-corelated relaxation due to DDxDSA 

		Parameters
		----------
		atom : paramagpy.protein.CustomAtom
			the active nuclear spin for which relaxation will be calculated
			must have attributes 'position' and 'gamma'
		atomPartner : paramagpy.protein.CustomAtom
			the coupling parnter nuclear spin
			must have method 'dipole_shift_tensor'

		Returns
		-------
		value : float
			the CCR differential line broadening in Hz
		"""
		dd_tensor = atomPartner.dipole_shift_tensor(atom.position)
		return self.ccr(atom.position, atom.gamma, dd_tensor)
			
	################################
	# Methods for RDC calculations #
	################################

	def rdc(self, vector, gammaProd):
		"""
		Calculate Residual Dipolar Coupling (RDC)

		Parameters
		----------
		vector : array of floats
			internuclear vector (x,y,z) in meters
		gammaProd : float
			the product of gyromagnetic ratios of spin A and B 
			where each has units of rad/s/T

		Returns
		-------
		rdc : float
			the RDC in Hz
		"""
		dist = np.linalg.norm(vector)
		pf = -(self.MU0 * gammaProd * self.HBAR) / (4 * np.pi * dist**5)
		rdc_radians = 3*pf * vector.dot(self.tensor_alignment).dot(vector) 
		return rdc_radians / (2*np.pi)

	def fast_rdc(self, vecarray, gammaProdArray):
		"""
		A vectorised version of :meth:`paramagpy.metal.Metal.rdc` method.

		This is generally used for speed in fitting RDC data

		Parameters
		----------
		vecarray : array with shape (n,3)
			array of internuclear vectors in meters
		gammaProdArray : array with shape (n,1)
			the products of gyromagnetic ratios of spins A and B 
			where each has units of rad/s/T

		Returns
		-------
		rdc_array : array with shape (n,1)
			the RDC values in Hz
		"""
		dist = np.linalg.norm(vecarray, axis=1)
		pf = -(self.MU0 * gammaProdArray * self.HBAR) / (4 * np.pi * dist**5)
		dot1 = np.einsum('ik,jk->ji', self.tensor_alignment, vecarray)
		dot2 = np.einsum('ij,ij->i', vecarray, dot1)
		rdc_radians = 3*pf * dot2
		return rdc_radians / (2*np.pi)

	def atom_rdc(self, atom1, atom2):
		"""
		Calculate the residual dipolar coupling between two atoms

		Parameters
		----------
		atom1 : biopython atom object
			must have 'position' and 'gamma' attribute
		atom1 : biopython atom object
			must have 'position' and 'gamma' attribute

		Returns
		-------
		rdc : float
			the RDC values in Hz
		"""
		vector = atom2.position - atom1.position
		gammaProd = atom1.gamma * atom2.gamma
		return self.rdc(vector, gammaProd)


	####################################
	# Methods for plotting isosurfaces #
	####################################

	def make_mesh(self, density=2, size=40.0, origin=None):
		"""
		Construct a 3D grid of points to map an isosurface

		This is contained in a cube

		Parameters
		----------
		density : int (optional)
			the points per Angstrom in the grid
		size : float (optional)
			the length of one edge of the cube

		Returns
		-------
		mesh : cubic grid array
			This has shape (n,n,n,3) where n is the number of points
			along one edge of the grid. Units are meters
		origin : array of floats, 
			the (x,y,z) location of mesh vertex
		low : array of ints, the integer location of the first 
			point in each dimension
		high : array of ints, the integer location of the last 
			point in each dimension
		points : array of ints, 
			the number of points along each dimension

		"""
		if origin is None:
			origin = self.position

		grid_origin = np.asarray(density * (self.position*1E10 - size/2.0), 
				dtype=int)
		low = grid_origin / float(density)
		high = low + size
		points = np.array([int(density*size)]*3) + 1
		domains = [1E-10*np.linspace(*i) for i in zip(low, high, points)]
		mesh = np.array(np.meshgrid(*domains, indexing='ij')).T
		return mesh, (grid_origin, low, high, points)

	def pcs_mesh(self, mesh):
		"""
		Calculate a PCS value at each location of cubic grid of points

		Parameters
		----------
		mesh : array with shape (n,n,n,3)
			a cubic grid as generated by the method <make_mesh>

		Returns
		-------
		pcs_mesh : array with shape (n,n,n,1)
			The same grid shape, with PCS values at the respective locations
		"""
		og_shape = mesh.shape[:3]
		pcs_mesh = self.fast_pcs(mesh.reshape(np.prod(og_shape),3))
		return pcs_mesh.reshape(*og_shape)

	def pre_mesh(self, mesh, gamma=2*np.pi*42.576E6, rtype='r2', 
		dsa=True, sbm=True):
		"""
		Calculate a PRE value at each location of cubic grid of points

		Parameters
		----------
		mesh : array with shape (n,n,n,3)
			a cubic grid as generated by the method <make_mesh>
		gamma : float
			the gyromagnetic ratio of the spin
		rtype : str
			either 'r1' or 'r2', the relaxation type
		dsa : bool (optional)
			when True (defualt), DSA or Curie spin relaxation is included
		sbm : bool (optional)
			when True (defualt), SBM spin relaxation is included

		Returns
		-------
		pre_mesh : array with shape (n,n,n,1)
			The same grid shape, with PRE values at the respective locations
		"""
		og_shape = mesh.shape[:3]
		flat = mesh.reshape(np.prod(og_shape),3)
		gamarr = np.ones(len(flat))*gamma
		pre_mesh = self.fast_pre(flat, gamarr, rtype=rtype, dsa=dsa, sbm=sbm)
		return pre_mesh.reshape(*og_shape)

	def write_pymol_script(self, isoval=1.0, surfaceName='isomap', 
		scriptName='isomap.pml', meshName='./isomap.pml.ccp4', pdbFile=None):
		"""
		Write a PyMol script to file which allows loading of the 
		isosurface file

		Parameters
		----------
		isoval : float (optional)
			the contour level of the isosurface
		surfaceName : str (optional)
			the name of the isosurface file within PyMol
		scriptName : str (optional)
			the name of the PyMol script to load the tensor isosurface
		meshName : str (optional)
			the name of the binary isosurface file
		pdbFile : str (optional)
			if not <None>, the file name of the PDB file to be loaded with
			the isosurface.
		"""
		posname = "pos_{}".format(surfaceName)
		negname = "neg_{}".format(surfaceName)
		oriname = "ori_{}".format(surfaceName)
		s = "# PyMOL macro for loading tensor isosurface from paramagpy\n"
		s += self.info()+'\n'
		s += "import os, pymol\n"
		s += "curdir = os.path.dirname(pymol.__script__)\n"
		s += "set normalize_ccp4_maps, off\n"
		s += "meshfile = os.path.join(curdir, '{}')\n".format(meshName)
		s += "cmd.load(meshfile, 'isomap', 1, 'ccp4')\n"
		s += "isosurface {}, isomap, {}\n".format(posname,isoval)
		s += "isosurface {}, isomap, {}\n".format(negname,-isoval)
		s += "set transparency, 0.5, {}\n".format(posname)
		s += "set transparency, 0.5, {}\n".format(negname)
		s += "set surface_color, blue, {}\n".format(posname)
		s += "set surface_color, red, {}\n".format(negname)
		s += "pseudoatom {}, pos={}\n".format(oriname,list(self.position*1E10))
		s += "show spheres, {}\n".format(oriname)
		s += "color pink, {}\n".format(oriname)
		if pdbFile:
			protName = ntpath.basename(pdbFile).replace('.pdb','')
			s += "cmd.load(os.path.join(curdir, '{}'),'{}')\n".format(
				pdbFile, protName)
			s += "show_as cartoon, {}\n".format(protName)
		with open(scriptName, 'w') as o:
			o.write(s)
			print("{} script written".format(scriptName))

	def write_isomap(self, mesh, bounds, fileName='isomap.pml.ccp4'):
		"""
		Write a PyMol script to file which allows loading of the 
		isosurface file

		Parameters
		----------
		mesh : 3D scalar np.ndarray of floats
			the scalar field of PCS or PRE values in a cubic grid
		bounds : tuple (origin, low, high, points)
			as generated by :meth:`paramagpy.metal.Metal.make_mesh`
		fileName : str (optional)
			the filename of the isosurface file
		"""
		mesh = np.asarray(mesh, np.float32)
		origin, low, high, points = bounds
		with open(fileName, 'wb') as o:
			for dim in points:
				o.write(struct.pack('i', dim)) # number points per dim
			o.write(struct.pack('i',2)) # mode 2: 32bit float
			for start in origin:
				o.write(struct.pack('i', start)) # start point of map
			for dim in points:
				intervals = dim - 1
				o.write(struct.pack('i', intervals)) # number intervals per dim
			for scale in (high - low):
				o.write(struct.pack('f', scale)) # cell dim scales
			for i in range(3):
				o.write(struct.pack('f', 90.0)) # lattice angles
			for i in (1,2,3):
				o.write(struct.pack('i', i)) # axis mappings (fast x->y->z slow)
			for value in (np.min(mesh), np.max(mesh), np.mean(mesh)):
				o.write(struct.pack('f', value)) # map min/avg/max values
			o.write(struct.pack('i',1)) # space group
			for x in range(24,257):
				o.write(struct.pack('i',0))  # fill other fields with zero
			o.write(mesh.tobytes()) # write data

			print("{} mesh written".format(fileName))

	def isomap(self, protein=None, isoval=1.0, **kwargs):
		mesh, bounds = self.make_mesh(**kwargs)
		pcs_mesh = self.pcs_mesh(mesh)
		self.write_isomap(pcs_mesh, bounds)
		self.write_pymol_script(isoval=isoval, pdbFile=protein)

	def save(self, fileName='tensor.txt'):
		with open(fileName, 'w') as o:
			o.write(self.info(comment=False))


def make_tensor(x, y, z, axial, rhombic, 
	alpha, beta, gamma, lanthanide=None, temperature=298.15):
	"""
	Make a ChiTensor isntance from given parameters.
	This is designed to use pdb coordinates (x, y, z) and euler angles
	from an output like Numbat.

	Parameters
	----------
	x, y, z : floats
		tensor position in pdb coordiante in Angstroms
	axial, rhombic : floats
		the tensor anisotropies in units 10^-32
	alpha, beta, gamma : floats
		the euler angles in degrees that maps the tensor
		to the pdb (I think?)

	Returns
	-------
	ChiTensor : object :class:`paramagpy.metal.Metal`
		a tensor object for calulating paramagnetic effects on 
		nuclear spins in the pdb coordinate
	"""

	t = Metal()
	if lanthanide:
		t.set_lanthanide(lanthanide)
	t.position = np.array([x, y, z])*1E-10
	t.axrh = np.array([axial, rhombic])*1E-32
	t.eulers = np.array([alpha, beta, gamma])*(np.pi/180.)
	return t


def load_tensor(fileName):
	"""
	Load a metal object from file
	"""
	with open(fileName) as o:
		params = []
		for line in o:
			try:
				variableUnits, value = line.split(':')
				variable, units = variableUnits.split('|')
				params.append([variable.strip(), float(value)])
			except ValueError:
				print("WARNING: Line ignored reading tensor:\n{}".format(line))
		t = Metal()
		t.set_params(params)
	return t

