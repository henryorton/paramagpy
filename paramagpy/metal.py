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
	>>> print(euler_to_matrix(eulers))
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
	>>> print(matrix_to_euler(matrix))
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


def unique_eulers(a, b, g):
	if a < 0.0 and g < 0.0:
		alpha = a + np.pi
		beta  = np.pi - b
		gamma = -g
	elif a < 0.0:
		alpha = a + np.pi
		beta  = np.pi - b
		gamma = np.pi - g
	elif g < 0.0:
		alpha = a
		beta  = b
		gamma = g + np.pi
	else:
		alpha = a
		beta  = b
		gamma = g
	return np.array([alpha, beta, gamma])


class Metal(object):
	"""
	An object for paramagnetic chi tensors and delta-chi tensors.
	This must be created by specifying position, euler angles and 
	eigenvalues only.
	"""

	# Gyromagnetic ratio of an electron
	MU0 = 4*np.pi*1E-7
	MUB = 9.274E-24
	K = 1.381E-23
	HBAR = 1.0546E-34
	GAMMA = 1.760859644E11

	# Stored values get scaled by this amount for the fitting algorithm
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
	}

	# J, g, T1e values for lanthanide series
	lanth_lib = OrderedDict([('Zero',(0., 0. , 0.)),
		('Ce', ( 5./2., 6./7. , 0.133E-12)),
		('Pr', ( 4./1., 4./5. , 0.054E-12)),
		('Nd', ( 9./2., 8./11., 0.210E-12)),
		('Pm', ( 4./1., 3./5. , np.nan   )),
		('Sm', ( 5./2., 2./7. , 0.074E-12)),
		('Eu', ( 2./1., 3./2. , 0.015E-12)),
		('Gd', ( 7./2., 2./1. , np.nan   )),
		('Tb', ( 6./1., 3./2. , 0.251E-12)),
		('Dy', (15./2., 4./3. , 0.240E-12)),
		('Ho', ( 8./1., 5./4. , 0.209E-12)),
		('Er', (15./2., 6./5. , 0.189E-12)),
		('Tm', ( 6./1., 7./6. , 0.268E-12)),
		('Yb', ( 7./2., 8./7. , 0.157E-12))]
	)

	# Template anisotropies [axial, rhombic]
	lanth_axrh = OrderedDict([('Zero',(0. , 0.)),
		('Ce', (  2.1,  0.7)),
		('Pr', (  3.4,  2.1)),
		('Nd', (  1.7,  0.4)),
		('Pm', (  0.0,  0.0)),
		('Sm', (  0.2,  0.1)),
		('Eu', (  2.4,  1.5)),
		('Gd', (  0.0,  0.0)),
		('Tb', ( 42.1, 11.2)),
		('Dy', ( 34.7, 20.3)),
		('Ho', ( 18.5,  5.8)),
		('Er', (-11.6, -8.6)),
		('Tm', (-21.9,-20.1)),
		('Yb', ( -8.3, -5.8))]
	)

	upper_coords = ((0,1,0,0,1),(0,1,1,2,2))
	lower_coords = ((0,1,1,2,2),(0,1,0,0,1))

	@staticmethod
	def pack_tensor_params(metals):
		params = tuple(metals[0].position*1E10)
		for m in metals:
			params += tuple(m.upper_triang*1E32)
		return params

	@staticmethod
	def unpack_tensor_params(params):
		num = (len(params)-3)//5
		pos = np.array(params[:3])*1E-10
		tensor_params = [np.array(params[3+5*i:3+5+5*i])*1E-32 
			for i in range(num)]
		return pos, tensor_params

	@classmethod
	def anisotropy_to_eigenvalues(cls, axial, rhombic):
		"""
		Calculate [dx,dy,dz] eigenvalues from axial and rhombic
		tensor anisotropies (axial and rhombic parameters).
		Calculations assume traceless tensor.

		Parameters
		----------
		axial : float
			the axial anisotropy of the tensor

		rhombic : float
			the rhombic anisotropy of the tensor

		Returns
		-------
		euler_angles : array of floats
			the euler angles [alpha,beta,gamma] in radians
			by ZYZ convention
		"""
		dx =  rhombic/2. - axial/3.
		dy = -rhombic/2. - axial/3.
		dz = (axial*2.)/3.
		return np.array([dx,dy,dz])
	
	@classmethod
	def eigenvalues_to_anisotropy(cls, dx, dy, dz):
		"""
		Calculate axial and rhombic tensor anisotropies from 
		eigenvalues dx,dy,dz

		Parameters
		----------
		dx, dy, dz : floats
			the eigenvalues of the tensor.
			These are the principle axis magnitudes

		Returns
		-------
		axial, rhombic : tuple of floats
			the tensor anisotropies
		"""
		axial = dz - (dx + dy) / 2.
		rhombic = dx - dy
		return np.array([axial, rhombic])

	def __init__(self, position=(0,0,0), eulers=(0,0,0), 
		axrh=(0,0), mueff=0.0, shift=0.0, temperature=298.15, t1e=0.0,
		B0=18.79, taur=0.0):
		"""
		Instantiate ChiTensor object

		Parameters
		----------
		lanthanide : str, optional
			the lanthanide name e.g. 'Tb'
			if not given, default is None and the object may only
			be used for calculating PCS
		position : array, optional
			the (x,y,z) position in meters. Default is (0,0,0)
			stored as a np.matrix object.
		eulers : array, optional
			the euler angles [alpha,beta,gamma] in radians 
			by ZYZ convention. Defualt is (0,0,0)
		eigenvalues : array, optional
			the principal axes magnitudes [x,y,z]. Decaulf is (0,0,0).
			If traceless <eigenvalues> and the argument <lanthanide> are
			provided, isotropic component is calculated and added 
			automatically. 
		"""
		self.position = np.array(position, dtype=float)
		self.eulers = np.array(eulers, dtype=float)
		self.axrh = np.array(axrh, dtype=float)
		self.mueff = mueff
		self.shift = shift
		self.temperature = temperature
		self.t1e = t1e
		self.B0 = B0
		self.taur = taur

	def __str__(self):
		return str(self.tensor)

	def __abs__(self):
		return sum(self.eigenvalues)/3.

	@property
	def tauc(self):
		return 1./(1./self.taur + 1./self.t1e)

	def copy(self):
		return self.__class__(
			position=tuple(self.position), 
			eulers=tuple(self.eulers),
			axrh=tuple(self.axrh), 
			mueff=self.mueff, 
			shift=self.shift, 
			temperature=self.temperature, 
			t1e=self.t1e,
			B0=self.B0, 
			taur=self.taur)

	def set_lanthanide(self, lanthanide, set_dchi=True):
		J, g, t1e = self.lanth_lib[lanthanide]
		self.t1e = t1e
		self.set_Jg(J, g)
		if set_dchi:
			ax, rh = self.lanth_axrh[lanthanide]
			self.axrh = np.array([ax,rh])*1E-32

	def set_Jg(self, J, g):
		self.mueff = g * self.MUB * (J*(J+1))**0.5

	def info(self, comment=True):
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
		return i

	def get_params(self, params):
		pars = []
		for param in params:
			scale = self.fit_scaling.get(param, 1.0)
			pars.append(scale * getattr(self, param))
		return pars

	def set_params(self, paramValues):
		for param, value in paramValues:
			scale = self.fit_scaling.get(param, 1.0)
			setattr(self, param, value/scale)

	@property
	def x(self):
		return self.position[0]
	@x.setter
	def x(self, value):
		self.position[0] = value
	@property
	def y(self):
		return self.position[1]
	@y.setter
	def y(self, value):
		self.position[1] = value
	@property
	def z(self):
		return self.position[2]
	@z.setter
	def z(self, value):
		self.position[2] = value
	@property
	def a(self):
		return self.eulers[0]
	@a.setter
	def a(self, value):
		self.eulers[0] = value
	@property
	def b(self):
		return self.eulers[1]
	@b.setter
	def b(self, value):
		self.eulers[1] = value
	@property
	def g(self):
		return self.eulers[2]
	@g.setter
	def g(self, value):
		self.eulers[2] = value
	@property
	def ax(self):
		return self.axrh[0]
	@ax.setter
	def ax(self, value):
		self.axrh[0] = value
	@property
	def rh(self):
		return self.axrh[1]
	@rh.setter
	def rh(self, value):
		self.axrh[1] = value
	@property
	def iso(self):
		return self.isotropy
	@iso.setter
	def iso(self, value):
		self.isotropy = value
	@property
	def B0_MHz(self):
		return self.B0 * 42.57747892
	@B0_MHz.setter
	def B0_MHz(self, value):
		self.B0 = value / 42.57747892

	@property
	def eigenvalues(self):
		return self.anisotropy_to_eigenvalues(*self.axrh) + self.isotropy

	@eigenvalues.setter
	def eigenvalues(self, newEigenvalues):
		self.axrh = self.eigenvalues_to_anisotropy(*newEigenvalues)
		self.isotropy = np.round(np.sum(newEigenvalues)/3., 40)

	@property
	def isotropy(self):
		return (self.MU0 * self.mueff**2) / (3*self.K*self.temperature)

	@isotropy.setter
	def isotropy(self, newIsotropy):
		if newIsotropy<0:
			raise ValueError("A tensor with negative isotropy is not allowed")
		self.mueff = ((newIsotropy*3*self.K*self.temperature) / self.MU0)**0.5

	@property
	def rotationMatrix(self):
		return euler_to_matrix(self.eulers)

	@rotationMatrix.setter
	def rotationMatrix(self, newMatrix):
		self.eulers = unique_eulers(*matrix_to_eulers(newMatrix))

	@property
	def tensor(self):
		R = self.rotationMatrix
		return R.dot(np.diag(self.eigenvalues)).dot(R.T)

	@tensor.setter
	def tensor(self, newTensor):
		eigenvals, eigenvecs = np.linalg.eigh(newTensor)
		eigs = zip(eigenvals, np.array(eigenvecs).T)
		iso = np.sum(eigenvals)/3.
		eigenvals, (x, y, z) = zip(*sorted(eigs, key=lambda x: abs(x[0]-iso)))
		eigenvecs = x * z.dot(np.cross(x,y)), y, z
		rotationMatrix = np.vstack(eigenvecs).T
		eulers = unique_eulers(*matrix_to_euler(rotationMatrix))
		self.eulers = np.array(eulers, dtype=float)
		self.eigenvalues = eigenvals

	@property
	def tensor_traceless(self):
		return self.tensor - np.identity(3)*self.isotropy

	@property
	def saupe_factor(self):
		return (self.B0**2) / (5 * self.MU0 * self.K * self.temperature)

	@property
	def tensor_saupe(self):
		return self.saupe_factor * self.tensor_traceless

	@tensor_saupe.setter
	def tensor_saupe(self, new_saupe_tensor):
		old_mueff = self.mueff
		self.tensor = new_saupe_tensor / self.saupe_factor
		self.mueff = old_mueff

	@property
	def upper_triang(self):
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
	def upper_triang_saupe(self):
		return self.tensor_saupe[self.upper_coords]

	@upper_triang.setter
	def upper_triang_saupe(self, elements):
		newTensor = np.zeros(9).reshape(3,3)
		newTensor[self.upper_coords] = elements
		newTensor[self.lower_coords] = elements
		newTensor[2,2] = - elements[0] - elements[1]
		self.tensor_saupe = newTensor

	def set_utr(self):
		self.tensor = self.tensor

	def dipole_shift_tensor(self, position):
		"""
		Calculate the chemical shift tensor at the given postition due to the
		paramagnetic dipole tensor field

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
		"""
		val = self.dipole_shift_tensor(position).trace()/3.
		return 1E6*val + self.shift


	def atom_pcs(self, atom, racs=False, rads=False):
		value = self.pcs(atom.position)
		if racs:
			value += self.racs(atom.csa)
		if rads:
			value += self.rads(atom.position)
		return value


	def atom_set_position(self, atom):
		self.position = atom.position

	def fast_pcs(self, posarray):
		"""
		Rapidly calculatew the psuedo-contact shift at `n` positions.
		This efficient algorithm calculates the PCSs for an array of 
		positions and is best used where speed is required for fitting.

		Parameters
		----------
		posarray : array of positions with shape (n,3)
			the position (x, y, z) in meters

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
		given postition. The partial alignment induced by an anisotropic 
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
		ds = self.fast_dipole_shift_tensor(posarray)
		rads = np.einsum('jk,ikl->ijl',self.tensor_saupe,ds)
		return 1E6*rads.trace(axis1=1,axis2=2)/3.

	def racs(self, csa):
		"""
		Calculate the residual anisotropic chemical shift at the 
		given postition. The partial alignment induced by an anisotropic 
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
		racs = self.tensor_saupe.dot(csa).trace()/3.
		return 1E6*racs

	def fast_racs(self, csaarray):
		racs = np.einsum('jk,ikl->ijl',self.tensor_saupe,csaarray)
		return 1E6*racs.trace(axis1=1,axis2=2)/3.

	def fast_dsa_r1(self, posarray, gammaarray, csaarray=0.0):
		ds = self.fast_dipole_shift_tensor(posarray)
		sis_para = self.fast_second_invariant_squared(ds + csaarray)
		if isinstance(csaarray, np.ndarray):
			sis_dia = self.fast_second_invariant_squared(csaarray)
		else:
			sis_dia = 0.0
		sis_eff = sis_para - sis_dia
		omegas = self.B0 * gammaarray
		pf = (2./15.)*sis_eff*omegas**2
		rate = pf * self.spec_dens(self.taur, omegas)
		return rate

	def fast_dsa_r2(self, posarray, gammaarray, csaarray=0.0):
		ds = self.fast_dipole_shift_tensor(posarray)
		xx, yy, zz = np.linalg.eigvals(ds+csaarray).T
		sis_para = self.fast_second_invariant_squared(ds + csaarray)
		if isinstance(csaarray, np.ndarray):
			sis_dia = self.fast_second_invariant_squared(csaarray)
		else:
			sis_dia = 0.0
		sis_eff = sis_para - sis_dia
		omegas = self.B0 * gammaarray
		pf = (1./45.)*sis_eff*omegas**2
		rate = pf * (4*self.spec_dens(self.taur, 0.    ) +
			         3*self.spec_dens(self.taur, omegas))
		return rate

	def fast_sbm_r1(self, posarray, gammaarray):
		distance = np.linalg.norm(posarray-self.position, axis=1)
		p1 = self.MU0/(4.*np.pi)
		p2 = (gammaarray * self.mueff)/distance**3
		rate = (2./15.)*(p1*p2)**2 * (
			3*self.spec_dens(self.tauc, self.B0*gammaarray) + 
			7*self.spec_dens(self.tauc, self.B0*self.GAMMA))
		return rate

	def fast_sbm_r2(self, posarray, gammaarray):
		distance = np.linalg.norm(posarray-self.position, axis=1)
		p1 = self.MU0/(4.*np.pi)
		p2 = (gammaarray * self.mueff)/distance**3
		rate = (1./15.)*(p1*p2)**2 * (
			 4*self.spec_dens(self.tauc, 0.) +
			 3*self.spec_dens(self.tauc, self.B0*gammaarray) + 
			13*self.spec_dens(self.tauc, self.B0*self.GAMMA))
		return rate

	def fast_pre(self, posarray, gammaarray, rtype, 
		dsa=True, sbm=True, csaarray=0.0):
		rates = 0.0
		if rtype=='r1':
			if dsa:
				rates += self.fast_dsa_r1(posarray, gammaarray, csaarray)
			if sbm:
				rates += self.fast_sbm_r1(posarray, gammaarray)
		elif rtype=='r2':
			if dsa:
				rates += self.fast_dsa_r2(posarray, gammaarray, csaarray)
			if sbm:
				rates += self.fast_sbm_r2(posarray, gammaarray)
		return rates
			

		
	@staticmethod
	def spec_dens(tau, omega):
		"""
		Calculate spectral density at omega

		Parameters
		----------
		omega : float
			the spectral density argument

		Returns
		-------
		value : float
			the value of the spectral denstiy at <omega>
		"""
		return tau/(1+(omega*tau)**2)

	@staticmethod
	def second_invariant_squared(tensor):
		"""
		Calculate the second invariant at some position
		due to the magnetic susceptibility

		Parameters
		----------
		position : array floats
			the position (x, y, z) in meters

		Returns
		-------
		secondInvariant : float
			the second invariant of the shift tensor
		"""
		aniso = tensor - tensor.trace()*np.identity(3)/3.
		eigenvals, eigenvecs = np.linalg.eig(aniso)
		x, y, z = eigenvals
		secondInvariantSquared = x*x + y*y + z*z - x*y - x*z - y*z
		if secondInvariantSquared.imag>0:
			raise ValueError("Imaginary second invariant")
		return secondInvariantSquared.real

	@staticmethod
	def fast_second_invariant_squared(tensorarray):
		"""
		Calculate the second invariant at some position
		due to the magnetic susceptibility

		Parameters
		----------
		position : array floats
			the position (x, y, z) in meters

		Returns
		-------
		secondInvariant : float
			the second invariant of the shift tensor
		"""
		xx, yy, zz = np.linalg.eigvals(tensorarray).real.T
		return xx*xx + yy*yy + zz*zz - xx*yy - xx*zz - yy*zz

	def dsa_r1(self, position, gamma, csa=0.0, ignorePara=False):
		"""
		Calculate R1 relaxation due to Curie Spin

		Parameters
		----------
		atom

		Returns
		-------
		value : float
			The R1 relaxation rate in /s
		"""
		# if method is 'dsaiso':
		# 	distance = np.linalg.norm(position-self.metal.position)
		# 	sec_invar_sq = self.metal.isotropic_second_invariant(distance)**2
		# sec_invar_sq = self.parse_curie_method(position, method, csa)
		omega = self.B0 * gamma
		if ignorePara:
			shieldingTensor = csa
		else:
			shieldingTensor = self.dipole_shift_tensor(position) + csa
		secondInvariantSquared = self.second_invariant_squared(shieldingTensor)
		pf = (2./15.)*secondInvariantSquared*omega**2
		rate = pf * self.spec_dens(self.taur, omega)
		return rate

	def dsa_r2(self, position, gamma, csa=0.0, ignorePara=False):
		"""
		Calculate R2 relaxation due to Curie Spin

		Parameters
		----------
		atom : 


		Returns
		-------
		value : float
			The R2 relaxation rate in /s
		"""

		# if method is 'dsaiso':
		# 	distance = np.linalg.norm(position-self.metal.position)
		# 	sec_invar_sq = self.metal.isotropic_second_invariant(distance)**2

		# sec_invar_sq = self.parse_curie_method(position, method, csa)
		omega = self.B0 * gamma
		if ignorePara:
			shieldingTensor = csa
		else:
			shieldingTensor = self.dipole_shift_tensor(position) + csa
		secondInvariantSquared = self.second_invariant_squared(shieldingTensor)
		pf = (1./45.)*secondInvariantSquared*omega**2
		rate = pf * (4*self.spec_dens(self.taur, 0.   ) +
			         3*self.spec_dens(self.taur, omega))
		return rate


	def sbm_r1(self, position, gamma):
		distance = np.linalg.norm(position-self.position)
		p1 = self.MU0/(4.*np.pi)
		p2 = (gamma * self.mueff)/distance**3
		rate = (2./15.)*(p1*p2)**2 * (
			3*self.spec_dens(self.tauc, self.B0*gamma) + 
			7*self.spec_dens(self.tauc, self.B0*self.GAMMA))
		return rate


	def sbm_r2(self, position, gamma):
		distance = np.linalg.norm(position-self.position)
		p1 = self.MU0/(4.*np.pi)
		p2 = (gamma * self.mueff)/distance**3
		rate = (1./15.)*(p1*p2)**2 * (
			 4*self.spec_dens(self.tauc, 0.) +
			 3*self.spec_dens(self.tauc, self.B0*gamma) + 
			13*self.spec_dens(self.tauc, self.B0*self.GAMMA))
		return rate

	def rdc(self, vector, gam1, gam2):
		"""
		Calculate Residual Dipolar Coupling (RDC)

		Parameters
		----------
		vector : [x,y,z] array of float
			internuclear vector in meters
		gam1 : float
			gyromagnetic ratio of spin 1 in rad/s/T
		gam2 : float
			gyromagnetic ratio of spin 2 in rad/s/T

		Returns
		-------
		rdc : float
			the RDC in Hz
		"""
		dist = np.linalg.norm(vector)
		pf = -(self.MU0 * gam1 * gam2 * self.HBAR) / (4 * np.pi * dist**5)
		rdc_radians = 3*pf * vector.dot(self.tensor_saupe).dot(vector) 
		return rdc_radians / (2*np.pi)


	def make_mesh(self, density=2, size=40.0):
		origin = np.asarray(density * (self.position*1E10 - size/2.0), dtype=int)
		low = origin / float(density)
		high = low + size
		points = np.array([int(density*size)]*3) + 1
		domains = [1E-10*np.linspace(*i) for i in zip(low, high, points)]
		mesh = np.array(np.meshgrid(*domains, indexing='ij')).T
		return mesh, (origin, low, high, points)

	def pcs_mesh(self, mesh):
		og_shape = mesh.shape[:3]
		pcs_mesh = self.fast_pcs(mesh.reshape(np.prod(og_shape),3))
		return pcs_mesh.reshape(*og_shape)

	def write_pymol_script(self, isoval=1.0, surfaceName='isomap', 
		scriptName='isomap.pml', meshName='./isomap.pml.ccp4', pdbFile=None):
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
			s += "cmd.load(os.path.join(curdir, '{}'),'{}')\n".format(pdbFile, protName)
			s += "show_as cartoon, {}\n".format(protName)
		with open(scriptName, 'w') as o:
			o.write(s)
			print("{} script written".format(scriptName))

	def write_isomap(self, mesh, bounds, fileName='isomap.pml.ccp4'):
		mesh = np.asarray(mesh, np.float32)
		origin, low, high, points = bounds
		with open(fileName, 'wb') as o:
			for dim in points:
				o.write(struct.pack('i', dim)) # number points / dim
			o.write(struct.pack('i',2)) # mode 2: 32bit float
			for start in origin:
				o.write(struct.pack('i', start)) # start point of map
			for dim in points:
				intervals = dim - 1
				o.write(struct.pack('i', intervals)) # number intervals / dim
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
	ChiTensor : object
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
	Docstring for load tensor
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

