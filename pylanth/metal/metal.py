import numpy as np
import struct



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
	matrix : numpy 3x3 matrix object
		the rotation matrix
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
	M : 3x3 array
		a rotation matrix

	Returns
	-------
	eulers : array of floats
		the euler angles [alpha,beta,gamma] in radians
		by ZYZ convention
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

	# J, g, T1e values for lanthanide series
	lanth_lib = {
		'Ce': ( 5./2., 6./7. , 0.133E-12),
		'Pr': ( 4./1., 4./5. , 0.054E-12),
		'Nd': ( 9./2., 8./11., 0.210E-12),
		'Pm': ( 4./1., 3./5. , None    ),
		'Sm': ( 5./2., 2./7. , 0.074E-12),
		'Eu': ( 2./1., 3./2. , 0.015E-12),
		'Gd': ( 7./2., 2./1. , None    ),
		'Tb': ( 6./1., 3./2. , 0.251E-12),
		'Dy': (15./2., 4./3. , 0.240E-12),
		'Ho': ( 8./1., 5./4. , 0.209E-12),
		'Er': (15./2., 6./5. , 0.189E-12),
		'Tm': ( 6./1., 7./6. , 0.268E-12),
		'Yb': ( 7./2., 8./7. , 0.157E-12)
	}

	upper_coords = [(0,1,0,0,1),(0,1,1,2,2)]
	lower_coords = [(0,1,1,2,2),(0,1,0,0,1)]

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
		tensor_params = [np.array(params[3+5*i:3+5+5*i])*1E-32 for i in range(num)]
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
		sortedEigs = cls.unique_eigenvalues(dx, dy, dz)
		return np.array(sortedEigs)
	
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
		dx, dy, dz = cls.unique_eigenvalues(dx, dy, dz)
		axial = dz - (dx + dy) / 2.
		rhombic = dx - dy
		return np.array([axial, rhombic])

	@classmethod
	def unique_anisotropy(cls, ax, rh):
		unsorteigs = cls.anisotropy_to_eigenvalues(ax, rh)
		eigenvalues = sorted(unsorteigs, key = lambda x: abs(x))
		return cls.eigenvalues_to_anisotropy(*eigenvalues)

	@classmethod
	def unique_eigenvalues(cls, dx, dy, dz):
		eigs = np.array([dx,dy,dz], dtype=float)
		iso = np.sum(eigs)/3.
		eigenvalues = sorted(eigs, key = lambda x: abs(x-iso))
		return eigenvalues

	@classmethod
	def make_tensor(cls, x, y, z, axial, rhombic, 
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
		if lanthanide:
			J, g, t1e = cls.lanth_lib[lanthanide]
			mueff = g * cls.MUB * (J*(J+1))**0.5
		else:
			mueff = 0.0
			t1e = None

		position = np.array([x, y, z])*1E-10
		axrh = cls.unique_anisotropy(axial, rhombic)*1E-32
		eulers = np.array([alpha, beta, gamma])*(np.pi/180.)
		return cls(position, eulers, axrh, mueff, temperature, t1e)
		

	def __init__(self, position=(0,0,0), eulers=(0,0,0), 
		axrh=(0,0), mueff=0.0, temperature=298.15, t1e=None,
		B0=None, taur=None):
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
		return self.__class__(tuple(self.position), tuple(self.eulers), 
			tuple(self.axrh), self.mueff, self.temperature, self.t1e)

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
		i += l.format('temp','K',self.temperature)
		if self.t1e:
			i += l.format('t1e','ps',self.t1e*1E12)
		else:
			i += l.format('t1e','ps',np.nan)
		return i

	@property
	def eigenvalues(self):
		return self.anisotropy_to_eigenvalues(*self.axrh) + self.isotropy

	@eigenvalues.setter
	def eigenvalues(self, newEigenvalues):
		isotropy = np.sum(newEigenvalues)/3.
		self.isotropy = isotropy
		self.axrh = self.eigenvalues_to_anisotropy(*newEigenvalues)

	@property
	def isotropy(self):
		return (self.MU0 * self.mueff**2) / (3*self.K*self.temperature)

	@isotropy.setter
	def isotropy(self, newIsotropy):
		newIsotropy = 1E-32*np.round(newIsotropy*1E32, 10)
		self.mueff = ((newIsotropy*3*self.K*self.temperature) / self.MU0)**0.5

	@property
	def rotationMatrix(self):
		return euler_to_matrix(self.eulers)

	@rotationMatrix.setter
	def rotationMatrix(self, newMatrix):
		self.eulers = matrix_to_eulers(newMatrix)

	@property
	def tensor(self):
		R = self.rotationMatrix
		return R.dot(np.diag(self.eigenvalues)).dot(R.T)

	@tensor.setter
	def tensor(self, newTensor):
		eigenvals, eigenvecs = np.linalg.eig(newTensor)
		eigs = zip(eigenvals, np.array(eigenvecs).T)
		iso = np.sum(eigenvals)/3.
		eigenvals, eigenvecs = zip(*sorted(eigs, key=lambda x: abs(x[0]-iso)))
		rotationMatrix = np.vstack(eigenvecs).T
		eulers = matrix_to_euler(rotationMatrix)
		self.eulers = np.array(eulers, dtype=float)
		self.eigenvalues = eigenvals

	@property
	def tensor_traceless(self):
		return self.tensor - np.identity(3)*self.isotropy

	@property
	def tensor_saupe(self):
		pf = (self.B0**2) / (15*self.MU0 * self.K * self.temperature)
		return pf * self.tensor_traceless

	@tensor_saupe.setter
	def tensor_saupe(self, new_saupe_tensor):
		pf =  (15*self.MU0 * self.K * self.temperature) / (self.B0**2)
		self.tensor = pf * new_saupe_tensor

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
			the peudo-contact shift in parts-per-million (ppm)
		"""
		return 1E6*self.dipole_shift_tensor(position).trace()/3.

	def pcs2(self, position):
		ds = self.dipole_shift_tensor(position)
		align = (2./3.)*self.tensor_saupe.dot(ds).trace()
		iso = ds.trace()/3.
		return 1E6*(align+iso)

	def fast_pcs(self, posarray):
		pos = posarray - self.position
		dist = np.linalg.norm(pos, axis=1)
		dot1 = np.einsum('ij,jk->ik', pos, self.tensor_traceless)
		dot2 = np.einsum('ij,ij->i', pos, dot1)
		return 1E6*(1./(4.*np.pi))*(dot2/dist**5)

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




	def curie_r1(self, position, gamma, csa=0.0, ignorePara=False):
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

	def curie_r2(self, position, gamma, csa=0.0, ignorePara=False):
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
		things
		tensor : ChiTensor object
			a paramagnetic tensor object from which
			<delta_tensor> 3x3 traceless matrix attribute must be present

		Returns
		-------
		value : float
			the RDC in Hz
		"""
		distance = np.linalg.norm(vector)
		numer = -self.HBAR * self.B0**2 * gam1 * gam2
		denom = 120. * self.K * self.temperature * np.pi**2
		preFactor = numer/denom
		p1 = (1./distance**5)*np.kron(vector,vector).reshape(3,3)
		p2 = (1./distance**3)*np.identity(3)
		return preFactor * ((3.*p1 - p2).dot(self.tensor_traceless)).trace()

	def rdc2(self, vector, gam1, gam2):
		dist = np.linalg.norm(vector)
		pf = -(self.MU0 * gam1 * gam2 * self.HBAR) / (4 * np.pi * dist**5)
		dd = 3*np.kron(vector,vector).reshape(3,3) - np.identity(3)*dist**2
		return pf * (dd.dot(self.tensor_saupe)).trace() / (2*np.pi)
















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
		pcs_mesh = self.fast_pcs(mesh.reshape(np.prod(og_shape),3)).reshape(*og_shape)
		return pcs_mesh

	def write_pymol_script(self, protein=None, isoval=1.0, fileName='isomap.pml'):
		s = "# PyMOL macro for loading tensor isosurface from pyrelax\n"
		s += self.info()+'\n'
		s += "set normalize_ccp4_maps, off\n"
		s += "load ./isomap.pml.ccp4, isomap, 1, ccp4\n"
		s += "isosurface isoPos, isomap, {}\n".format(isoval)
		s += "isosurface isoNeg, isomap, -{}\n".format(isoval)
		s += "set transparency, 0.5, isoPos\n"
		s += "set transparency, 0.5, isoNeg\n"
		s += "set surface_color, blue, isoPos\n"
		s += "set surface_color, red, isoNeg\n"
		s += "pseudoatom isoOrig, pos=[{},{},{}]\n".format(*self.position*1E10)
		s += "show spheres, isoOrig\n"
		s += "color pink, isoOrig\n"
		if protein:
			s += "load {}\n".format(protein.name)
			s += "show_as cartoon, {}".format("".join(protein.name.split('.')[:-1]))
		with open(fileName, 'w') as o:
			o.write(s)


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

		return None


	def isomap(self, protein=None, isoval=1.0, **kwargs):
		mesh, bounds = self.make_mesh(**kwargs)
		pcs_mesh = self.pcs_mesh(mesh)
		self.write_isomap(pcs_mesh, bounds)
		self.write_pymol_script(protein, isoval)




# def simple_rdc(vector, gamma1, gamma2, B0, ax, rh, k, T, hbar):
# 	dist = np.linalg.norm(vector)
# 	x, y, z = vector
# 	pf = -(1./(4*np.pi)) * ((B0**2)/(15*k*T)) * ((gamma1*gamma2*hbar)/(2*np.pi*dist**3))
# 	val = ax * ((2*z**2-x**2-y**2)/(dist**2)) + (3./2.)*rh * ((x**2-y**2)/(dist**2))
# 	return pf*val



# x = 0.0E-10
# y = 0.0E-10
# z = 1.0E-10
# vec = np.array([x,y,z])

# gamma = 2*np.pi*42.576E6

# eulers = np.random.random(size=3)
# # eulers = np.array([0,0,0.])
# rot = euler_to_matrix(eulers)

# ax = -5.385
# rh = -2.365
# a, b, g = eulers*(180./np.pi)

# t = Metal.make_tensor(0,0,0,ax,rh,a,b,g,'Er')
# t.B0 = 14.1


# print(t.rdc(vec, gamma, gamma))
# print(t.rdc2(vec, gamma, gamma))
# ax, rh = t.axrh
# vecr = rot.T.dot(vec)
# print(simple_rdc(vecr, gamma, gamma, t.B0, ax, rh, t.K, t.temperature, t.HBAR))





	# # def second_invariant(self, position, add=0.0):
	# # 	"""
	# # 	Calculate the second invariant at some position
	# # 	due to the magnetic susceptibility

	# # 	Parameters
	# # 	----------
	# # 	position : array floats
	# # 		the position (x, y, z) in meters

	# # 	Returns
	# # 	-------
	# # 	secondInvariant : float
	# # 		the second invariant of the shift tensor
	# # 	"""
	# # 	ds = self.dipole_shift_tensor(position) + add
	# # 	dsa = ds - ds.trace()*np.identity(3)/3.
	# # 	eigenvals, eigenvecs = np.linalg.eig(dsa)
	# # 	x, y, z = eigenvals
	# # 	secondInvariantSquared = x*x + y*y + z*z - x*y - x*z - y*z
	# # 	return abs(secondInvariantSquared)**0.5


	# def isotropic_second_invariant(self, distance):
	# 	"""
	# 	Calculate the second invariant at some position
	# 	due to an isotropic magnetic susceptibility

	# 	Parameters
	# 	----------
	# 	distance : float
	# 		the distance in meters from the paramagnetic centre

	# 	Returns
	# 	-------
	# 	value : float
	# 		the isotropic second invariant
	# 	"""
	# 	return (3. * sum(self.eigenvalues)/3.) / (4. * np.pi * distance**3)
















