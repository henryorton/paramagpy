from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom, DisorderedAtom
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.Polypeptide import standard_aa_names
import numpy as np



def rotation_matrix(axis, theta):
	"""Return the rotation matrix associated with counterclockwise 
	rotation about the given axis by theta radians.

	Parameters
	----------
	axis : array of floats
		the [x,y,z] axis for rotation.

	Returns
	-------
	matrix : numpy 3x3 matrix object
		the rotation matrix
	"""
	axis = np.array(axis)
	axis /= np.linalg.norm(axis)
	a = np.cos(theta/2.0)
	b, c, d = -axis*np.sin(theta/2.0)
	aa, bb, cc, dd = a*a, b*b, c*c, d*d
	bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
	return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
					 [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
					 [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


class CustomAtom(Atom):

	MU0 = 4*np.pi*1E-7
	HBAR = 1.0546E-34

	gyro_lib = {
		'H': 2*np.pi*42.576E6,
		'N': 2*np.pi*-4.316E6,
		'C': 2*np.pi*10.705E6} # rad/s/T

	csa_lib = {
		'H': (np.array([-5.8 , 0.0  ,5.8 ])*1E-6, 8. *(np.pi/180.)),
		'N': (np.array([-62.8,-45.7 ,108.5])*1E-6, 19.*(np.pi/180.)),
		'C': (np.array([-86.5 ,11.8, 74.7])*1E-6, 38.*(np.pi/180.))}

	"""docstring for CustomAtom"""
	def __init__(self, *arg, **kwargs):
		super().__init__(*arg, **kwargs)
		self.coord = np.asarray(self.coord, dtype=np.float64)
		self.gamma = self.gyro_lib.get(self.element, 0.0)
		self._csa = None

	def __repr__(self):
		return "<Atom {0:3d}-{1:}>".format(self.parent.id[1], self.name)

	def top(self):
		return self.parent.parent.parent.parent

	@property
	def position(self):
		return self.coord*1E-10

	@position.setter
	def position(self, value):
		self.coord = value*1E10

	@property
	def csa(self):
		"""
		Get the CSA tensor at the nuclear position
		This uses the geometry of neighbouring atoms
		and a standard library from Bax J. Am. Chem. Soc. 2000

		Returns
		-------
		matrix : 3x3 array
			the CSA tensor in the PDB frame
			if appropriate nuclear positions are not
			available <None> is returned.
		"""

		if self._csa is not None:
			return self._csa

		def norm(x):
			return x/np.linalg.norm(x)

		res = self.parent
		resid = res.id
		respid = resid[0], resid[1]-1, resid[2]
		resnid = resid[0], resid[1]+1, resid[2]
		resp = res.parent.child_dict.get(respid)
		resn = res.parent.child_dict.get(resnid)

		pas, beta = self.csa_lib.get(self.name, (None,None))
		if resp:
			Hcond = self.element=='H', 'N' in res ,'C' in resp, beta
			Ncond = self.element=='N', 'H' in res ,'C' in resp, beta
		else:
			Hcond = (None,)
			Ncond = (None,)
		if resn:
			Ccond = self.element=='C', 'H' in resn,'N' in resn, beta
		else:
			Ccond = (None,)

		if all(Hcond):
			NC_vec = resp['C'].coord - res['N'].coord
			NH_vec =  res['H'].coord - res['N'].coord
			z = norm(np.cross(NC_vec, NH_vec))
			R = rotation_matrix(-z,beta)
			x = norm(R.dot(NH_vec))
			y = norm(np.cross(z, x))

		elif all(Ncond):
			NC_vec = resp['C'].coord - res['N'].coord
			NH_vec =  res['H'].coord - res['N'].coord
			y = norm(np.cross(NC_vec, NH_vec))
			R = rotation_matrix(-y,beta)
			z = norm(R.dot(NH_vec))
			x = norm(np.cross(y, z))

		elif all(Ccond):
			CN_vec = resn['N'].coord -  res['C'].coord
			NH_vec = resn['H'].coord - resn['N'].coord
			x = norm(np.cross(NH_vec, CN_vec))
			R = rotation_matrix(x,beta)
			z = norm(R.dot(CN_vec))
			y = norm(np.cross(z, x))

		else:
			return np.zeros(9).reshape(3,3)
		transform = np.vstack([x, y, z]).T
		tensor = transform.dot(np.diag(pas)).dot(transform.T)
		return tensor

	@csa.setter
	def csa(self, newTensor):
		if newTensor is None:
			self._csa = None
			return
		try:
			assert newTensor.shape == (3,3)
		except (AttributeError, AssertionError):
			print("The specified CSA tensor does not have the correct format")
			raise
		self._csa = newTensor


	def dipole_shift_tensor(self, position):
		pos = np.array(position, dtype=float) - self.position
		distance = np.linalg.norm(pos)
		preFactor = (self.MU0 * self.gamma * self.HBAR * 0.5) / (4.*np.pi)
		p1 = (1./distance**5)*np.kron(pos,pos).reshape(3,3)
		p2 = (1./distance**3)*np.identity(3)
		return (preFactor * (3.*p1 - p2))


# a = CustomAtom('H',[0,0,0],0,0,0,'H',0)

# pos = np.array([0,0,1.])*1E-10
# field = a.dipole_shift_tensor(pos).dot(np.array([0.,0.,1.]))
# omega = field * CustomAtom.gyro_lib['H']
# print(a.dipole_shift_tensor(pos)*1E6)



class CustomStructure(Structure):
	"""docstring for CustomStructure"""
	def __init__(self, *arg, **kwargs):
		super().__init__(*arg, **kwargs)

	# def select_models(self, models):
	# 	new = self.copy()
	# 	if type(models)!=list:
	# 		models = [models]
	# 	for i in list(new.child_dict):
	# 		if i not in models:
	# 			new.detach_child(i)
	# 	return new

	def parse(self, dataValues, models=None):
		used = set([])
		data = []

		if type(models)==int:
			chains = self[models].get_chains()
		elif type(models) in (list, tuple):
			chains = []
			for m in models:
				chains += self[m].get_chains()
		else:
			chains = self.get_chains()

		if dataValues.dtype in ('PCS', 'PRE'):
			for chain in chains:
				for key in dataValues:
					seq, name = key
					if seq in chain:
						resi = chain[seq]
						if name in resi:
							atom = resi[name]
							data.append((atom, *dataValues[key]))
							used.add(key)

		elif dataValues.dtype in ('RDC', 'CCR'):
			for chain in chains:
				for key in dataValues:
					(seq1, name1), (seq2, name2) = key
					if seq1 in chain and seq2 in chain:
						resi1 = chain[seq1]
						resi2 = chain[seq2]
						if name1 in resi1 and name2 in resi2:
							atom1 = resi1[name1]
							atom2 = resi2[name2]
							data.append((atom1, atom2, *dataValues[key]))
							used.add(key)

		unused = set(dataValues) - used
		if unused:
			message = "WARNING: Some values were not parsed to {}:"
			print(message.format(self.id))
			print(list(unused))
		return data
			


class CustomStructureBuilder(StructureBuilder):
	"""docstring for CustomStructureBuilder"""
	def __init__(self, *arg, **kwargs):
		super().__init__(*arg, **kwargs)

	def init_structure(self, structure_id):
		self.structure = CustomStructure(structure_id)

	def init_atom(self, name, coord, b_factor, occupancy, altloc, fullname,
				  serial_number=None, element=None):
		"""Create a new Atom object.
		Arguments:
		 - name - string, atom name, e.g. CA, spaces should be stripped
		 - coord - Numeric array (Float0, size 3), atomic coordinates
		 - b_factor - float, B factor
		 - occupancy - float
		 - altloc - string, alternative location specifier
		 - fullname - string, atom name including spaces, e.g. " CA "
		 - element - string, upper case, e.g. "HG" for mercury
		"""
		residue = self.residue
		# if residue is None, an exception was generated during
		# the construction of the residue
		if residue is None:
			return
		# First check if this atom is already present in the residue.
		# If it is, it might be due to the fact that the two atoms have atom
		# names that differ only in spaces (e.g. "CA.." and ".CA.",
		# where the dots are spaces). If that is so, use all spaces
		# in the atom name of the current atom.
		if residue.has_id(name):
			duplicate_atom = residue[name]
			# atom name with spaces of duplicate atom
			duplicate_fullname = duplicate_atom.get_fullname()
			if duplicate_fullname != fullname:
				# name of current atom now includes spaces
				name = fullname
				warnings.warn("Atom names %r and %r differ "
							  "only in spaces at line %i."
							  % (duplicate_fullname, fullname,
								 self.line_counter),
							  PDBConstructionWarning)
		self.atom = CustomAtom(name, coord, b_factor, occupancy, altloc,
						 fullname, serial_number, element)
		if altloc != " ":
			# The atom is disordered
			if residue.has_id(name):
				# Residue already contains this atom
				duplicate_atom = residue[name]
				if duplicate_atom.is_disordered() == 2:
					duplicate_atom.disordered_add(self.atom)
				else:
					# This is an error in the PDB file:
					# a disordered atom is found with a blank altloc
					# Detach the duplicate atom, and put it in a
					# DisorderedAtom object together with the current
					# atom.
					residue.detach_child(name)
					disordered_atom = DisorderedAtom(name)
					residue.add(disordered_atom)
					disordered_atom.disordered_add(self.atom)
					disordered_atom.disordered_add(duplicate_atom)
					residue.flag_disordered()
					warnings.warn("WARNING: disordered atom found "
								  "with blank altloc before line %i.\n"
								  % self.line_counter,
								  PDBConstructionWarning)
			else:
				# The residue does not contain this disordered atom
				# so we create a new one.
				disordered_atom = DisorderedAtom(name)
				residue.add(disordered_atom)
				# Add the real atom to the disordered atom, and the
				# disordered atom to the residue
				disordered_atom.disordered_add(self.atom)
				residue.flag_disordered()
		else:
			# The atom is not disordered
			residue.add(self.atom)




def load_pdb(fileName, ident=None):
	if not ident:
		ident = fileName
	parser = PDBParser(structure_builder=CustomStructureBuilder())
	return parser.get_structure(ident, fileName)
















