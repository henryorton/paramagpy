from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom, DisorderedAtom
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.Polypeptide import standard_aa_names
import numpy as np

# Datatypes for numpy structured arrays after parsing 
EXPDATA_DTYPE = {
	'PCS':np.dtype([
				('mdl', int   ),
				('atm', object),
				('exp', float ),
				('cal', float ),
				('err', float ),
				('idx', int   )]),
	'RDC':np.dtype([
				('mdl',  int   ),
				('atm',  object),
				('atx',  object),
				('exp',  float ),
				('cal',  float ),
				('err',  float ),
				('idx',  int   )]),
	'PRE':np.dtype([
				('mdl', int   ),
				('atm', object),
				('exp', float ),
				('cal', float ),
				('err', float ),
				('idx', int   )]),
	'CCR':np.dtype([
				('mdl',  int   ),
				('atm',  object),
				('atx',  object),
				('exp',  float ),
				('cal',  float ),
				('err',  float ),
				('idx',  int   )])}

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
	For example: [7,4,7,9,9,10,1] -> [2,1,2,3,3,4,0]

	Parameters
	----------
	indices : array-like integers
		a list of integers

	Returns
	-------
	new_indices : array-like integers
		the mapped integers with smallest size
	"""
	_, new_indices = np.unique(indices, return_inverse=True)
	return new_indices


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
		if self.parent:
			parentID = "{:3d}".format(self.parent.id[1])
		else:
			parentID = ""
		return "<Atom {0:}-{1:}>".format(parentID, self.name)

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
		"""
		Calculate the magnetic field shielding tensor at the given postition
		due to the nuclear dipole

		Assumes nuclear spin 1/2

		Parameters
		----------
		position : array floats
			the position (x, y, z) in meters

		Returns
		-------
		dipole_shielding_tensor : 3x3 array
			the tensor describing magnetic shielding at the given position
		"""
		pos = np.array(position, dtype=float) - self.position
		distance = np.linalg.norm(pos)
		preFactor = (self.MU0 * self.gamma * self.HBAR * 0.5) / (4.*np.pi)
		p1 = (1./distance**5)*np.kron(pos,pos).reshape(3,3)
		p2 = (1./distance**3)*np.identity(3)
		return (preFactor * (3.*p1 - p2))


class CustomStructure(Structure):
	"""This is an overload hack of the BioPython Structure object"""
	def __init__(self, *arg, **kwargs):
		super().__init__(*arg, **kwargs)

	def parse(self, dataValues, models=None):
		"""
		Associate experimental data with atoms of the PDB file
		This method takes a DataContainer instance from the 
		dataparse module

		Parameters
		----------
		dataValues : DataContainer instance
			a dictionary containing the experimental values

		Returns
		-------
		dataArray : numpy structured array
			the returned array has a row for each relevant atom
			in the PDB file. The columns contain model,
			experimental/calculated data, errors and indexes.
		"""
		if type(models)==int:
			mods = [self[models]]
		elif type(models) in (list, tuple):
			mods = (self[m] for m in models)
		else:
			mods = self.get_models()

		data = []
		used = set([])

		if dataValues.dtype in ('PCS', 'PRE'):
			for m in mods:
				for chain in m:
					for key in dataValues:
						seq, name = key
						if seq in chain:
							resi = chain[seq]
							if name in resi:
								a = resi[name]
								exp, err = dataValues[key]
								idx = a.serial_number
								tmp = (m.id, a, exp, np.nan, err, idx)
								data.append(tmp)
								used.add(key)

		if dataValues.dtype in ('RDC', 'CCR'):
			for m in mods:
				for chain in m:
					for key in dataValues:
						(seq1, name1), (seq2, name2) = key
						if seq1 in chain and seq2 in chain:
							resi1 = chain[seq1]
							resi2 = chain[seq2]
							if name1 in resi1 and name2 in resi2:
								a1 = resi1[name1]
								a2 = resi2[name2]
								exp, err = dataValues[key]
								idx1 = a1.serial_number
								idx2 = a2.serial_number
								if dataValues.dtype=='RDC':
									idx = unique_pairing(idx1, idx2)
								elif dataValues.dtype=='CCR':
									idx = cantor_pairing(idx1, idx2)
								tmp = (m.id, a1, a2, exp, np.nan, err, idx)
								data.append(tmp)
								used.add(key)

		unused = set(dataValues) - used
		if unused:
			message = "WARNING: Some values were not parsed to {}:"
			print(message.format(self.id))
			print(list(unused))
		arr = np.array(data, dtype=EXPDATA_DTYPE.get(dataValues.dtype))
		arr['idx'] = clean_indices(arr['idx'])
		return arr
			


class CustomStructureBuilder(StructureBuilder):
	"""This is an overload hack of BioPython's CustomStructureBuilder"""
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
	"""
	Read PDB from file into biopython structure object

	Parameters
	----------
	fileName : str
		the path to the file
	ident : str (optional)
		the desired identity of the structure object

	Returns
	-------
	values : :class:`paramagpy.protein.CustomStructure`
		a structure object containing the atomic coordinates
	"""
	if not ident:
		ident = fileName
	parser = PDBParser(structure_builder=CustomStructureBuilder())
	return parser.get_structure(ident, fileName)



class PyMolScript(object):
	"""
	A PyMol helper class for constructing a PyMol script
	that will load PDB files and density map files
	"""
	def __init__(self):
		self.preamble = "# PyMOL script generated by paramagpy\n"
		self.preamble += "import os, pymol\n"
		self.preamble += "curdir = os.path.dirname(pymol.__script__)\n"
		self.preamble += "set normalize_ccp4_maps, off\n"

		self.lines = []

	def __add__(self, value):
		self.lines += [value]
		return self

	def add_atom(self, position, name, colour='pink', size=1.0, label=True):
		self += "pseudoatom {name}, pos=[{:.4f}, {:.4f}, {:.4f}]".format(
			name=name, *position)
		self += "show spheres, {name}".format(name=name)
		self += "color {colour}, {name}".format(colour=colour, name=name)
		self += "set sphere_scale, {size}, {name}".format(size=size, name=name)
		if label:
			self += "label {name}, '{name}'".format(name=name)

	def add_pdb(self, path, name=None, showAs="cartoon"):
		if name is None:
			name = path
		self += "cmd.load(os.path.join(curdir, '{path}'),'{name}')".format(
			path=path, name=name)
		self += "show_as {show_as}, {name}".format(
			show_as=showAs, name=name)

	def add_map(self, path, name, isoVals, colours, 
		surfaceType='isosurface', transparency=0.5):
		if surfaceType == 'isosurface':
			surfaceColour = 'surface_color'
		elif surfaceType == 'isomesh':
			surfaceColour = 'mesh_color'
		elif surfaceType == 'isodot':
			surfaceColour = 'dot_color'
		else:
			surfaceColour = ''

		self += "cmd.load(os.path.join(curdir, '{mapPath}'), '{mapName}', 1, 'ccp4')".format(
			mapPath=path, mapName=name)

		for isoVal, colour in zip(isoVals, colours):
			self += "{surfaceType} {mapName}{isoValue:+}, {mapName}, {isoValue}".format(
				surfaceType=surfaceType, mapName=name, isoValue=isoVal)
			self += "set {surfaceColour}, {colour}, {mapName}{isoValue:+}".format(
				surfaceColour=surfaceColour, colour=colour, mapName=name, isoValue=isoVal)
			self += "set transparency, {transparency}, {mapName}{isoValue:+}".format(
				transparency=transparency, mapName=name, isoValue=isoVal)

	def write(self, fileName):
		with open(fileName, 'w') as f:
			f.write(self.preamble)
			f.write('\n'.join(self.lines))



