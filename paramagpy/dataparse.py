import numpy as np
from collections import OrderedDict
from pprint import pformat

class DataContainer(OrderedDict):
	"""docstring for DataContainer"""

	# datatypes = {
	# 	'PCS':np.dtype([
	# 				('mdl', int   ),
	# 				('use', bool  ),
	# 				('atm', object),
	# 				('cal', float ),
	# 				('exp', float ),
	# 				('err', float ),
	# 				('idx', int   )]),
	# 	'RDC':np.dtype([
	# 				('mdl', int   ),
	# 				('use', bool  ),
	# 				('atm', object),
	# 				('atn', object),
	# 				('cal', float ),
	# 				('exp', float ),
	# 				('err', float ),
	# 				('idx', int   )]),
	# 	'PRE':np.dtype([
	# 				('mdl', int   ),
	# 				('use', bool  ),
	# 				('atm', object),
	# 				('cal', float ),
	# 				('exp', float ),
	# 				('err', float ),
	# 				('idx', int   )])}

	def __init__(self, *args, **kwargs):
		dtype = kwargs.pop('dtype', None)
		super().__init__(*args, **kwargs)
		self.dtype = dtype

	def __str__(self):
		return pformat(self)

	# def parse(self, protein, models=None, all_atoms=False):
	# 	usedKeys = set([])
	# 	data = []

	# 	if type(models)==int:
	# 		chains = protein[models].get_chains()
	# 	elif type(models) in (list, tuple):
	# 		chains = []
	# 		for m in models:
	# 			chains += protein[m].get_chains()
	# 		print(chains)
	# 	else:
	# 		chains = protein.get_chains()

	# 	if all_atoms:
	# 		if self.dtype=='RDC':
	# 			raise TypeError("Cannot fetch all atoms for dtype RDC")
	# 		atoms = []
	# 		for chain in chains:
	# 			mdl = chain.parent.id
	# 			for atm in chain.get_atoms():
	# 				_, _, _, (_,seq,_), (name, _) = atm.get_full_id()
	# 				key = seq, name
	# 				if key in self:
	# 					exp, err = self[key]
	# 					use = True
	# 					usedKeys.add(key)
	# 				else:
	# 					exp, err = np.nan, np.nan, 
	# 					use = False

	# 				row = mdl, use, atm, np.nan, exp, err, atm.serial_number
	# 				data.append(row)

	# 	elif self.dtype in ('PCS', 'PRE'):
	# 		for chain in chains:
	# 			for key, (exp, err) in self.items():
	# 				seq, name = key
	# 				if seq in chain:
	# 					resi = chain[seq]
	# 					if name in resi:
	# 						mdl = chain.parent.id
	# 						use = True
	# 						atm = resi[name]
	# 						cal = np.nan
	# 						idx = atm.serial_number
	# 						row = mdl, use, atm, cal, exp, err, idx
	# 						data.append(row)
	# 						usedKeys.add(key)

	# 	elif self.dtype=='RDC':
	# 		for chain in chains:
	# 			for key, (exp, err) in self.items():
	# 				(seq1, name1), (seq2, name2) = key
	# 				if seq1 in chain and seq2 in chain:
	# 					resi1 = chain[seq1]
	# 					resi2 = chain[seq2]
	# 					if name1 in resi1 and name2 in resi2:
	# 						mdl = chain.parent.id
	# 						use = True
	# 						atm = resi1[name1]
	# 						atn = resi2[name2]
	# 						cal = np.nan
	# 						idx = unique_pairing(atm.serial_number, 
	# 							atn.serial_number)
	# 						row = mdl, use, atm, atn, cal, exp, err, idx
	# 						data.append(row)
	# 						usedKeys.add(key)

	# 	unused = set(self) - usedKeys
	# 	if unused:
	# 		message = "WARNING: Some values were not parsed to {}:"
	# 		print(message.format(protein.id))
	# 		print(list(unused))
	# 	return np.array(data, dtype=self.datatypes[self.dtype]).view(np.recarray)


def read_pcs(fileName):
	values = DataContainer(dtype='PCS')
	with open(fileName) as o:
		for line in o:
			if line.strip().startswith("#"):
				continue
			seq, name, value, error = line.split()
			key = int(seq), name
			values[key] = float(value), float(error)
	return values


def read_rdc(fileName):
	values = DataContainer(dtype='RDC')
	with open(fileName) as o:
		for line in o:
			if line.strip().startswith("#"):
				continue
			seq1, name1, seq2, name2, value, error = line.split()
			key = frozenset([(int(seq1), name1), (int(seq2), name2)])
			values[key] = float(value), float(error)
	return values


def read_pre(fileName):
	values = DataContainer(dtype='PRE')
	with open(fileName) as o:
		for line in o:
			if line.strip().startswith("#"):
				continue
			seq, name, value, error = line.split()
			key = int(seq), name
			values[key] = float(value), float(error)
	return values


def blank_array(dtype, length):
	return np.recarray(length, DataContainer.datatypes[dtype])

