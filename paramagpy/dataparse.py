import numpy as np
from collections import OrderedDict
from pprint import pformat

class DataContainer(OrderedDict):
	"""
	A dictionary-like container for storing 
	PCS, RDC, PRE and CCR data
	Has an additional attribute 'dtype' to define datatype
	"""
	def __init__(self, *args, **kwargs):
		dtype = kwargs.pop('dtype', None)
		super().__init__(*args, **kwargs)
		self.dtype = dtype

	def __str__(self):
		return pformat(self)


def read_pcs(fileName):
	"""
	Read pseudo contact shift values from file.
	The returned object is a dicationary.
	They keys are tuples of (sequence, atomName)
	The values are tuples of (value, error)

	Parameters
	----------
	fileName : str
		the path to the file

	Returns
	-------
	values : :class:`paramagpy.dataparse.DataContainer`
		a dictionary containing the parsed data

	Examples
	--------
	>>> values = paramagpy.dataparse.read_pcs("calbindin_Er_HN_PCS_errors.npc")
	>>> for v in values.items():
	...     print(v)
	... 
	((2, 'H'), (-0.04855485, 0.0016))
	((2, 'N'), (-0.03402764, 0.0009))
	((4, 'H'), (-0.18470315, 0.0004))
	...
	((75, 'H'), (0.19553661, 0.0005))
	((75, 'N'), (0.17840666, 0.0004))
	"""
	values = DataContainer(dtype='PCS')
	with open(fileName) as o:
		for line in o:
			try:
				if line.strip().startswith("#"):
					continue
				seq, name, value, error = line.split()
				key = int(seq), name
				values[key] = float(value), float(error)
			except ValueError:
				print("Line ignored while reading file: {}\n{}".format(
					fileName, line))
	return values


def read_rdc(fileName):
	"""
	Read residual dipolar coupling values from file.
	The returned object is a dicationary.
	They keys are frozensets of tuples of the form:
	frozenset({(sequence1, atomName1), (sequence2, atomName2)})
	The frozenset only allows unordered unique atom identification pairs
	The values are tuples of (value, error)

	Parameters
	----------
	fileName : str
		the path to the file

	Returns
	-------
	values : :class:`paramagpy.dataparse.DataContainer`
		a dictionary containing the parsed data

	Examples
	--------
	>>> values = paramagpy.dataparse.read_rdc("ubiquitin_a28c_c1_Tb_HN.rdc")
	>>> for v in values.items():
	...     print(v)
	... 
	(frozenset({(2, 'N'), (2, 'H')}), (-2.35, 0.32))
	(frozenset({(3, 'N'), (3, 'H')}), (-4.05, 0.38))
	(frozenset({(4, 'H'), (4, 'N')}), (-3.58, 0.42))
	...
	(frozenset({(73, 'N'), (73, 'H')}), (-0.47, 0.75))
	(frozenset({(76, 'H'), (76, 'N')}), (0.14, 0.3))
	"""
	values = DataContainer(dtype='RDC')
	with open(fileName) as o:
		for line in o:
			try:
				if line.strip().startswith("#"):
					continue
				seq1, name1, seq2, name2, value, error = line.split()
				key = frozenset([(int(seq1), name1), (int(seq2), name2)])
				values[key] = float(value), float(error)
			except ValueError:
				print("Line ignored while reading file: {}\n{}".format(
					fileName, line))
	return values


def read_pre(fileName):
	"""
	Read paramagnetic relaxation enhancement values from file.
	The returned object is a dicationary.
	They keys are tuples of (sequence, atomName)
	The values are tuples of (value, error)

	Parameters
	----------
	fileName : str
		the path to the file

	Returns
	-------
	values : :class:`paramagpy.dataparse.DataContainer`
		a dictionary containing the parsed data

	Examples
	--------
	see :func:`paramagpy.dataparse.read_pcs` which has the same file structure
	"""
	values = DataContainer(dtype='PRE')
	with open(fileName) as o:
		for line in o:
			try:
				if line.strip().startswith("#"):
					continue
				seq, name, value, error = line.split()
				key = int(seq), name
				values[key] = float(value), float(error)
			except ValueError:
				print("Line ignored while reading file: {}\n{}".format(
					fileName, line))
	return values


def read_ccr(fileName):
	"""
	Read cross-correlated relaxation values from file.
	These are typically Curie-spin cross Dipole-dipole relaxation rates
	The returned object is a dicationary.
	They keys are tuples of the form:
	((sequence1, atomName1), (sequence2, atomName2))
	Note that the first column is for the active nucleus undergoing 
	relaxation and the second column is for the partner spin.
	The values are tuples of (value, error)

	Parameters
	----------
	fileName : str
		the path to the file

	Returns
	-------
	values : :class:`paramagpy.dataparse.DataContainer`
		a dictionary containing the parsed data

	Examples
	--------
	see :func:`paramagpy.dataparse.read_rdc` which has the similar file structure
	"""
	values = DataContainer(dtype='CCR')
	with open(fileName) as o:
		check1 = set([])
		check2 = set([])
		for line in o:
			try:
				if line.strip().startswith("#"):
					continue
				seq1, name1, seq2, name2, value, error = line.split()
				key = (int(seq1), name1), (int(seq2), name2)
				values[key] = float(value), float(error)
				check1.add(name1)
				check2.add(name2)
			except ValueError:
				print("Line ignored while reading file: {}\n{}".format(
					fileName, line))
		if len(check1 & check2) > 0:
			print("WARNING: Varied atom ordering detected in file")
	return values


