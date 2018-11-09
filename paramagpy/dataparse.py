import numpy as np
from collections import OrderedDict
from pprint import pformat

class DataContainer(OrderedDict):
	"""docstring for DataContainer"""

	def __init__(self, *args, **kwargs):
		dtype = kwargs.pop('dtype', None)
		super().__init__(*args, **kwargs)
		self.dtype = dtype

	def __str__(self):
		return pformat(self)


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

