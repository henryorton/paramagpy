import numpy as np
from pprint import pprint


class DataContainer(dict):
	"""docstring for DataContainer"""
	def __init__(self, *args, **kwargs):
		dtype = kwargs.pop('dtype', None)
		super().__init__(*args, **kwargs)
		self.dtype = dtype

def read_pcs(fileName):
	npc_values = DataContainer(dtype='pcs')
	with open(fileName) as o:
		for line in o:
			seq, name, value, error = line.split()
			key = int(seq), name
			npc_values[key] = float(value)

	return npc_values

def read_rdc(fileName):
	pass




