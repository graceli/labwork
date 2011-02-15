import numpy
import tables

def convert_to_numpy(table, dtype=numpy.float32):
	numpy_array = table.read().view(dtype=dtype).reshape(-1, len(table[0]))
	return numpy_array