import rowtypes
import csv
import os

# Strategy pattern to implement parsing of multiple file formats
class FileFormat(object):
	"""docstring for FileFormat"""
	def __init__(self, arg):
		super(FileFormat, self).__init__()
		self.arg = arg
	
	def parse(self):
		"""docstring for parse"""
		pass

	def get_as_matrix(self):
		"""docstring for get_as_matrix"""
		pass

	
class XVGFile:	
	"""This class represents the flat xvg file containing data pertaining to a type of analysis (an output of a gromacs analysis tool). The data from the xvg file is represented as an M by N-tuple, where M is the number of rows in the xvg file and N is the number of columns in the row"""
	
	def __init__(self):
		self.data = []

	def __del__(self):
		del self.data

	def parse(self, rowtype, fixed, filename):
		# parse out information from file name for the 'fixed' columns
		# refactored meta data extraction from filename into indexing system

		# read in the contents from the flat xvg file
		data = []
		colsDescription = rowtypes.Description(rowtype)
		colNamesInOrder = colsDescription._v_names[fixed:]	
		r = csv.reader(open(filename), delimiter=' ', skipinitialspace=True)
		numappended = 0

		print "XVGFile.parse():"
		print colNamesInOrder
		print fixed
		print filename
		
		for line in r:
			if self._find('#', line) or self._find('@', line) or not line:
				continue
			else:			
				#reconstruct dictionary for row
				rowdict = {}
				for i in range(0, len(colNamesInOrder)):
					key = colNamesInOrder[i]
					rowdict[key] = line[i]
  				
				newrow = self._converted(colsDescription._v_types, colNamesInOrder, rowdict)

				numappended+=1
				data.append(tuple(newrow))

		return data

	def _find(self, str, list):
		for item in list:
			if item != None:
				if item.find(str) > -1:
					return True

		return False
	
	def _converted(self, vtypes, colnames, line):
		# pre: colnames are assumed to be in order
		#      vtypes is a dictionary of table column names and their 
		#	   data types as defined in PyTable
		# post: newrow is a table row with table columns ordered 
		#       exactly as provided in colnames

		newrow = [0]*len(colnames)
		for key in vtypes.keys():
			# find position of the current key (table column name)
			# in the in order list of column names
			pos = -1 
			if key in colnames:
				pos = colnames.index(key) 
			
			if pos > -1:	
				if vtypes[key] == 'int32':
					newrow[pos] = int(line[key])
				elif vtypes[key]== 'float32':
					newrow[pos] = float(line[key])
				else:
					# note inserting value of unknown type
					# might cause an error
					newrow[pos] = line[key]

		return newrow
		
	def getData(self):
		return self.data

