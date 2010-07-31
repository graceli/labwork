import rowtypes
import csv
import os

class XVGFile:
	#This class represents the flat xvg file containing data pertaining to a type of analysis
	#(an output of a gromacs analysis tool)
	#The data from the xvg file is represented as an M by N-tuple, where M is the number of rows in the xvg file
	#and N is the number of columns in the row
	def parse(self, rowtype, filename):
		stripped = os.path.basename(filename)
		parts = stripped.split("_")
		
		self.partnum = parts[0][4:]
		self.replicanum = parts[1][7:]
		
		#read in the contents from the flat xvg file
		self.data = []
		colsDescription = rowtypes.Description(rowtype)
		colNamesInOrder = colsDescription._v_names	

		print "colNamesInOrder:", colNamesInOrder
	
		r=csv.DictReader(open(filename), colNamesInOrder, delimiter=' ', skipinitialspace=True)

		for line in r:
			print line
			if "#" in line.values() or "@" in line.values():
				print "parsed out", line
				continue
			else:			
				print "inserting", line
				newrow = self._converted(colsDescription._v_types, colNamesInOrder, line)
				self.data.append(tuple(newrow))
				print "appending", newrow
				 
		return self.data
	
	def _converted(self, vtypes, colnames, line):
		# pre: colnames are assumed to be in order
		#      vtypes is a dictionary of table column names and their 
		#	   data types as defined in PyTable
		# post: newrow is a table row with table columns ordered 
		#       exactly as provided in colnames

		newrow = [0]*len(colnames)
		print "initial newrow", newrow
		for key in vtypes.keys():
			# find position of the current key (table column name)
			# in the in order list of column names
			pos = colnames.index(key) 
			print key, pos, line[key]
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



if __name__ == "__main__":
	a=XVGFile()
	#data = a.parse(rowtypes.DefaultTable, "data/default/part1/part1_replica131.xvg")
	#print data

	data = a.parse(rowtypes.RGTable, "data/rg/part1/part1_replica131.xvg")
	print data
	
