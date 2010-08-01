import rowtypes
import csv
import os

class XVGFile:
	#This class represents the flat xvg file containing data pertaining to a type of analysis
	#(an output of a gromacs analysis tool)
	#The data from the xvg file is represented as an M by N-tuple, where M is the number of rows in the xvg file
	#and N is the number of columns in the row
	def parse(self, rowtype, fixed, filename):
		stripped = os.path.basename(filename)
		parts = stripped.split("_")
		
		self.partnum = parts[0][4:]
		self.replicanum = parts[1][7:]
		
		#read in the contents from the flat xvg file
		self.data = []
		colsDescription = rowtypes.Description(rowtype)
		colNamesInOrder = colsDescription._v_names[fixed:]	

		#print "colNamesInOrder:", colNamesInOrder
	
		#r=csv.DictReader(open(filename), colNamesInOrder, delimiter=' ', skipinitialspace=True)
		r = csv.reader(open(filename), delimiter=' ', skipinitialspace=True)
		#r = open(filename, 'r')

		numappended = 0
		for line in r:
			#print "line=", line
			if self._find('#', line) or self._find('@', line):
				#print "parsed out", line
				continue
			else:			
				#print "inserting", line
				#reconstructo dictionary for row
				rowdict = {}
				for i in range(0, len(colNamesInOrder)):
					key = colNamesInOrder[i]
					rowdict[key] = line[i]
  				
				rowdict['temp'] = 0 	
				rowdict['replicanum'] = int(self.replicanum)
				rowdict['partnum'] = int(self.partnum)
				rowdict['seqnum'] = 8

				newrow = [ rowdict['temp'], rowdict['replicanum'], rowdict['seqnum'], rowdict['partnum'] ]
				#print "inserting line as dict", rowdict	
				readinrow = self._converted(colsDescription._v_types, colNamesInOrder, rowdict)
				newrow.extend(readinrow)	
				#print newrow
				#print "appending new row (as a tuple)", newrow
				numappended+=1
				self.data.append(tuple(newrow))

		print "total number of lines appended=", numappended 

		return self.data

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
		#print "initial newrow", newrow
		for key in vtypes.keys():
			# find position of the current key (table column name)
			# in the in order list of column names
			pos = -1 
			if key in colnames:
				pos = colnames.index(key) 
			
			if pos > -1:	
				#print key, pos, line[key]
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

	data = a.parse(rowtypes.RGTable, 4, "data/rg/part1/part1_replica131_rg.xvg")
	#print data
	
