import rowtypes
import csv
import os

class XVGFile:
	#This class represents the flat xvg file containing data pertaining to a type of analysis
	#(an output of a gromacs analysis tool)
	#The data from the xvg file is represented as an M by N-tuple, where M is the number of rows in the xvg file
	#and N is the number of columns in the row
	def __init__(self):
		self.data = []

	def __del__(self):
		del self.data

	def parse(self, rowtype, fixed, filename):
		#parse out information from file name for the 'fixed' columns
		stripped = os.path.basename(filename)
		self.replicanum, self.seqnum, self.temp,garb = stripped[2:].split(".")
		
		#read in the contents from the flat xvg file
		data = []

		colsDescription = rowtypes.Description(rowtype)
		colNamesInOrder = colsDescription._v_names[fixed:]	

		print colNamesInOrder

		#r=csv.DictReader(open(filename), colNamesInOrder, delimiter=' ', skipinitialspace=True)
		r = csv.reader(open(filename), delimiter=' ', skipinitialspace=True)

		numappended = 0
		for line in r:
			if self._find('#', line) or self._find('@', line):
				continue
			else:			
				print line
				#print "inserting", line
				#reconstruct dictionary for row
				rowdict = {}
				for i in range(0, len(colNamesInOrder)):
					print i
					key = colNamesInOrder[i]
					rowdict[key] = line[i]
  				
				rowdict['temp'] = int(self.temp)
				rowdict['replicanum'] = int(self.replicanum)
#				rowdict['partnum'] = int(self.partnum)
				rowdict['seqnum'] = int(self.seqnum)
				
				newrow = [ rowdict['temp'], rowdict['replicanum'], rowdict['seqnum'] ]
				readinrow = self._converted(colsDescription._v_types, colNamesInOrder, rowdict)
				newrow.extend(readinrow)	
				print "appending new row (as a tuple)", newrow
				print "size",len(newrow)
				numappended+=1
				data.append(tuple(newrow))
				print "data so far is", data


		print "total number of lines appended=", numappended 

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

