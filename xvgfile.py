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
		r=csv.DictReader(open(filename), rowtype.keys(), delimiter=' ', skipinitialspace=True)
		for line in r:
			if line.values()[0] == "#" or line.values()[0] == "@":
			#	print "parsed out", line
				continue
			else:			
				newrow = self._converted(rowtype, line)
				self.data.append(tuple(newrow))
			#	print "appending", newrow
				 
		return self.data
	
	def _converted(self,rowtype, row):
		newrow = []
		for key in rowtype.keys():
		#	print key, rowtype[key]
			if rowtype[key] == rowtypes.Int32Col():
				newrow.append(int(row[key]))
			#	print "converted", key, rowtype[key], "to int"
				
			elif rowtype[key] == rowtypes.Float32Col():
				newrow.append(float(row[key]))
			#	print "converted", key, rowtype[key], "to float"
			else:
				newrow.append(row[key])

		return newrow
		
	def getData(self):
		return self.data



if __name__ == "__main__":
	a=XVGFile()
	data = a.parse(rowtypes.DefaultTable, "part1_replica131.xvg")
	print data
	
