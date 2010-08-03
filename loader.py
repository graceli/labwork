# -*- coding: utf-8 -*-
import result
import xvgfile
import rowtypes	
import os

# This class appends new data into underlying tables
# does not yet implement /dev/shm

class Loader:
	def __init__(self, location):
		# location is the execution path of this class
		self.location = location
		self._xvgfile = xvgfile.XVGFile()
		self._result = result.Result(location)
	

	def load(self, filename, analysisName, TableStructure, fixed):
		# filename is the name of the analysis file
		# analysisName is the name of the analysis eg. 'rg'
		# TableStructure is the rowtype object to be created and inserted into the table
		# fixed number of fixed data columns in table (ie. info is not read in from 'filename'
		
		if analysisName.find(".") != -1:
			name,ext=analysisName.split(".")
		else:
			name = analysisName

		analysisRoot = os.path.join(self.location, name)
		xvgfilepath = os.path.join(analysisRoot, filename)
		
		#print "Loader.load: parsing", xvgfilepath
		
		data = self._xvgfile.parse(TableStructure, fixed, xvgfilepath)
		table = self._result.addToTable(data, group='Protein', tableName=analysisName, tableStruct=TableStructure)


#if __name__ == "__main__":
#	a=Loader()
#	table = a.load('rg',rowtypes.RGTable, 4, 'data/rg')

