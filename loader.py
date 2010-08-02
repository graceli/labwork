# -*- coding: utf-8 -*-
import result
import xvgfile
import rowtypes	
import os

# This class appends new data into underlying tables
# does not yet implement /dev/shm

class Loader:
	def __init__(self, location):
		print "working in", location

		self.location = location
		self._xvgfile = xvgfile.XVGFile()
		self._result = result.Result(location)
	

	def load(self, filename, analysisName, TableStructure, fixed):
		analysisRoot = os.path.join(self.location, analysisName)

		print "your analysisRoot is", analysisRoot

		xvgfilepath = os.path.join(analysisRoot, filename)

		print "parsing file", xvgfilepath

		data = self._xvgfile.parse(TableStructure, fixed, xvgfilepath)

		print data
		print "adding", len(data), "to table"

		table = self._result.addToTable(data, group='Protein', tableName=analysisName, tableStruct=TableStructure)

# why does it need to return table?
#		return table

if __name__ == "__main__":
	a=Loader()
	table = a.load('rg',rowtypes.RGTable, 4, 'data/rg')

