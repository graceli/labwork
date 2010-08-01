# -*- coding: utf-8 -*-
import result
import xvgfile
import rowtypes	
import os

# This class appends new data into underlying tables
# does not yet implement /dev/shm

class Loader:
	def __init__(self, target):
		print "writing to", target	

		self.target = target
		self._xvgfile = xvgfile.XVGFile()
		self._result = result.Result(target)

	def load(self, analysisName, TableStructure, fixed):
		analysisRoot = os.path.join(self.target, analysisName)

		print "your analysisRoot is", analysisRoot

		dirList = os.listdir(analysisRoot)

		# note I don't like assertions errors
		#assert len(dironly) > 0, "%s does not have any subdirs" % analysisRoot 

		for file in dirList:
			path = os.path.join(analysisRoot, file)
			print "parsing file", path
			data = self._xvgfile.parse(TableStructure, fixed, path)
			print data
			print "adding", len(data), "to table"
			table = self._result.addToTable(data, group='Protein', tableName=analysisName, tableStruct=TableStructure)

# why does it need to return table?
#		return table

if __name__ == "__main__":
	a=Loader()
	table = a.load('rg',rowtypes.RGTable, 4, 'data/rg')

