# -*- coding: utf-8 -*-
import result
import xvgfile
import rowtypes	
import os

# This class appends new data into underlying tables

class Loader:
	def __init__(self):
		self._xvgfile = xvgfile.XVGFile()
		self._result = result.Result()

	def load(self, analysisName, TableStructure, fixed, analysisRoot):
		
		print "your analysisRoot is", analysisRoot

		dirList = os.listdir(analysisRoot)
		#dironly = filter(os.path.isdir, dirList)

		# note I don't like assertions errors
		#assert len(dironly) > 0, "%s does not have any subdirs" % analysisRoot 

		for dir in dirList:
			filesList = os.listdir(os.path.join(analysisRoot,dir))
			for file in filesList:
				path = os.path.join(analysisRoot, os.path.join(dir,file))
				print "parsing file", path
				data = self._xvgfile.parse(TableStructure, fixed, path)
				#print data
				print "adding", len(data), "to table"
				table = self._result.addToTable(data, group='Protein', tableName=analysisName, tableStruct=TableStructure)
		return table

if __name__ == "__main__":
	a=Loader()
	table = a.load('rg',rowtypes.RGTable, 4, 'data/rg')

