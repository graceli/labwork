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

	def load(self, analysisName, TableStructure, analysisRoot):
		
		print "your analysisRoot is", analysisRoot

		dirList = os.listdir(analysisRoot)
		#dironly = filter(os.path.isdir, dirList)

		# note I don't like assertions errors
		#assert len(dironly) > 0, "%s does not have any subdirs" % analysisRoot 

		for dir in dirList:
			filesList = os.listdir(os.path.join(analysisRoot,dir))
			#print filesList
			for file in filesList:
				path = os.path.join(analysisRoot, os.path.join(dir,file))
				#print path
				data = self._xvgfile.parse(TableStructure,path)

				print "Data to be added to table:", data

				table = self._result.addToTable(data, group='Protein', tableName=analysisName, tableStruct=TableStructure)
		return table

if __name__ == "__main__":
	a=Loader()
	table = a.load('rg',rowtypes.RGTable, 'data/rg')

