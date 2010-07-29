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
		dirList = os.listdir(analysisRoot)
		print dirList	
		for dir in dirList:
			filesList = os.listdir(os.path.join(analysisRoot,dir))
			print filesList
			for file in filesList:
				path = os.path.join(analysisRoot, os.path.join(dir,file))
				print path
				data = self._xvgfile.parse(TableStructure,path)
				table = self._result.addToTable(data, group='Protein', tableName=analysisName, tableStruct=TableStructure)
		return table

if __name__ == "__main__":
	a=Loader()
	table = a.load('default',rowtypes.DefaultTable, 'default')

